#!/usr/bin/env python
# coding: utf-8

# # Molecular Dynamics Analysis
# 

# **PBC**

# In[1]:


get_ipython().system('echo "4 0" | gmx trjconv -s step5_1.tpr -f step5_1.xtc -n index_new.ndx -pbc mol -center -ur compact -o traj_pbc.xtc')


# **fit rot+trans**

# In[2]:


get_ipython().system('echo "4 0" | gmx trjconv -s step5_1.tpr -f traj_pbc.xtc -n index_new.ndx -fit rot+trans -o traj_fit.xtc')


# In[1]:


import subprocess

# '1 | 14' означает объединить (OR) группу 1 (Protein) и 14 (NAD)
ndx_commands = """
1 | 14
r 170
r 263
r 109
r 110
18 | 19
16 | 17
r 14 15 22 86 90 93 97 108 109 110 111 115 118 145 146 147 148 149 168 169 170 171 172 173 174 175 177 178 195 199 200 201 202 204 228 238 242 243 244 245 246 247 248 260 261 262 263 264 265 266 273 276 277 281 283 284 286 287 291 314 315 316 317 318 319 321 357 445
name 36 pocket
r 93
27 | 29
name 38 complex
q
"""

# 2. Запускаем gmx make_ndx
process = subprocess.run(
    ["gmx", "make_ndx", "-f", "step5_1.gro", "-n", "index_new.ndx", "-o", "index_updated.ndx"],
    input=ndx_commands,      # Передаем наши команды
    text=True,               # Говорим, что работаем со строками, а не байтами
    capture_output=True      # Перехватываем вывод, чтобы он не засорял фоновый терминал
)


# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns # Для красивой матрицы

# Настраиваем стиль графиков для статьи
plt.rcParams.update({'font.size': 12, 'font.family': 'sans-serif'})


# ### Основные функции для анализа

# In[3]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def load_gro_atom_mapping(gro_filename):
    """Парсит .gro файл и возвращает словарь {номер_атома: 'НомерОстатка_ИмяОстатка@ИмяАтома'}"""
    atom_map = {}
    print(f"Парсим структуру из {gro_filename}...")
    with open(gro_filename, 'r') as f:
        lines = f.readlines()
    # Пропускаем 2 строки заголовка и 1 строку вектора ячейки в конце
    for line in lines[2:-1]:
        if len(line) < 20: continue
        try:
            res_num   = line[0:5].strip()
            res_name  = line[5:10].strip()
            atom_name = line[10:15].strip()
            atom_num  = int(line[15:20].strip())
            atom_map[atom_num] = f"{res_num}{res_name}@{atom_name}"
        except ValueError:
            continue
    print(f"Успешно загружено атомов: {len(atom_map)}")
    return atom_map

def load_hbond_pairs_from_ndx(ndx_filename):
    """Парсит hbond.ndx и возвращает словарь {индекс_связи: (атом_водорода, атом_акцептора)}."""
    pairs = {}
    current_idx = 0
    in_target_section = False
    with open(ndx_filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('['):
                # Ищем секции с водородными связями
                if 'hbonds' in line.lower(): 
                    in_target_section = True
                else: 
                    in_target_section = False
                continue

            if in_target_section:
                parts = line.split()
                # GROMACS выводит: Донор(0) Водород(1) Акцептор(2)
                # Берем мостик H-A
                if len(parts) >= 3:
                    pairs[current_idx] = (int(parts[1]), int(parts[2]))
                    current_idx += 1
    return pairs

def parse_xpm(filename):
    """Парсер для конвертации XPM в Numpy Array."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    data_lines = [line.split('"')[1] for line in lines if line.startswith('"') and len(line.split('"')) > 1]
    meta = data_lines[0].split()
    y_dim = int(meta[1])
    num_colors = int(meta[2])
    matrix_lines = data_lines[1 + num_colors :]
    binary_matrix = np.zeros((y_dim, len(matrix_lines[0])), dtype=int)
    for i, line in enumerate(matrix_lines):
        for j, char in enumerate(line):
            if char != ' ':
                binary_matrix[i, j] = 1
    return binary_matrix


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def analyze_hbond_matrix(gro_file="step5_1.gro", ndx_file="hbond.ndx", xpm_file="hbmat.xpm", cutoff=20.0):
    """
    Загружает результаты gmx hbond, фильтрует по occupancy, 
    рисует отсортированную хитмапу и выводит текстовый отчет.

    Параметры:
    - gro_file: путь к файлу координат (.gro)
    - ndx_file: путь к индексному файлу (.ndx)
    - xpm_file: путь к матрице существования связей (.xpm)
    - cutoff: порог отсечения по времени существования (в процентах)

    Возвращает:
    - df_report: pandas.DataFrame с результатами
    """

    # 1. Загружаем все файлы (предполагается, что функции парсинга уже определены в ячейке выше)
    gro_atoms = load_gro_atom_mapping(gro_file)
    atom_mapping = load_hbond_pairs_from_ndx(ndx_file)
    hb_matrix_full = parse_xpm(xpm_file)

    # 2. Считаем базовый occupancy
    occupancy = np.mean(hb_matrix_full, axis=1) * 100
    print(f"Данные из {xpm_file} успешно загружены. Фильтрация по порогу >{cutoff}%...")

    valid_indices = [i for i, occ in enumerate(occupancy) if occ > cutoff]

    if len(valid_indices) == 0:
        print(f"⚠️ Ни одна связь не превысила порог в {cutoff}%.")
        return pd.DataFrame()

    # Обрезаем матрицу под выбранный порог
    hb_matrix_filtered = hb_matrix_full[valid_indices, :]

    report_data = []
    yticklabels_labels = [] 

    for idx, i in enumerate(valid_indices):
        occ = occupancy[i]
        h_atom_idx, a_atom_idx = atom_mapping.get(i, (None, None))

        # Вытаскиваем красивые лейблы из .gro
        h_label = gro_atoms.get(h_atom_idx, f"Atom_{h_atom_idx}")
        a_label = gro_atoms.get(a_atom_idx, f"Atom_{a_atom_idx}")

        report_data.append({
            "bond_idx": i,
            "occupancy": occ,
            "hydrogen_label": h_label,
            "acceptor_label": a_label
        })
        yticklabels_labels.append(f"Idx {i} ({h_label} ... {a_label})")

    df_report = pd.DataFrame(report_data)

    # 3. Синхронная сортировка: чтобы график и таблица шли от самых частых связей к редким
    sort_order = np.argsort(df_report["occupancy"])[::-1]
    df_report = df_report.iloc[sort_order].reset_index(drop=True)
    hb_matrix_filtered = hb_matrix_filtered[sort_order, :]
    yticklabels_labels = [yticklabels_labels[j] for j in sort_order]

    # 4. Визуализация хитмапы
    plt.figure(figsize=(14, 2 + len(valid_indices) * 0.4)) 
    ax = plt.gca()
    sns.heatmap(hb_matrix_filtered, cmap="Blues", cbar=False, yticklabels=yticklabels_labels, ax=ax)
    plt.title(f"H-Bond Existence Map (Occupancy > {cutoff}%)")
    plt.xlabel("Frame")
    plt.ylabel("Hydrogen Bond (H ... A)")
    plt.yticks(rotation=0) 
    plt.tight_layout()
    plt.show()

    # 5. Текстовый отчет
    print("\n" + "="*75)
    print(f"{f'FINAL HYDROGEN BONDS REPORT (Occupancy > {cutoff}%)':^75}")
    print("="*75)
    print(f"{'Index':<6} | {'Occupancy':<9} | {'Hydrogen Group':<20} <...> {'Acceptor Group':<20}")
    print("-"*75)
    for _, row in df_report.iterrows():
        print(f"{int(row['bond_idx']):<6} | {row['occupancy']:>7.1f}% | {row['hydrogen_label']:<20} <...> {row['acceptor_label']:<20}")
    print("="*75)

    return df_report



import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
from IPython.display import display # На всякий случай импортируем display для Jupyter

def analyze_protein_ligand_hbonds(u, protein_resids, ligand_resnames, d_a_cutoff=3.5, d_h_a_angle_cutoff=150):
    """
    Анализирует водородные связи между заданным массивом остатков белка и списком лигандов.

    Параметры:
    - u: MDAnalysis.Universe объект
    - protein_resids: список номеров остатков белка (например, [93, 109, 110])
    - ligand_resnames: список имен лигандов (например, ["4GA", "0GA", "ROH"])
    - d_a_cutoff: порог расстояния Донор-Акцептор (по умолчанию 3.5 Å)
    - d_h_a_angle_cutoff: порог угла D-H-A (по умолчанию 150 градусов)

    Возвращает:
    - df_hb: Pandas DataFrame с результатами
    """

    print(f"=== Старт анализа водородных связей ===")
    print(f"Остатки белка: {protein_resids}")
    print(f"Лиганды: {ligand_resnames}\n")

    # Автоматически собираем синтаксические маски для MDAnalysis
    protein_query = " or ".join([f"resid {r}" for r in protein_resids])
    ligand_query  = " or ".join([f"resname {l}" for l in ligand_resnames])

    # Строим строгие маски для атомов (избегаем element, ищем по name)
    hydrogens_prot = f"({protein_query}) and name H*"
    acceptors_lig  = f"({ligand_query}) and (name O* or name N*)"

    hydrogens_lig  = f"({ligand_query}) and name H*"
    acceptors_prot = f"({protein_query}) and (name O* or name N*)"

    donors_prot = f"({protein_query}) and (name O* or name N*)"
    donors_lig  = f"({ligand_query}) and (name O* or name N*)"


    print(f"1. Запуск Раунда А (Белок ---> Лиганд)...")
    hb_round_A = HydrogenBondAnalysis(
        universe=u,
        donors_sel=donors_prot, 
        hydrogens_sel=hydrogens_prot,
        acceptors_sel=acceptors_lig,
        d_a_cutoff=d_a_cutoff,
        d_h_a_angle_cutoff=d_h_a_angle_cutoff,
        update_selections=False
    )
    hb_round_A.run()

    print(f"2. Запуск Раунда Б (Лиганд ---> Белок)...")
    hb_round_B = HydrogenBondAnalysis(
        universe=u,
        donors_sel=donors_lig,
        hydrogens_sel=hydrogens_lig, 
        acceptors_sel=acceptors_prot,
        d_a_cutoff=d_a_cutoff,
        d_h_a_angle_cutoff=d_h_a_angle_cutoff,
        update_selections=False
    )
    hb_round_B.run()

    # 3. Объединяем и обрабатываем результаты
    print("3. Сбор и сортировка результатов...")
    total_frames = len(u.trajectory)
    bond_counts = {}

    for hb_instance in [hb_round_A, hb_round_B]:
        if hb_instance.results.hbonds is not None and len(hb_instance.results.hbonds) > 0:
            for bond in hb_instance.results.hbonds:
                frame, d_idx, h_idx, a_idx, dist, ang = bond
                pair = (int(d_idx), int(a_idx))
                if pair not in bond_counts:
                    bond_counts[pair] = 0
                bond_counts[pair] += 1

    # Переводим индексы атомов в красивые названия
    report = []
    for (d_idx, a_idx), freq in bond_counts.items():
        donor = u.atoms[d_idx]
        acceptor = u.atoms[a_idx]
        occ = (freq / total_frames) * 100

        report.append({
            'Donor': f"{donor.resid}{donor.resname}@{donor.name}",
            'Acceptor': f"{acceptor.resid}{acceptor.resname}@{acceptor.name}",
            'Frequency': freq,
            'Occupancy (%)': occ
        })

    df_hb = pd.DataFrame(report)

    # 4. Вывод финального отчета
    if not df_hb.empty:
        df_hb = df_hb.sort_values(by='Occupancy (%)', ascending=False).reset_index(drop=True)
        print("\n" + "="*85)
        print(f"{'ARRAY HYDROGEN BONDS REPORT':^85}")
        print("="*85)
        print(f"{'Occupancy':<9} | {'Donor (Отдает H)':<30} ---> {'Acceptor (Принимает H)':<30}")
        print("-"*85)
        for _, row in df_hb.iterrows():
            print(f"{row['Occupancy (%)']:>7.1f}% | {row['Donor']:<30} ---> {row['Acceptor']:<30}")
        print("="*85 + "\n")

        # Выводим топ-20 в красивом формате Jupyter
        display(df_hb.head(20))
    else:
        print(f"\n❌ Между заданным массивом остатков белка и лигандами водородных связей не найдено.")

    return df_hb


import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm.auto import tqdm

# ====================================================================
# 1. ОПРЕДЕЛЕНИЕ ФУНКЦИЙ АНАЛИЗА И ВИЗУАЛИЗАЦИИ
# ====================================================================

def parse_gro_mapping(gro_path="step5_1.gro"):
    """
    Парсит .gro файл по фиксированной ширине колонок 
    и возвращает словарь { 'ResidResName@AtomName': 0_based_index }
    """
    gro_mapping = {}
    atom_idx_0_based = 0

    with open(gro_path, "r") as f:
        lines = f.readlines()
        # Пропускаем заголовок и последнюю строку с боксом
        for line in lines[2:-1]:
            if len(line.strip()) < 20:
                continue
            res_num = line[0:5].strip()
            res_name = line[5:10].strip()
            atom_name = line[10:15].strip()

            key = f"{res_num}{res_name}@{atom_name}"
            gro_mapping[key] = atom_idx_0_based
            atom_idx_0_based += 1

    print(f"   [GRO] Успешно проиндексировано {len(gro_mapping)} атомов из {gro_path}")
    return gro_mapping


def calculate_distances_by_indices(u, gro_mapping, bonds_list):
    """
    Проходит по траектории MDAnalysis ОДИН раз и считает расстояния 
    для всех пар, найденных в маппинге.
    """
    # Фильтруем список и собираем только существующие в .gro пары
    valid_pairs = []
    for d_label, a_label in bonds_list:
        if d_label not in gro_mapping or a_label not in gro_mapping:
            print(f"⚠️ Предупреждение: Пара {d_label} -> {a_label} отсутствует в .gro! Пропускаем.")
            continue
        valid_pairs.append((gro_mapping[d_label], gro_mapping[a_label], f"{d_label} ... {a_label}"))

    if not valid_pairs:
        print("❌ Нет валидных пар для расчета!")
        return None, None

    # Инициализация структур данных
    distances_dict = {label: [] for _, _, label in valid_pairs}
    times = []

    print(f"   [MD] Расчет расстояний для {len(valid_pairs)} связей по траектории...")
    for ts in tqdm(u.trajectory, desc="Траектория"):
        times.append(ts.time)
        for idx1, idx2, label in valid_pairs:
            p1 = u.atoms[idx1].position
            p2 = u.atoms[idx2].position
            distances_dict[label].append(np.linalg.norm(p1 - p2))

    # Переводим время в наносекунды
    times = np.array(times) / 1000.0
    return times, distances_dict


def plot_individual_bonds(times, distances_dict, title_prefix="Cluster"):
    """
    Строит индивидуальный отдельный график для каждой связи в виде сетки subplots.
    Включает сырые данные, сглаженную линию и линию отсечки 3.5 Å.
    """
    labels = list(distances_dict.keys())
    n_plots = len(labels)

    # Задаем сетку: 1 колонка, а строк столько, сколько у нас связей
    fig, axes = plt.subplots(n_plots, 1, figsize=(14, 3.5 * n_plots), sharex=True)

    # Если связь всего одна, axes не будет массивом, превращаем его в список вручную
    if n_plots == 1:
        axes = [axes]

    fig.suptitle(f"Individual Hydrogen Bond Dynamics: {title_prefix}", fontsize=16, fontweight='bold', y=1.01)

    for i, label in enumerate(labels):
        ax = axes[i]
        raw_dist = np.array(distances_dict[label])

        # Сглаживание скользящим средним (окно 10 фреймов)
        smoothed_dist = pd.Series(raw_dist).rolling(window=10, min_periods=1).mean()

        # Рисуем сырые полупрозрачные данные на фоне тонкой линией
        ax.plot(times, raw_dist, color='gray', alpha=0.25, linewidth=0.8, label="Raw data")
        # Рисуем красивую сглаженную линию основного тренда
        ax.plot(times, smoothed_dist, color='tab:blue', alpha=0.9, linewidth=1.5, label="Smoothed trend")

        # Наносим красный пунктир границы H-связи (3.5 Ангстрема)
        ax.axhline(y=3.5, color='tab:red', linestyle='--', linewidth=1.5, label="Cutoff (3.5 Å)")

        # Оформление каждого отдельного сабплота
        ax.set_title(label, fontsize=11, fontweight='semibold', loc='left', color='dimgray')
        ax.set_ylabel("Distance (Å)", fontsize=10)
        ax.set_ylim(2.0, 5.5)
        ax.grid(True, linestyle=':', alpha=0.5)
        if i == 0:
            ax.legend(loc="upper right", fontsize=9, frameon=True)

    # Добавляем подпись оси X только самому нижнему графику
    axes[-1].set_xlabel("Time (ns)", fontsize=12)

    plt.tight_layout()
    plt.show()


def df_to_bond_list(df, min_occupancy=0.0):
    """
    Преобразует DataFrame с результатами водородных связей в список кортежей (Donor, Acceptor).

    Параметры:
    - df: Pandas DataFrame (результат функции analyze_protein_ligand_hbonds)
    - min_occupancy: минимальный процент Occupancy. Связи ниже этого порога не попадут в список.

    Возвращает:
    - bond_list: Список кортежей, готовый для графиков.
    """
    if df is None or df.empty:
        return []

    # Отфильтровываем слабые связи (по умолчанию min_occupancy=0.0, то есть берем все)
    filtered_df = df[df['Occupancy (%)'] >= min_occupancy]

    # Склеиваем две колонки в список кортежей с помощью zip
    bond_list = list(zip(filtered_df['Donor'], filtered_df['Acceptor']))

    return bond_list




# In[6]:


get_ipython().system('echo "25 36" |gmx hbond -f 4_traj_fit.xtc -s step5_1.tpr -n index_updated.ndx -num hbnum.xvg -dist hbdist.xvg -hbm hbmat.xpm -hbn hbond.ndx -nomerge')


# In[7]:


# Базовый запуск (с дефолтными именами файлов и порогом 20%)
results_df = analyze_hbond_matrix()

# Запуск для других файлов или с другим порогом
# results_df = analyze_hbond_matrix(gro_file="step5_2.gro", xpm_file="hbmat_replica2.xpm", cutoff=40.0)


# In[8]:


import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

# File paths
tpr_file = "step5_1.tpr"
xtc_file = "4_traj_fit.xtc"
ndx_file = "index_updated.ndx"

# Load Universe
u = mda.Universe(tpr_file, xtc_file)
print(f"Loaded Universe with {len(u.trajectory)} frames")


# ## 1. RMSD Analysis
# 
# We calculate RMSD for the Protein and index groups 26, 27, and 28, aligned to the Protein backbone.

# In[9]:


# Manual parser for GROMACS .ndx files
def parse_ndx(filename):
    groups = {}
    current_group = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('['):
                current_group = line.strip('[] ').strip()
                groups[current_group] = []
            elif current_group and line:
                groups[current_group].extend([int(x) - 1 for x in line.split()])
    return groups

ndx_groups = parse_ndx(ndx_file)

# Selections
selections = {
    "Protein": "protein",
    "Group 26 (non-Water)": "index " + " ".join(map(str, ndx_groups["non-Water"])),
    "Group 27 (ROH_4GA_0GA_PHO)": "index " + " ".join(map(str, ndx_groups["M6P"])),
    "Group 28 (Complex)": "index " + " ".join(map(str, ndx_groups["Protein_NAD"]))
}

rmsd_results = {}

for name, sel in selections.items():
    print(f"Calculating RMSD for {name}...")
    R = rms.RMSD(u, u, select=sel, superposition_selection="protein", ref_frame=0)
    R.run()
    rmsd_results[name] = R.results.rmsd


# ## 2. Bond Distance Analysis
# 
# We calculate distances over time for the following pairs:
# - ASH 170 (OD2) -- 4GA 448 (O4)
# - TYD 263 (OH) -- 0GA 449 (H2)

# In[10]:


# Atom indices (0-based)
bond1_idx = (2648, 7039) # ASH 170 OD2 - 4GA 448 O4
bond2_idx = (4108, 7043) # TYD 263 OH - 0GA 449 H2

b1_dist = []
b2_dist = []
times = []

print("Calculating bond distances...")
for ts in tqdm(u.trajectory):
    p1 = u.atoms[bond1_idx[0]].position
    p2 = u.atoms[bond1_idx[1]].position
    p3 = u.atoms[bond2_idx[0]].position
    p4 = u.atoms[bond2_idx[1]].position

    b1_dist.append(np.linalg.norm(p1 - p2))
    b2_dist.append(np.linalg.norm(p3 - p4))
    times.append(ts.time)

b1_dist = np.array(b1_dist)
b2_dist = np.array(b2_dist)
times = np.array(times)


# ## 3. Hydrogen Bond Analysis
# 
# Analyzing hydrogen bonds between Group 27 (ROH_4GA_0GA_PHO) and Group 28 (Complex).

# In[11]:


# =================== НАСТРОЙКИ ПОД ТВОЮ СИСТЕМУ ===================
res_protein_num = 93     # Номер твоего остатка (ASH/ASP)
ligand_resname = "PHO"   # Имя твоего лиганда
# ==================================================================


# In[12]:


# ======================= НАСТРОЙКИ ПОИСКА =======================
# 1. Сюда просто вписывай массив номеров остатков белка, которые нужно проверить
protein_resids = [93, 109, 110, 170, 263, 171, 86, 262, 147, 445]  

# 2. Список твоих лигандов (уже забит по твоему запросу)
ligand_resnames = ["4GA", "0GA", "ROH"]  
# ================================================================


# In[13]:


my_resids = [93, 109, 110, 170, 263, 171, 86, 262, 147, 445, 261]
my_ligands = ["PHO"]

results_df = analyze_protein_ligand_hbonds(u, my_resids, my_ligands)


# In[22]:


my_resids = [93, 109, 110, 170, 263, 171, 86, 262, 147, 445, 261,318, 108]
my_ligands = ["4GA", "0GA", "ROH"]

results_df = analyze_protein_ligand_hbonds(u, my_resids, my_ligands)


# In[23]:


# 1. Получаем твои 2 датафрейма (например, разбив лиганды на группы)
df_pho = analyze_protein_ligand_hbonds(u, my_resids, ["PHO"])
df_sugars = analyze_protein_ligand_hbonds(u, my_resids, ["4GA", "0GA", "ROH"])

# 2. Автоматически генерируем списки для графиков! 
# (добавим отсечку в 1%, чтобы не рисовать совсем уж случайные касания)
bonds_group_1 = df_to_bond_list(df_pho, min_occupancy=5.0)
bonds_group_2 = df_to_bond_list(df_sugars, min_occupancy=5.0)

print(f"Собрано {len(bonds_group_1)} связей для первой группы и {len(bonds_group_2)} для второй.")


# In[19]:


mapping = parse_gro_mapping("step5_1.gro")

times_g1, dists_g1 = calculate_distances_by_indices(u, mapping, bonds_group_1)
if dists_g1:
    plot_individual_bonds(times_g1, dists_g1, title_prefix="PHO & NAD Cluster")

times_g2, dists_g2 = calculate_distances_by_indices(u, mapping, bonds_group_2)
if dists_g2:
    plot_individual_bonds(times_g2, dists_g2, title_prefix="4GA, 0GA & ROH Cluster")


# In[21]:


import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
from tqdm import tqdm # Импортируем tqdm для Jupyter

# ======================= НАСТРОЙКИ ПОИСКА =======================
protein_resids = [93, 109, 110, 170, 263,108]  
water_resname = "WAT"  # Имя остатка воды (SOL, WAT, HOH)

# Порог отсечения транзитной воды (в процентах)
cutoff_occ = 10.0 
# ================================================================

protein_query = " or ".join([f"resid {r}" for r in protein_resids])

donors_prot    = f"({protein_query}) and (name O* or name N*)"
hydrogens_prot = f"({protein_query}) and name H*"
acceptors_prot = f"({protein_query}) and (name O* or name N*)"

donors_wat    = f"resname {water_resname} and (name O* or name N*)"
hydrogens_wat = f"resname {water_resname} and name H*"
acceptors_wat = f"resname {water_resname} and (name O* or name N*)"

print("1. Запуск Раунда А (Белок отдает H ---> Воде)...")
hb_round_A = HydrogenBondAnalysis(
    universe=u, donors_sel=donors_prot, hydrogens_sel=hydrogens_prot,
    acceptors_sel=acceptors_wat, d_a_cutoff=3.5, d_h_a_angle_cutoff=150, update_selections=False
)
# verbose=True автоматически включит tqdm прогресс-бар от MDAnalysis!
hb_round_A.run(verbose=True) 

print("\n2. Запуск Раунда Б (Вода отдает H ---> Белку)...")
hb_round_B = HydrogenBondAnalysis(
    universe=u, donors_sel=donors_wat, hydrogens_sel=hydrogens_wat,
    acceptors_sel=acceptors_prot, d_a_cutoff=3.5, d_h_a_angle_cutoff=150, update_selections=False
)
hb_round_B.run(verbose=True)

print("\n3. Идентификация конкретных молекул воды...")
total_frames = len(u.trajectory)
bond_counts = {}

# Обрабатываем Раунд А с прогресс-баром
if hb_round_A.results.hbonds is not None and len(hb_round_A.results.hbonds) > 0:
    for bond in tqdm(hb_round_A.results.hbonds, desc="Парсинг Раунда А"):
        frame, d_idx, h_idx, a_idx, dist, ang = bond
        p_idx = int(d_idx)
        w_idx = int(a_idx)
        w_resid = u.atoms[w_idx].resid 

        key = (p_idx, w_resid, 'Донор (отдает H воде)')
        if key not in bond_counts: bond_counts[key] = set()
        bond_counts[key].add(int(frame))

# Обрабатываем Раунд Б с прогресс-баром
if hb_round_B.results.hbonds is not None and len(hb_round_B.results.hbonds) > 0:
    for bond in tqdm(hb_round_B.results.hbonds, desc="Парсинг Раунда Б"):
        frame, d_idx, h_idx, a_idx, dist, ang = bond
        w_idx = int(d_idx)
        p_idx = int(a_idx)
        w_resid = u.atoms[w_idx].resid 

        key = (p_idx, w_resid, 'Акцептор (принимает H от воды)')
        if key not in bond_counts: bond_counts[key] = set()
        bond_counts[key].add(int(frame))

# Формируем отчет (этот этап обычно занимает миллисекунды, тут tqdm не нужен)
print("Формирование финальной таблицы...")
report = []
for (p_idx, w_resid, role), frames in bond_counts.items():
    occ = (len(frames) / total_frames) * 100

    if occ >= cutoff_occ:
        atom = u.atoms[p_idx]
        report.append({
            'Occupancy (%)': occ,
            'Protein Atom': f"{atom.resid}{atom.resname}@{atom.name}",
            'Water ID': f"{water_resname}_{w_resid}",
            'Role': role
        })

df_hb = pd.DataFrame(report)

if not df_hb.empty:
    df_hb = df_hb.sort_values(by=['Protein Atom', 'Occupancy (%)'], ascending=[True, False]).reset_index(drop=True)

    print("\n" + "="*85)
    print(f"{f'ORDERED WATER MOLECULES (Occupancy > {cutoff_occ}%)':^85}")
    print("="*85)
    print(f"{'Occupancy':<9} | {'Protein Atom':<20} | {'Water ID':<15} | {'Role':<30}")
    print("-"*85)

    for _, row in df_hb.iterrows():
        print(f"{row['Occupancy (%)']:>7.1f}% | {row['Protein Atom']:<20} | {row['Water ID']:<15} | {row['Role']:<30}")
    print("="*85)

    display(df_hb)
else:
    print(f"\n❌ Структурных молекул воды с Occupancy > {cutoff_occ}% не обнаружено. Попробуй снизить cutoff_occ.")


# ## 4. Combined Plots

# In[24]:


all_plots = list(rmsd_results.items()) + [
    ("Bond: ASH 170 (OD2) -- 4GA 448 (O4)", b1_dist),
    ("Bond: TYD 263 (OH) -- 0GA 449 (H2)", b2_dist),
]

fig, axes = plt.subplots(len(all_plots), 1, figsize=(12, 3 * len(all_plots)), sharex=False)

for i, (title, data) in enumerate(all_plots):
    ax = axes[i]
    if "Bond" in title and "Hydrogen" not in title:
        ax.plot(times, data, color='firebrick')
        ax.set_ylabel("Distance (\u00c5)")
    elif "Hydrogen" in title:
        ax.plot(times, data, color='forestgreen')
        ax.set_ylabel("Count")
    else:
        ax.plot(data[:, 1], data[:, 2], color='navy')
        ax.set_ylabel("RMSD (\u00c5)")

    ax.set_title(title)
    ax.grid(True, linestyle='--', alpha=0.6)

axes[-1].set_xlabel("Time (ps)")
plt.tight_layout()
plt.show()


# In[ ]:


get_ipython().system('echo "0" | gmx trjconv -f traj_fit.xtc -s step5_1.tpr -n index_updated.ndx -skip 20 -o trajectory.pdb')


# In[ ]:


get_ipython().system('echo "38" | gmx trjconv -f traj_fit.xtc -s step5_1.tpr -n index_updated.ndx -o trajectory_no_h20.pdb')


# In[29]:


from xmlrpc.client import ServerProxy
from IPython.display import Image
import os, sys


# In[30]:


import __main__

__main__.pymol_argv = [ 'pymol', '-x' ]

### Если вывод в графическое окно тормозит или не нужен, то:
##__main__.pymol_argv = [ 'pymol', '-cp' ]


import pymol

pymol.finish_launching()

from pymol import cmd


# In[31]:


cmd.do('''
load trajectory.pdb, my_traj
''')


# In[32]:


cmd.do('''
remove solvent and not resi 3846
''')


# In[ ]:


cmd.do('''
delete all
''')


# In[ ]:


cmd.do('''
load trajectory_no_h20.pdb, my_traj
''')


# In[33]:


import re
from pymol import cmd


def parse_atom(atom_string):
    """
    Разбивает строку '445NAD@O\'N2' или '448_4GA@O2' на селекции для атома и целого остатка.
    """
    res_part, atom_name = atom_string.split('@')

    # Регулярка: (\d+) ловит цифры (resi)
    # _? ловит подчеркивание, если оно есть (опционально)
    # ([A-Za-z0-9]+) ловит буквы/цифры в конце (resn)
    match = re.match(r"^(\d+)_?([A-Za-z0-9]+)$", res_part)

    # Ловушка для ошибок парсинга: теперь скрипт точно скажет, на какой строке он споткнулся
    if not match:
        raise ValueError(f"Ошибка парсинга: регулярное выражение не смогло обработать '{res_part}' (из строки '{atom_string}')")

    resi, resn = match.groups()

    # Кавычки вокруг atom_name спасают нас от ошибок парсинга символа штриха (')
    atom_sel = f"resi {resi} and resn {resn} and name \"{atom_name}\""
    res_sel = f"resi {resi} and resn {resn}"

    return atom_sel, res_sel

def build_interactions(group_name, bonds):
    residues_to_select = set()
    dist_objects = []

    for i, (a1_str, a2_str) in enumerate(bonds):
        atom1_sel, res1_sel = parse_atom(a1_str)
        atom2_sel, res2_sel = parse_atom(a2_str)

        # Сохраняем остатки в сет, чтобы избежать дубликатов
        residues_to_select.add(res1_sel)
        residues_to_select.add(res2_sel)

        # Строим связь-измерение (аналог Wizard -> Measurement)
        dist_name = f"{group_name}_dist_{i+1}"
        cmd.distance(dist_name, atom1_sel, atom2_sel)
        dist_objects.append(dist_name)

    # 1. Создаем объект-выделение из всех участвующих остатков
    if residues_to_select:
        group_sel_string = " or ".join(residues_to_select)
        cmd.select(f"{group_name}_residues", group_sel_string)

    # 2. Группируем созданные измерения (distances) в одну папку для удобного переключения
    if dist_objects:
        cmd.group(f"{group_name}_measurements", " ".join(dist_objects))

# Запускаем обработку обеих групп
build_interactions("Group1", bonds_group_1)
build_interactions("Group2", bonds_group_2)


# In[34]:


from pymol import cmd

def build_mn_cluster():
    # 1. Создаем общее выделение для всех участников кластера (остатки + лиганды + металл)
    cmd.select("Mn_complex_residues", "resi 169 or resi 200 or resn 0GA or resn MN or resn NAD")

    # 2. Определяем селекцию для центрального иона марганца
    mn_sel = "resn MN"

    # 3. Задаем целевые атомы для связей (расстояний)
    # Словарь форматов: "Имя_объекта_измерения": "PyMOL_селекция_целевого_атома"
    coordination_targets = {
        "Mn_S_link": "(resi 169 or resi 200) and name SG",
        "Mn_NE2_link": "(resi 169 or resi 200) and name NE2",
        "Mn_O2_sugar": "(resn 0GA or resn NAD) and name O2",
        "Mn_O3_sugar": "(resn 0GA or resn NAD) and name O3",
        "Mn_O7N_link": "(resn 0GA or resn NAD) and name O7N"
    }

    dist_objects = []

    # Строим измерения от марганца ко всем целям
    for dist_name, target_sel in coordination_targets.items():
        # cmd.distance создает объект пунктирной линии между двумя селекциями
        cmd.distance(dist_name, mn_sel, target_sel)
        dist_objects.append(dist_name)

    # 4. Группируем все созданные связи в одну папку для удобного переключения
    if dist_objects:
        cmd.group("Mn_coordination_bonds", " ".join(dist_objects))

    # Базовая настройка отображения, чтобы сразу видеть красоту


# Запускаем построение
build_mn_cluster()


# In[35]:


from pymol import cmd

def build_catalytic_contacts():
    # 1. Создаем выделение для остатков 170, 263 и сахаров, чтобы их было удобно отображать
    cmd.select("Catalytic_residues", "resi 170 or resi 263 or resi 448 or resi 449")

    # 2. Строим связь: TYR 263 (OH) <---> 0GA 449 (H2)
    cmd.distance("Tyr263_0GA_link", "resi 263 and name OH", "resi 449 and resn 0GA and name H2")

    # 3. Строим связь: ASH 170 (HD2) <---> 4GA 448 (O4)
    # Явно указываем resn ASH для надежности
    cmd.distance("Ash170_4GA_link", "resi 170 and resn ASH and name HD2", "resi 448 and resn 4GA and name O4")

    # 4. Собираем эти связи в отдельную папку в интерфейсе
    cmd.group("Ligand_H_bonds", "Tyr263_0GA_link Ash170_4GA_link")

    # 5. Визуализация: показываем остатки в виде стиков и красим пунктиры в желтый цвет
# Запускаем функцию
build_catalytic_contacts()


# In[36]:


from pymol import cmd

def build_water_glu_contact_fixed():

    wat_sel = f"(resi 3846 and resn WAT)"
    glu109_sel = "(resi 109 and resn GLU)"  # Оставили эту селекцию СТРОГО для связи

    # НОВОЕ: Создаем переменную для отображения, куда добавили 110-й остаток
    view_sel = "(resi 109 or resi 110)" 

    # Обновляем команду выделения: теперь в группе вода, 109 и 110
    cmd.select("Wat_Glu_residues_Fixed", f"{wat_sel} or {view_sel}")

    # Связь строим по-старому, используя glu109_sel
    cmd.distance("Wat3846_Glu109_link_Fixed", 
                 f"{wat_sel} and (name O or name OW)", 
                 f"{glu109_sel} and (name OE1 or name OE2)")

    cmd.group("Solvent_contacts_Fixed", "Wat3846_Glu109_link_Fixed")

    # Отрисовываем стиками сразу оба остатка (109 и 110), используя view_sel


# Запускаем скрипт
build_water_glu_contact_fixed()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




