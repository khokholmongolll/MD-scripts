import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import os

# =============================================================================
# Конфигурация
# =============================================================================
XVG_FILE = "step5_1_pullf.xvg"  # Имя вашего файла
SMOOTHING_WINDOW = 51           # Размер окна для сглаживания (должен быть нечетным)
SMOOTHING_POLYORDER = 3         # Порядок полинома для сглаживания
PULLING_VELOCITY_NM_PS = None   # Если известна скорость вытягивания (в нм/пс), укажите её для оценки работы (например, 0.01)

# =============================================================================
# 1. Парсинг .xvg файла
# =============================================================================
def parse_xvg(filename):
    """Читает .xvg файл, игнорируя строки заголовков (@ и #)"""
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Файл {filename} не найден. Убедитесь, что он находится в текущей директории.")

    time = []
    force = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # Пропускаем строки заголовков GROMACS
            if line.startswith(('@', '#')) or not line:
                continue

            parts = line.split()
            if len(parts) >= 2:
                time.append(float(parts[0]))
                force.append(float(parts[1])) # Обычно 2-я колонка - это сила в kJ/(mol*nm)

    return np.array(time), np.array(force)

# =============================================================================
# 2. Физический анализ
# =============================================================================
def analyze_force(time_ps, force_kj_mol_nm):
    """Выполняет расчет физических характеристик в единицах GROMACS"""
    
    # Базовая статистика (в родных единицах kJ/(mol*nm))
    mean_force = np.mean(force_kj_mol_nm)
    max_force = np.max(force_kj_mol_nm)
    min_force = np.min(force_kj_mol_nm)
    std_force = np.std(force_kj_mol_nm)

    # Время достижения максимальной силы
    time_at_max = time_ps[np.argmax(force_kj_mol_nm)]

    # Сглаживание кривой силы (фильтр Савицкого-Голея) для выделения тренда
    # Убедимся, что длина окна не превышает длину данных и является нечетной
    window = min(SMOOTHING_WINDOW, len(force_kj_mol_nm) if len(force_kj_mol_nm) % 2 != 0 else len(force_kj_mol_nm) - 1)
    if window < SMOOTHING_POLYORDER + 2:
        window = SMOOTHING_POLYORDER + 2

    smoothed_force = savgol_filter(force_kj_mol_nm, window_length=window, polyorder=SMOOTHING_POLYORDER)

    # Оценка механической работы (если известна скорость)
    # Work = Integral(F * dx) = Integral(F * v * dt)
    work_kj_mol = None
    if PULLING_VELOCITY_NM_PS is not None:
        # Интегрирование по времени (трапецеидальный метод)
        # dt в пс, v в нм/пс -> dx в нм. F в kJ/(mol*nm). F * dx = kJ/mol
        dt = np.diff(time_ps)
        f_mid = (force_kj_mol_nm[:-1] + force_kj_mol_nm[1:]) / 2
        dx = PULLING_VELOCITY_NM_PS * dt
        work_kj_mol = np.sum(f_mid * dx)

    return {
        'time_ps': time_ps,
        'force_kj_mol_nm': force_kj_mol_nm,
        'smoothed_force_kj_mol_nm': smoothed_force,
        'mean_force_kj_mol_nm': mean_force,
        'max_force_kj_mol_nm': max_force,
        'min_force_kj_mol_nm': min_force,
        'std_force_kj_mol_nm': std_force,
        'time_at_max_ps': time_at_max,
        'work_kj_mol': work_kj_mol
    }

# =============================================================================
# 3. Визуализация и вывод результатов
# =============================================================================
def plot_results(results):
    """Строит график силы от времени в единицах GROMACS"""
    plt.figure(figsize=(10, 6))

    # Исходные зашумленные данные (полупрозрачные)
    plt.plot(results['time_ps'], results['force_kj_mol_nm'],
             color='lightgray', label='Исходная сила (шум)', alpha=0.7)

    # Сглаженная кривая
    plt.plot(results['time_ps'], results['smoothed_force_kj_mol_nm'],
             color='blue', linewidth=2, label='Сглаженная сила (Savitzky-Golay)')

    # Отметка максимальной силы
    plt.axhline(results['max_force_kj_mol_nm'], color='red', linestyle='--', linewidth=1,
                label=f'Макс. сила: {results["max_force_kj_mol_nm"]:.2f} kJ/(mol·nm)')
    plt.axvline(results['time_at_max_ps'], color='red', linestyle=':', linewidth=1)

    plt.title('Анализ силы вытягивания (COM Pulling Force)', fontsize=14)
    plt.xlabel('Время (ps)', fontsize=12)
    plt.ylabel('Сила (kJ/(mol·nm))', fontsize=12)
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('pulling_force_analysis.png', dpi=300)
    print("График сохранен как 'pulling_force_analysis.png'")
    plt.show()

def print_report(results):
    """Выводит текстовый отчет о физическом анализе в единицах GROMACS"""
    print("\n" + "="*60)
    print("ФИЗИЧЕСКИЙ АНАЛИЗ СИЛЫ ВЫТЯГИВАНИЯ (GROMACS UNITS)")
    print("="*60)
    print(f"Всего точек данных:      {len(results['time_ps'])}")
    print(f"Длительность симуляции:  {results['time_ps'][-1] - results['time_ps'][0]:.2f} ps")
    print("-" * 60)
    print(f"Средняя сила:            {results['mean_force_kj_mol_nm']:>8.2f} kJ/(mol·nm)")
    print(f"Максимальная сила:       {results['max_force_kj_mol_nm']:>8.2f} kJ/(mol·nm)  (при t = {results['time_at_max_ps']:.2f} ps)")
    print(f"Минимальная сила:        {results['min_force_kj_mol_nm']:>8.2f} kJ/(mol·nm)")
    print(f"Стд. отклонение (шум):   {results['std_force_kj_mol_nm']:>8.2f} kJ/(mol·nm)")

    if results['work_kj_mol'] is not None:
        print("-" * 60)
        print(f"Оценочная механическая работа: {results['work_kj_mol']:.4f} kJ/mol")
        print(f"  (при скорости вытягивания {PULLING_VELOCITY_NM_PS} nm/ps)")
    else:
        print("-" * 60)
        print("Оценка работы не выполнена: не указана скорость вытягивания (PULLING_VELOCITY_NM_PS).")

    print("="*60 + "\n")

# =============================================================================
# 4. Запуск
# =============================================================================
if __name__ == "__main__":
    try:
        print(f"Чтение файла: {XVG_FILE} ...")
        time_ps, force_kj_mol_nm = parse_xvg(XVG_FILE)

        print("Выполнение физического анализа...")
        results = analyze_force(time_ps, force_kj_mol_nm)

        print_report(results)
        plot_results(results)

    except Exception as e:
        print(f"Ошибка при анализе: {e}")
