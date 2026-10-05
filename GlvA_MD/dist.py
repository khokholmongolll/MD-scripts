import numpy as np
import matplotlib.pyplot as plt
import os

# =============================================================================
# Конфигурация
# =============================================================================
# ВНИМАНИЕ: Здесь должен быть файл с РАССТОЯНИЕМ, а не с силой!
# Обычно это step5_1_pullx.xvg или имя, которое вы задали в pull-coord1-output
XVG_FILE = "step5_1_pullx.xvg" 
TEMPERATURE_K = 300.0

# Физические константы
KB = 0.008314462618 # kJ / (mol * K)
KT = KB * TEMPERATURE_K # kJ / mol

# =============================================================================
# 1. Парсинг .xvg файла (расстояния)
# =============================================================================
def parse_xvg_distance(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Файл {filename} не найден. Проверьте имя файла с расстояниями (обычно pullx.xvg).")
    
    time = []
    distance = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(('@', '#')) or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                time.append(float(parts[0]))
                distance.append(float(parts[1])) # 2-я колонка - расстояние в нм
                
    return np.array(time), np.array(distance)

# =============================================================================
# 2. Анализ
# =============================================================================
def analyze_thermal_fluctuations(time_ps, distance_nm):
    mean_dist = np.mean(distance_nm)
    std_dist = np.std(distance_nm) # Это и есть амплитуда тепловых флуктуаций (RMSF)
    min_dist = np.min(distance_nm)
    max_dist = np.max(distance_nm)
    
    # Оценка "жесткости" естественного потенциала удержания (из флуктуаций)
    # Согласно теореме о равнораспределении: 0.5 * k_eff * <dx^2> = 0.5 * k_B T
    # Следовательно, эффективная жесткость окружения k_eff = k_B T / variance
    variance = std_dist ** 2
    if variance > 0:
        k_eff = KT / variance # в кДж / (mol * nm^2)
    else:
        k_eff = 0
        
    return {
        'mean_dist_nm': mean_dist,
        'std_dist_nm': std_dist,
        'min_dist_nm': min_dist,
        'max_dist_nm': max_dist,
        'k_eff': k_eff
    }

# =============================================================================
# 3. Вывод
# =============================================================================
def print_report(results):
    print("\n" + "="*60)
    print("АНАЛИЗ ЕСТЕСТВЕННЫХ ТЕПЛОВЫХ ФЛУКТУАЦИЙ РАССТОЯНИЯ")
    print("="*60)
    print(f"Среднее расстояние:      {results['mean_dist_nm']:.4f} нм")
    print(f"Минимальное расстояние:  {results['min_dist_nm']:.4f} нм")
    print(f"Максимальное расстояние: {results['max_dist_nm']:.4f} нм")
    print("-" * 60)
    print(f"Амплитуда флуктуаций (RMSF): {results['std_dist_nm']:.4f} нм ({results['std_dist_nm']*10:.2f} Å)")
    
    # Переводим флуктуацию в энергетический контекст
    print(f"\nЭнергетическая интерпретация:")
    print(f"  k_B T при {TEMPERATURE_K} K = {KT:.3f} кДж/моль")
    print(f"  Эффективная жесткость окружения (k_eff): {results['k_eff']:.1f} кДж/(моль·нм²)")
    print("  (Это показывает, насколько 'жестко' белок/растворитель удерживают контакт сами по себе)")
    print("="*60 + "\n")

def plot_results(time_ps, distance_nm, results):
    plt.figure(figsize=(10, 5))
    plt.plot(time_ps, distance_nm, color='blue', alpha=0.5, linewidth=0.5, label='Расстояние')
    plt.axhline(results['mean_dist_nm'], color='red', linestyle='--', label=f"Среднее: {results['mean_dist_nm']:.3f} нм")
    
    # Область тепловых флуктуаций +/- 1 std
    plt.fill_between(time_ps, 
                     results['mean_dist_nm'] - results['std_dist_nm'], 
                     results['mean_dist_nm'] + results['std_dist_nm'], 
                     color='red', alpha=0.2, label='±1 σ (тепловые флуктуации)')
    
    # Линия активации flat-bottom (для справки)
    plt.axhline(0.45, color='green', linestyle=':', linewidth=2, label='Граница flat-bottom (0.45 нм)')
    
    plt.title('Тепловые флуктуации расстояния между группами', fontsize=14)
    plt.xlabel('Время (пс)', fontsize=12)
    plt.ylabel('Расстояние (нм)', fontsize=12)
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('thermal_distance_fluctuations.png', dpi=300)
    print("График сохранен как 'thermal_distance_fluctuations.png'")
    plt.show()

# =============================================================================
# 4. Запуск
# =============================================================================
if __name__ == "__main__":
    try:
        print(f"Чтение файла расстояний: {XVG_FILE} ...")
        time_ps, distance_nm = parse_xvg_distance(XVG_FILE)
        results = analyze_thermal_fluctuations(time_ps, distance_nm)
        print_report(results)
        plot_results(time_ps, distance_nm, results)
    except Exception as e:
        print(f"Ошибка: {e}")
        print("Убедитесь, что файл step5_1_pullx.xvg существует в текущей папке.")
