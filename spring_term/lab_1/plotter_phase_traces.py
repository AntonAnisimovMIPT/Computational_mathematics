import matplotlib.pyplot as plt
import numpy as np

# Значения e
e_values = [0.01, 0.1, 0.2, 0.5, 0.7, 0.9, 1, 2, 3, 4, 5]

# Построение графиков для каждого e
for e in e_values:
    # Загрузка данных
    data = np.loadtxt(f'./results/data/phase_traces/runge_kutta4_e{e}.txt')
    t, x, z = data[:, 0], data[:, 1], data[:, 2]

    # Построение фазовой траектории
    plt.figure()
    plt.plot(x, z, label=f'e = {e}')
    plt.xlabel('x')
    plt.ylabel('z')
    plt.title(f'Фазовая траектория (e = {e})')
    plt.legend()
    plt.grid()
    plt.savefig(f'./results/plots/phase_traces/phase_trajectory_e{e}.png')
    plt.close()