import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Параметры задачи
eps = 0.01
x0 = 1 # Пример начального значения x
y0 = 1  # Пример начального значения y
a10 = 0.001  # Пример начального значения a1
a20 = 10  # Начальное значение a2
T = 2000  # Время интегрирования

# Система ОДУ
def system(t, variables, eps):
    x, y, a1, a2 = variables
    dxdt = x * (2 * a1 - 0.5 * x - (a1**2) * (a2**(-2)) * y)
    dydt = y * (2 * a2 - (a1**(-2)) * (a2**2) * x - 0.5 * y)
    da1dt = eps * (2 - 2 * a1 * (a2**(-2)) * y)
    da2dt = eps * (2 - 2 * (a1**(-2)) * a2 * x)
    return [dxdt, dydt, da1dt, da2dt]

# Начальные условия
initial_conditions = [x0, y0, a10, a20]

# Временной интервал
t_span = (0, T)

# Решение системы ОДУ
sol = solve_ivp(system, t_span, initial_conditions, args=(eps,), method='Radau', t_eval=np.linspace(0, T, 10000))

# Построение графиков
plt.figure(figsize=(12, 8))

# График x(t)
plt.subplot(2, 2, 1)
plt.plot(sol.t, sol.y[0], label='x(t)')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.legend()

# График y(t)
plt.subplot(2, 2, 2)
plt.plot(sol.t, sol.y[1], label='y(t)', color='orange')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.legend()

# График a1(t)
plt.subplot(2, 2, 3)
plt.plot(sol.t, sol.y[2], label='a1(t)', color='green')
plt.xlabel('t')
plt.ylabel('a1(t)')
plt.legend()

# График a2(t)
plt.subplot(2, 2, 4)
plt.plot(sol.t, sol.y[3], label='a2(t)', color='red')
plt.xlabel('t')
plt.ylabel('a2(t)')
plt.legend()

plt.tight_layout()
plt.savefig("output.png")
