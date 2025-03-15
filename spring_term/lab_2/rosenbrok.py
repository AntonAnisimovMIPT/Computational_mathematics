import numpy as np
import matplotlib.pyplot as plt

def system(vars, eps):
    x, y_val, a1, a2 = vars
    dxdt = x * (2*a1 - 0.5*x - (a1**2)/(a2**2)*y_val)
    dydt = y_val * (2*a2 - (a2**2)/(a1**2)*x - 0.5*y_val)
    da1dt = eps * (2 - 2*a1*y_val/(a2**2))
    da2dt = eps * (2 - 2*a2*x/(a1**2))
    return np.array([dxdt, dydt, da1dt, da2dt])

def jacobian_system(vars, eps):
    x, y_val, a1, a2 = vars
    J = np.zeros((4, 4))
    
    # Производные для dxdt
    J[0,0] = 2*a1 - x - (a1**2)/(a2**2)*y_val
    J[0,1] = -x*(a1**2)/(a2**2)
    J[0,2] = x*(2 - 2*a1*y_val/(a2**2))
    J[0,3] = 2*x*a1**2*y_val/(a2**3)
    
    # Производные для dydt
    J[1,0] = -y_val*(a2**2)/(a1**2)
    J[1,1] = 2*a2 - y_val - (a2**2)/(a1**2)*x
    J[1,2] = 2*y_val*a2**2*x/(a1**3)
    J[1,3] = y_val*(2 - 2*a2*x/(a1**2))
    
    # Производные для da1dt
    J[2,1] = -eps*2*a1/(a2**2)
    J[2,2] = -eps*2*y_val/(a2**2)
    J[2,3] = eps*4*a1*y_val/(a2**3)
    
    # Производные для da2dt
    J[3,0] = -eps*2*a2/(a1**2)
    J[3,2] = eps*4*a2*x/(a1**3)
    J[3,3] = -eps*2*x/(a1**2)
    
    return J

# Реализация метода Розенброка-Ваннера 2-го порядка
def rosenbrock2_step(y_n, h, eps):
    gamma = 1.7071067811865476  # 1 + 1/√2
    alpha21 = 1.0
    gamma21 = -gamma
    m1 = 0.8786796564403576
    m2 = 0.2928932188134524
    
    J = jacobian_system(y_n, eps)
    I = np.eye(4)
    A = I - h * gamma * J
    
    # Вычисление стадий
    f_n = system(y_n, eps)
    k1 = np.linalg.solve(A, h * f_n)
    
    y_temp = y_n + alpha21 * k1
    f_temp = system(y_temp, eps)
    rhs = h * (f_temp + gamma21 * k1 / h)
    k2 = np.linalg.solve(A, rhs)
    
    y_next = y_n + m1 * k1 + m2 * k2
    return y_next

# Реализация метода Розенброка-Ваннера 3-го порядка (ROS3P)
def rosenbrock3_step(y_n, h, eps):
    gamma = 0.435866521508
    a21 = 1.0
    a31 = 1.0
    a32 = 0.0
    c21 =-1.382799922
    c31 = 0.327682820
    c32 = 0.523979068
    m1 = 0.476612548
    m2 = 0.098511733
    m3 = 0.424875719
    
    J = jacobian_system(y_n, eps)
    I = np.eye(4)
    A = I / gamma - h * J
    
    # Стадия 1
    f_n = system(y_n, eps)
    k1 = np.linalg.solve(A, f_n)
    
    # Стадия 2
    y_temp = y_n + a21 * h * k1
    f_temp = system(y_temp, eps)
    rhs = f_temp + c21 * k1 / h
    k2 = np.linalg.solve(A, rhs)
    
    # Стадия 3
    y_temp = y_n + a31 * h * k1 + a32 * h * k2
    f_temp = system(y_temp, eps)
    rhs = f_temp + c31 * k1 / h + c32 * k2 / h
    k3 = np.linalg.solve(A, rhs)
    
    y_next = y_n + h * (m1 * k1 + m2 * k2 + m3 * k3)
    return y_next

def solver(nach_usl, T, h, eps, method='row2'):
    N = int(T/h)
    t = np.linspace(0, T, N+1)
    sol = np.zeros((N+1, len(nach_usl)))
    sol[0] = nach_usl
    y_current = nach_usl.copy()
    
    step_func = rosenbrock2_step if method == 'row2' else rosenbrock3_step
    
    for i in range(N):
        y_next = step_func(y_current, h, eps)
        sol[i+1] = y_next
        y_current = y_next
        if i % 100 == 0:
            print(f"Step {i}/{N}")
    
    return t, sol

# Параметры
eps = 0.000001
nach_usl = np.array([20.0, 20.0, 0.005, 10.0])
T = 2000 
h = 0.001  

# Решение методом ROW2
t_row2, sol_row2 = solver(nach_usl, T, h, eps, method='row2')

# Решение методом ROW3
t_row3, sol_row3 = solver(nach_usl, T, h, eps, method='row3')

# Построение графиков для ROW2
plt.figure(figsize=(12, 8))
plt.suptitle("Метод Розенброка-Ваннера 2-го порядка")

plt.subplot(2, 2, 1)
plt.plot(t_row2, sol_row2[:, 0], 'b-')
plt.xlabel("t")
plt.ylabel("x")
plt.title("x(t)")

plt.subplot(2, 2, 2)
plt.plot(t_row2, sol_row2[:, 1], 'r-')
plt.xlabel("t")
plt.ylabel("y")
plt.title("y(t)")

plt.subplot(2, 2, 3)
plt.plot(t_row2, sol_row2[:, 2], 'g-')
plt.xlabel("t")
plt.ylabel("a1")
plt.title("a1(t)")

plt.subplot(2, 2, 4)
plt.plot(t_row2, sol_row2[:, 3], 'm-')
plt.xlabel("t")
plt.ylabel("a2")
plt.title("a2(t)")

plt.tight_layout()
plt.savefig("rosenbrock2.png")

# Построение графиков для ROW3
plt.figure(figsize=(12, 8))
plt.suptitle("Метод Розенброка-Ваннера 3-го порядка")

plt.subplot(2, 2, 1)
plt.plot(t_row3, sol_row3[:, 0], 'b-')
plt.xlabel("t")
plt.ylabel("x")
plt.title("x(t)")

plt.subplot(2, 2, 2)
plt.plot(t_row3, sol_row3[:, 1], 'r-')
plt.xlabel("t")
plt.ylabel("y")
plt.title("y(t)")

plt.subplot(2, 2, 3)
plt.plot(t_row3, sol_row3[:, 2], 'g-')
plt.xlabel("t")
plt.ylabel("a1")
plt.title("a1(t)")

plt.subplot(2, 2, 4)
plt.plot(t_row3, sol_row3[:, 3], 'm-')
plt.xlabel("t")
plt.ylabel("a2")
plt.title("a2(t)")

plt.tight_layout()
plt.savefig("rosenbrock3.png")