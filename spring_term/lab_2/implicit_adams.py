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
    
    # Для x*(2*a1 - 0.5*x - (a1**2)/(a2**2)*y)
    J[0, 0] = 2*a1 - x - (a1**2)/(a2**2)*y_val                 # d/dx
    J[0, 1] = - x*(a1**2)/(a2**2)                                # d/dy
    J[0, 2] = x * (2 - (2*a1*y_val)/(a2**2))                     # d/da1
    J[0, 3] = 2*x*(a1**2)*y_val/(a2**3)                          # d/da2

    # Для y*(2*a2 - (a2**2)/(a1**2)*x - 0.5*y)
    J[1, 0] = - y_val*(a2**2)/(a1**2)                        
    J[1, 1] = 2*a2 - (a2**2)/(a1**2)*x - y_val                  
    J[1, 2] = y_val*(2*x*(a2**2)/(a1**3))                       
    J[1, 3] = y_val*(2 - (2*a2*x)/(a1**2))                     

    # Для eps*(2 - 2*a1*y/(a2**2))
    J[2, 0] = 0                                               
    J[2, 1] = -eps*(2*a1)/(a2**2)                              
    J[2, 2] = -eps*(2*y_val)/(a2**2)                            
    J[2, 3] = eps*(4*a1*y_val)/(a2**3)                         

    # Для eps*(2 - 2*a2*x/(a1**2))
    J[3, 0] = -eps*(2*a2)/(a1**2)                             
    J[3, 1] = 0                                               
    J[3, 2] = eps*(4*a2*x)/(a1**3)                            
    J[3, 3] = -eps*(2*x)/(a1**2)                             
    
    return J

def adams_moulton_step(y_n, h, eps, tol=1e-8, max_iter=20):
    y_next_guess = y_n + h * system(y_n, eps)  # Начальное приближение (явный Эйлер)
    
    for _ in range(max_iter):
        f_next = system(y_next_guess, eps)
        F = y_next_guess - y_n - 0.5*h*(system(y_n, eps) + f_next)
        
        if np.linalg.norm(F) < tol:
            break
        
        J = np.eye(4) - 0.5*h * jacobian_system(y_next_guess, eps)
        delta = np.linalg.solve(J, -F)
        y_next_guess += delta
    else:
        print("Ньютон не сошелся за", max_iter, "итераций")
    
    return y_next_guess

def solver(nach_usl, T, h, eps):
    N = int(T/h)
    t = np.linspace(0, T, N+1)
    sol = np.zeros((N+1, len(nach_usl)))
    sol[0] = nach_usl
    y_current = nach_usl.copy()
    
    for i in range(N):
        y_next = adams_moulton_step(y_current, h, eps)
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

t, sol = solver(nach_usl, T, h, eps)

# Построение графиков
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(t, sol[:, 0], 'b-')
plt.xlabel("t")
plt.ylabel("x")
plt.title("x(t)")

plt.subplot(2, 2, 2)
plt.plot(t, sol[:, 1], 'r-')
plt.xlabel("t")
plt.ylabel("y")
plt.title("y(t)")

plt.subplot(2, 2, 3)
plt.plot(t, sol[:, 2], 'g-')
plt.xlabel("t")
plt.ylabel("a1")
plt.title("a1(t)")

plt.subplot(2, 2, 4)
plt.plot(t, sol[:, 3], 'm-')
plt.xlabel("t")
plt.ylabel("a2")
plt.title("a2(t)")

plt.tight_layout()
plt.savefig("implicit_adams.png")