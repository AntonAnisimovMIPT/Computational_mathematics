import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

##
## уравнение вида du/dt + c * du/dx = 0,  x ∈ [0, L], t ∈ [0, T]
## c начальным условием u(x, 0) = u_0(x) и периодическими граничными
## имеет решение вида u(x, t) = u_0(x - ct)
##

L = 10.0      
T = 0.2       
c = 1.0      

def analytical_solution(x, t):
    return np.exp(-0.5 * ((x - c * t - L/2) / 0.5)**2)

def initial_condition(x):
    return analytical_solution(x, 0)

def upwind_scheme(u, c, tau, h):
    u_new = u.copy()
    u_new[1:] = u[1:] - c * tau / h * (u[1:] - u[:-1])
    u_new[0] = u[0] - c * tau / h * (u[0] - u[-1])  # Периодическое ГУ
    return u_new

def ftcs_scheme(u, c, tau, h):
    u_new = u.copy()
    u_new[1:-1] = u[1:-1] - c * tau / (2 * h) * (u[2:] - u[:-2])
    # Периодические ГУ
    u_new[0] = u[0] - c * tau / (2 * h) * (u[1] - u[-1])
    u_new[-1] = u[-1] - c * tau / (2 * h) * (u[0] - u[-2])
    return u_new

def lax_wendroff_scheme(u, c, tau, h):
    u_new = u.copy()
    sigma = c * tau / h
    sigma2 = sigma**2
    u_new[1:-1] = (u[1:-1] - sigma/2 * (u[2:] - u[:-2]) + 
                   sigma2/2 * (u[2:] - 2*u[1:-1] + u[:-2]))
    # Периодические ГУ
    u_new[0] = (u[0] - sigma/2 * (u[1] - u[-1]) + 
                sigma2/2 * (u[1] - 2*u[0] + u[-1]))
    u_new[-1] = (u[-1] - sigma/2 * (u[0] - u[-2]) + 
                 sigma2/2 * (u[0] - 2*u[-1] + u[-2]))
    return u_new

Nx = 100    
h = L / (Nx - 1) 
CFL = 0.5      
tau = CFL * h / c 
Nt = int(T / tau) 

x = np.linspace(0, L, Nx)
u_analytical = initial_condition(x)
u_upwind = u_analytical.copy()
u_ftcs = u_analytical.copy()
u_lw = u_analytical.copy()

plt.figure(figsize=(12, 8))
plt.title("Сравнение численных схем с аналитическим решением")
plt.xlabel("x")
plt.ylabel("u")

u_upwind_history = [u_upwind.copy()]
u_ftcs_history = [u_ftcs.copy()]
u_lw_history = [u_lw.copy()]
u_analytical_history = [u_analytical.copy()]


for n in range(Nt):
    t = (n + 1) * tau
    
    u_analytical = analytical_solution(x, t)
    
    u_upwind = upwind_scheme(u_upwind, c, tau, h)
    u_ftcs = ftcs_scheme(u_ftcs, c, tau, h)
    u_lw = lax_wendroff_scheme(u_lw, c, tau, h)
    
    u_upwind_history.append(u_upwind.copy())
    u_ftcs_history.append(u_ftcs.copy())
    u_lw_history.append(u_lw.copy())
    u_analytical_history.append(u_analytical.copy())


plt.plot(x, u_analytical_history[-1], 'k-', label='Аналитическое решение', linewidth=2)
plt.plot(x, u_upwind_history[-1], 'r--', label='уголок')
plt.plot(x, u_ftcs_history[-1], 'b-.', label='прямоугольник')
plt.plot(x, u_lw_history[-1], 'g:', label='явная 4х точечная')
plt.legend()
plt.grid()
plt.savefig('comparison.png', dpi=300)
plt.close()

Nx_list = [20, 40, 120, 140, 200, 400, 450, 500, 600] 
errors_upwind = []
errors_ftcs = []
errors_lw = []

for Nx in Nx_list:
    h = L / (Nx - 1)
    tau = CFL * h / c
    Nt = int(T / tau)
    
    x = np.linspace(0, L, Nx)
    u_analytical = analytical_solution(x, T)
    
    u_upwind = initial_condition(x)
    for n in range(Nt):
        u_upwind = upwind_scheme(u_upwind, c, tau, h)
    error = np.linalg.norm(u_upwind - u_analytical) / np.sqrt(Nx)
    errors_upwind.append(error)
    
    u_ftcs = initial_condition(x)
    for n in range(Nt):
        u_ftcs = ftcs_scheme(u_ftcs, c, tau, h)
    error = np.linalg.norm(u_ftcs - u_analytical) / np.sqrt(Nx)
    errors_ftcs.append(error)
    
    u_lw = initial_condition(x)
    for n in range(Nt):
        u_lw = lax_wendroff_scheme(u_lw, c, tau, h)
    error = np.linalg.norm(u_lw - u_analytical) / np.sqrt(Nx)
    errors_lw.append(error)

h_list = [L / (Nx - 1) for Nx in Nx_list]
plt.figure(figsize=(10, 6))
plt.loglog(h_list, errors_upwind, 'ro-', label='уголок')
plt.loglog(h_list, errors_ftcs, 'bs-', label='прямоугольник')
plt.loglog(h_list, errors_lw, 'g^-', label='явная 4х точечная')

h_theory = np.array(h_list)
plt.loglog(h_theory, 10*h_theory, 'k--', label='O(h)')
plt.loglog(h_theory, 10*h_theory**2, 'k:', label='O(h²)')

plt.xlabel('Шаг сетки h')
plt.ylabel('Норма ошибки')
plt.title('Исследование сходимости численных схем')
plt.legend()
plt.grid(True, which="both", ls="-")
plt.savefig('convergence.png', dpi=300)
plt.close()