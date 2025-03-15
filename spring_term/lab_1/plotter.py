import matplotlib.pyplot as plt

def plot_results(filename, title, output_filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    t = []
    x = []
    z = []
    
    for line in lines:
        parts = line.split()
        t.append(float(parts[0]))
        x.append(float(parts[1]))
        z.append(float(parts[2]))
    
    plt.figure(figsize=(10, 6))
    plt.plot(t, x, label='x(t)')
    plt.plot(t, z, label='z(t)')
    plt.title(title)
    plt.xlabel('Time')
    plt.ylabel('Values')
    plt.legend()
    plt.grid(True)
    
    plt.savefig(output_filename)
    plt.close()

plot_results('./results/data/runge_kutta1.txt', 'Runge-Kutta1 Method', './results/plots/runge_kutta1.png')
plot_results('./results/data/runge_kutta2.txt', 'Runge-Kutta2 Method', './results/plots/runge_kutta2.png')
plot_results('./results/data/runge_kutta3.txt', 'Runge-Kutta3 Method', './results/plots/runge_kutta3.png')
plot_results('./results/data/runge_kutta4.txt', 'Runge-Kutta4 Method', './results/plots/runge_kutta4.png')


plot_results('./results/data/adams2.txt', 'Adams2 Method', './results/plots/adams2.png')
plot_results('./results/data/adams3.txt', 'Adams3 Method', './results/plots/adams3.png')
plot_results('./results/data/adams4.txt', 'Adams4 Method', './results/plots/adams4.png')

plot_results('./results/data/bdf2.txt', 'BDF2 Method', './results/plots/bdf2.png')
plot_results('./results/data/bdf3.txt', 'BDF3 Method', './results/plots/bdf3.png')
plot_results('./results/data/bdf4.txt', 'BDF4 Method', './results/plots/bdf4.png')