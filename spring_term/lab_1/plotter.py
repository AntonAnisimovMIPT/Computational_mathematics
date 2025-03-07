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

plot_results('./results/data/runge_kutta.txt', 'Runge-Kutta Method', './results/plots/runge_kutta.png')
plot_results('./results/data/adams.txt', 'Adams Method', './results/plots/adams.png')
plot_results('./results/data/fdb.txt', 'FDB Method', './results/plots/fdb.png')