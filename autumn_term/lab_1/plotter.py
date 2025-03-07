import os
import matplotlib.pyplot as plt
import pandas as pd

def plot_data(directory, function_name, output_dir):
    
    csv_files = [f for f in os.listdir(directory) if f.endswith('.csv')]
    
    plt.figure(figsize=(8, 6))
    
    for csv_file in csv_files:
        filepath = os.path.join(directory, csv_file)
        data = pd.read_csv(filepath)
        
        label = csv_file.replace('.csv', '')
        plt.plot(data['h'], data['error'], label=label)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Step size (h)')
    plt.ylabel('Absolute error')
    plt.title(f'Error vs Step Size for {function_name}')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    plot_filename = os.path.join(output_dir, f'{function_name}_error_plot.png')
    plt.savefig(plot_filename)
    plt.close()

def main():
    base_dir = 'calculated_data'  
    output_dir = 'plots' 

    if not os.path.exists(base_dir):
        print(f"Ошибка: Путь {base_dir} не существует.")
        return
    
    function_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    
    for function_dir in function_dirs:
        function_path = os.path.join(base_dir, function_dir)
        function_name = function_dir  
        plot_data(function_path, function_name, output_dir)
        print(f'График для {function_name} сохранен.')

if __name__ == "__main__":
    main()
