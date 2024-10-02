import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_residuals(csv_path, output_dir):
    filename = os.path.basename(csv_path).replace(".csv", "")
    
    data = pd.read_csv(csv_path)
    
    plt.figure()
    plt.plot(data['iteration'], data['residual'], marker='o', label=filename)
    plt.xlabel('Iteration')
    plt.ylabel('Residual')
    plt.title(f'Residuals Plot - {filename}')
    plt.grid(True)
    plt.legend()
    
    output_path = os.path.join(output_dir, f"{filename}.png")
    plt.savefig(output_path)
    plt.close()

def process_all_csv(input_dir, output_dir):
    if not os.path.exists(input_dir):
        print(f"Input directory {input_dir} does not exist.")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    csv_files = [f for f in os.listdir(input_dir) if f.endswith(".csv")]

    if not csv_files:
        print(f"No CSV files found in {input_dir}.")
        return

    for csv_file in csv_files:
        csv_path = os.path.join(input_dir, csv_file)
        print(f"Processing {csv_path}...")
        plot_residuals(csv_path, output_dir)
    
    print(f"All graphs saved to {output_dir}.")

if __name__ == "__main__":
    input_directory = "./../plots_data"  
    output_directory = "."  
    process_all_csv(input_directory, output_directory)
