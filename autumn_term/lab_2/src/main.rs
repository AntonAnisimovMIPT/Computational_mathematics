mod methods;

use methods::plotting::save_residuals_to_csv;

use methods::gauss_pivot::gauss_pivot_solve;
use methods::lu::lu_solve;
use methods::jacobi::jacobi_solve;
use methods::seidel::seidel_solve;
use methods::upper_relaxation::upper_relaxation_solve;
use methods::gradient_descent::gradient_descent_solve;
use methods::minimal_residuals::minimal_residuals_solve;
use methods::bicgstab::bicgstab_solve;

use nalgebra::{DMatrix, DVector};

use std::fs::File;
use std::io::{self, Write};

// т.к. решал пункт д)
fn create_custom_matrix(n: usize, a: f64) -> DMatrix<f64> {
    let mut matrix = DMatrix::zeros(n, n);

    for i in 0..n {
        for j in 0..n {
            if i == j {
                matrix[(i, j)] = a;  
            } else if (i as isize - j as isize).abs() <= 2 {
                matrix[(i, j)] = 1.0;  
            }
        }
    }

    matrix
}

// т.к. решал пункт д)
fn create_vector_f(n: usize) -> DVector<f64> {
    let mut f = DVector::zeros(n);

    for i in 0..n {
        f[i] = (i + 1) as f64;  
    }

    f
}

fn save_results_to_file(file_path: &str, method_name: &str, solution: &DVector<f64>) -> io::Result<()> {
    let mut file = File::options().append(true).create(true).open(file_path)?; 
    writeln!(file, "{} method result:\n{}", method_name, solution)?;
    Ok(())
}

fn main() {

    let n = 100;
    let a = 10.0;
    let matrix = create_custom_matrix(n, a);
    let f = create_vector_f(n);

    let results_file = "./../../results/solutions.txt"; 

    let x = gauss_pivot_solve(&matrix, &f);
    save_results_to_file(results_file, "Gauss", &x).expect("Error while writing Gauss method results");

    let x = lu_solve(&matrix, &f);
    save_results_to_file(results_file, "LU", &x).expect("Error while writing LU method results");

    let tolerance = 1e-6;
    let max_iterations = 1000000;

    let x = jacobi_solve(&matrix, &f, tolerance, max_iterations);
    save_results_to_file(results_file, "Jacobi", &x).expect("Error while writing Jacobi method results");

    let x = seidel_solve(&matrix, &f, tolerance, max_iterations);
    save_results_to_file(results_file, "Seidel", &x).expect("Error while writing Seidel method results");

    let omega = 1.00;
    let x = upper_relaxation_solve(&matrix, &f, omega, tolerance, max_iterations);
    save_results_to_file(results_file, "Upper relaxation", &x).expect("Error while writing Upper relaxation method results");

    let learning_rate = 0.01;
    let (x, residuals) = gradient_descent_solve(&matrix, &f, learning_rate, tolerance, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/gradient_descent.csv");
    save_results_to_file(results_file, "Gradient descent", &x).expect("Error while writing Gradient descent method results");

    let (x, residuals) = minimal_residuals_solve(&matrix, &f, tolerance, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/minimals_residual.csv");
    save_results_to_file(results_file, "Minimal residual", &x).expect("Error while writing Minimal residual method results");

    let (x, residuals) = bicgstab_solve(&matrix, &f, tolerance, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/bicgstab.csv");
    save_results_to_file(results_file, "BiCGSTAB", &x).expect("Error while writing BiCGSTAB method results");
}