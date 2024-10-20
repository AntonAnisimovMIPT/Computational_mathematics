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

use std::fs;
use std::io::Write;
use std::error::Error;
use std::fs::{OpenOptions, create_dir_all};
use std::path::Path;

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

fn compute_residual_norm(matrix: &DMatrix<f64>, x: &DVector<f64>, f: &DVector<f64>) -> f64 {
    let residual = f - matrix * x;
    residual.norm()
}

fn save_solution_to_csv(file_path: &str, method_name: &str, solution: &DVector<f64>) -> Result<(), Box<dyn Error>> {

    if let Some(parent) = Path::new(file_path).parent() {
        create_dir_all(parent)?;
    }

    let mut file = OpenOptions::new().append(true).create(true).open(file_path)?;

    writeln!(file, "{}", method_name)?;

    let solution_line = solution.iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>()
        .join(", ");
    
    writeln!(file, "{}", solution_line)?;

    Ok(())
}

const MACHINE_EPSILON: f64 = 1e-12;

fn main() {

    let n = 100;
    let a = 10.0;
    let matrix = create_custom_matrix(n, a);
    let f = create_vector_f(n);

    let results_file = "./../../results/solutions.csv";
    if Path::new(results_file).exists() {
        fs::remove_file(results_file).expect("Failed to delete previous results file");
    }

    let x = gauss_pivot_solve(&matrix, &f);
    let residual_norm = compute_residual_norm(&matrix, &x, &f);
    assert!(residual_norm <= MACHINE_EPSILON, "Residual too large for Gauss method");
    save_solution_to_csv(results_file, "Gauss", &x).expect("Error while writing Gauss method results");

    let x = lu_solve(&matrix, &f);
    let residual_norm = compute_residual_norm(&matrix, &x, &f);
    assert!(residual_norm <= MACHINE_EPSILON, "Residual too large for LU method");
    save_solution_to_csv(results_file, "LU", &x).expect("Error while writing LU method results");

    let tolerance_for_iters = 1e-12;
    let max_iterations = 1000000;

    let (x, residuals) = jacobi_solve(&matrix, &f, tolerance_for_iters, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/jacobi.csv");
    save_solution_to_csv(results_file, "Jacobi", &x).expect("Error while writing Jacobi method results");

    let (x, residuals) = seidel_solve(&matrix, &f, tolerance_for_iters, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/seidel.csv");
    save_solution_to_csv(results_file, "Seidel", &x).expect("Error while writing Seidel method results");

    let omega = 1.00;
    let (x, residuals) = upper_relaxation_solve(&matrix, &f, omega, tolerance_for_iters, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/upper_relaxation.csv");
    save_solution_to_csv(results_file, "Upper relaxation", &x).expect("Error while writing Upper relaxation method results");

    let learning_rate = 0.01;
    let (x, residuals) = gradient_descent_solve(&matrix, &f, learning_rate, tolerance_for_iters, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/gradient_descent.csv");
    save_solution_to_csv(results_file, "Gradient descent", &x).expect("Error while writing Gradient descent method results");

    let (x, residuals) = minimal_residuals_solve(&matrix, &f, tolerance_for_iters, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/minimals_residual.csv");
    save_solution_to_csv(results_file, "Minimal residual", &x).expect("Error while writing Minimal residual method results");

    let (x, residuals) = bicgstab_solve(&matrix, &f, tolerance_for_iters, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/bicgstab.csv");
    save_solution_to_csv(results_file, "BiCGSTAB", &x).expect("Error while writing BiCGSTAB method results");

}
