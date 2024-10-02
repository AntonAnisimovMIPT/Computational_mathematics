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

fn main() {

    let matrix = DMatrix::from_row_slice(3, 3, &[
        4.0, 1.0, 1.0,
        2.0, 7.0, 1.0,
        1.0, -3.0, 12.0
    ]);

    let f = DVector::from_row_slice(&[5.0, 4.0, 10.0]);

    let x = gauss_pivot_solve(&matrix, &f);
    println!("Gauss method: x = {}", x);

    let x = lu_solve(&matrix, &f);
    println!("LU method: x = {}", x);

    let tolerance = 1e-6;
    let max_iterations = 1000000;
    let x = jacobi_solve(&matrix, &f, tolerance, max_iterations);
    println!("Jacobi method: x = {}", x);

    let x = seidel_solve(&matrix, &f, tolerance, max_iterations);
    println!("Seidel method: x = {}", x);

    let omega = 1.00;
    let x = upper_relaxation_solve(&matrix, &f, omega, tolerance, max_iterations);
    println!("Upper relaxation method: x = {}", x);

    let learning_rate = 0.01;
    let (x, residuals) = gradient_descent_solve(&matrix, &f, learning_rate, tolerance, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/gradient_descent.csv");
    println!("Gradient descent method: x = {}", x);

    let (x, residuals) = minimal_residuals_solve(&matrix, &f, tolerance, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/minimals_residual.csv");
    println!("Minimal residual method: x = {}", x);

    let (x, residuals) = bicgstab_solve(&matrix, &f, tolerance, max_iterations);
    save_residuals_to_csv(&residuals, "./../../plots_data/bicgstab.csv");
    println!("BiCGSTAB method: x = {}", x);
}
