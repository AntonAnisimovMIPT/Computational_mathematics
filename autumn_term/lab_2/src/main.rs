mod methods;

use methods::gauss_pivot::gauss_pivot_solve;
use methods::lu::lu_solve;
use methods::jacobi::jacobi_solve;
use methods::seidel::seidel_solve;
use methods::upper_relaxation::upper_relaxation_solve;

use nalgebra::{DMatrix, DVector};

fn main() {

    // Создание матрицы 3x3
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
}
