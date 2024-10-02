use nalgebra::{DMatrix, DVector};

pub fn jacobi_solve(matrix: &DMatrix<f64>, b: &DVector<f64>, tolerance: f64, max_iterations: usize) -> DVector<f64> {
    let n = matrix.nrows();
    let mut x = DVector::zeros(n); 
    let mut x_new = x.clone(); 

    for iteration in 0..max_iterations {
        for i in 0..n {
            let mut sum = 0.0;

            for j in 0..n {
                if i != j {
                    sum += matrix[(i, j)] * x[j];
                }
            }

            x_new[i] = (b[i] - sum) / matrix[(i, i)];
        }

        if (&x_new - &x).norm() < tolerance {
            println!("The Jacobi method converged in {} iterations", iteration + 1);
            return x_new;
        }

        x = x_new.clone();
    }

    println!("Warning!!! The Jacobi method did not converge for {} iterations.", max_iterations);
    x_new
}