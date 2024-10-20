use nalgebra::{DMatrix, DVector};

pub fn seidel_solve(matrix: &DMatrix<f64>, b: &DVector<f64>, tolerance: f64, max_iterations: usize) -> (DVector<f64>, Vec<f64>) {
    let n = matrix.nrows();
    let mut x = DVector::zeros(n);
    let mut residuals = Vec::new(); 

    for iteration in 0..max_iterations {
        let mut x_new = x.clone(); 

        for i in 0..n {
            let mut sum1 = 0.0;
            let mut sum2 = 0.0;

            for j in 0..i {
                sum1 += matrix[(i, j)] * x_new[j];
            }

            for j in i + 1..n {
                sum2 += matrix[(i, j)] * x[j];
            }

            x_new[i] = (b[i] - sum1 - sum2) / matrix[(i, i)];
        }

        let residual = b - matrix * &x_new;
        residuals.push(residual.norm());

        if (&x_new - &x).norm() < tolerance {
            println!("The Seidel method converged in {} iterations", iteration + 1);
            return (x_new, residuals);
        }

        x = x_new;
    }

    println!("Warning!!! The Seidel method did not converge for {} iterations.", max_iterations);
    (x, residuals) 
}
