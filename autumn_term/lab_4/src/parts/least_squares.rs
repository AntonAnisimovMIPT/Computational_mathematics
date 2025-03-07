pub struct LeastSquaresPolynomial {
    coeffs: Vec<f64>,
}

impl LeastSquaresPolynomial {
    pub fn new(x: &[f64], y: &[f64], degree: usize) -> Self {
        let mut a = vec![vec![0.0; degree + 1]; degree + 1];
        let mut b = vec![0.0; degree + 1];
        
        for i in 0..=degree {
            for j in 0..=degree {
                a[i][j] = x.iter().map(|&xi| xi.powi((i + j) as i32)).sum();
            }
            b[i] = x.iter().zip(y.iter()).map(|(&xi, &yi)| yi * xi.powi(i as i32)).sum();
        }

        let coeffs = gaussian_elimination(a, b);
        LeastSquaresPolynomial { coeffs }
    }

    pub fn evaluate(&self, x: f64) -> f64 {
        self.coeffs
            .iter()
            .enumerate()
            .map(|(i, &c)| c * x.powi(i as i32))
            .sum()
    }
}

fn gaussian_elimination(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    
    for i in 0..n {
        let mut max_row = i;
        for k in i + 1..n {
            if a[k][i].abs() > a[max_row][i].abs() {
                max_row = k;
            }
        }
        a.swap(i, max_row);
        b.swap(i, max_row);

        for k in i + 1..n {
            let factor = a[k][i] / a[i][i];
            for j in i..n {
                a[k][j] -= factor * a[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        x[i] = b[i] / a[i][i];
        for k in 0..i {
            b[k] -= a[k][i] * x[i];
        }
    }

    x
}
