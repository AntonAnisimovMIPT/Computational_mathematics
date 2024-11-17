pub fn divided_differences(x: &[f64], y: &[f64]) -> Vec<Vec<f64>> {
    let n = x.len();
    let mut table = vec![vec![0.0; n]; n];

    for i in 0..n {
        table[i][0] = y[i];
    }

    for j in 1..n {
        for i in 0..(n - j) {
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x[i + j] - x[i]);
        }
    }

    table
}


pub fn newton_polynomial(x: &[f64], target: f64, diff_table: &[Vec<f64>]) -> f64 {
    let n = x.len();
    let mut result = diff_table[0][0]; 
    let mut product_term = 1.0;

    for i in 1..n {
        product_term *= target - x[i - 1];
        result += diff_table[0][i] * product_term;
    }

    result
}
