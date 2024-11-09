pub fn divided_differences(years: &Vec<f64>, populations: &Vec<f64>) -> Vec<Vec<f64>> {
    let n = years.len();
    let mut diff_table = vec![vec![0.0; n]; n];

    for i in 0..n {
        diff_table[i][0] = populations[i];
    }

    for j in 1..n {
        for i in 0..(n - j) {
            diff_table[i][j] = (diff_table[i + 1][j - 1] - diff_table[i][j - 1]) /
                               (years[i + j] - years[i]);
        }
    }

    diff_table
}

pub fn newton_polynomial(years: &Vec<f64>, x: f64, diff_table: &Vec<Vec<f64>>) -> f64 {
    let mut result = diff_table[0][0];
    let mut product = 1.0;

    for i in 1..years.len() {
        product *= x - years[i - 1];
        result += diff_table[0][i] * product;
    }

    result
}
