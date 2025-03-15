use crate::runge_kutta4;

pub fn adams2<F>(
    f: F,
    x0: f64,
    z0: f64,
    t0: f64,
    t_end: f64,
    h: f64,
    e: f64,
) -> Vec<(f64, f64, f64)>
where
    F: Fn(f64, f64, f64) -> (f64, f64),
{
    let mut result = Vec::new();

    let mut rk_result = runge_kutta4(&f, x0, z0, t0, t0 + h, h, e);
    result.extend(rk_result);

    let mut t = t0 + h;
    let mut x = result[result.len() - 1].1;
    let mut z = result[result.len() - 1].2;

    while t <= t_end {
        let (fx, fz) = f(x, z, t);

        if result.len() >= 2 {
            let prev_f1 = f(result[result.len() - 1].1, result[result.len() - 1].2, result[result.len() - 1].0);
            let prev_f2 = f(result[result.len() - 2].1, result[result.len() - 2].2, result[result.len() - 2].0);

            x += (h / 2.0) * (3.0 * fx - prev_f2.0);
            z += (h / 2.0) * (3.0 * fz - prev_f2.1);
        }

        t += h;
        result.push((t, x, z));
    }

    result
}

pub fn adams3<F>(
    f: F,
    x0: f64,
    z0: f64,
    t0: f64,
    t_end: f64,
    h: f64,
    e: f64,
) -> Vec<(f64, f64, f64)>
where
    F: Fn(f64, f64, f64) -> (f64, f64),
{
    let mut result = Vec::new();

    let mut rk_result = runge_kutta4(&f, x0, z0, t0, t0 + 2.0 * h, h, e);
    result.extend(rk_result);

    let mut t = t0 + 2.0 * h;
    let mut x = result[result.len() - 1].1;
    let mut z = result[result.len() - 1].2;

    while t <= t_end {
        let (fx, fz) = f(x, z, t);

        if result.len() >= 3 {
            let prev_f1 = f(result[result.len() - 1].1, result[result.len() - 1].2, result[result.len() - 1].0);
            let prev_f2 = f(result[result.len() - 2].1, result[result.len() - 2].2, result[result.len() - 2].0);
            let prev_f3 = f(result[result.len() - 3].1, result[result.len() - 3].2, result[result.len() - 3].0);

            x += (h / 12.0) * (23.0 * fx - 16.0 * prev_f2.0 + 5.0 * prev_f3.0);
            z += (h / 12.0) * (23.0 * fz - 16.0 * prev_f2.1 + 5.0 * prev_f3.1);
        }

        t += h;
        result.push((t, x, z));
    }

    result
}

pub fn adams4<F>(
    f: F,
    x0: f64,
    z0: f64,
    t0: f64,
    t_end: f64,
    h: f64,
    e: f64,
) -> Vec<(f64, f64, f64)>
where
    F: Fn(f64, f64, f64) -> (f64, f64),
{
    let mut result = Vec::new();

    let mut rk_result = runge_kutta4(&f, x0, z0, t0, t0 + 3.0 * h, h, e);
    result.extend(rk_result);

    let mut t = t0 + 3.0 * h;
    let mut x = result[result.len() - 1].1;
    let mut z = result[result.len() - 1].2;

    while t <= t_end {
        let (fx, fz) = f(x, z, t);

        if result.len() >= 4 {
            let prev_f1 = f(result[result.len() - 1].1, result[result.len() - 1].2, result[result.len() - 1].0);
            let prev_f2 = f(result[result.len() - 2].1, result[result.len() - 2].2, result[result.len() - 2].0);
            let prev_f3 = f(result[result.len() - 3].1, result[result.len() - 3].2, result[result.len() - 3].0);
            let prev_f4 = f(result[result.len() - 4].1, result[result.len() - 4].2, result[result.len() - 4].0);

            x += (h / 24.0) * (55.0 * fx - 59.0 * prev_f2.0 + 37.0 * prev_f3.0 - 9.0 * prev_f4.0);
            z += (h / 24.0) * (55.0 * fz - 59.0 * prev_f2.1 + 37.0 * prev_f3.1 - 9.0 * prev_f4.1);
        }

        t += h;
        result.push((t, x, z));
    }

    result
}