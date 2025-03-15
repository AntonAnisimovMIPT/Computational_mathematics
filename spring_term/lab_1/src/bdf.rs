use crate::runge_kutta4;

pub fn bdf2<F>(
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
            let prev_x = result[result.len() - 2].1;
            let prev_z = result[result.len() - 2].2;

            x = (4.0 / 3.0) * x - (1.0 / 3.0) * prev_x + (2.0 / 3.0) * h * fx;
            z = (4.0 / 3.0) * z - (1.0 / 3.0) * prev_z + (2.0 / 3.0) * h * fz;
        }

        t += h;
        result.push((t, x, z));
    }

    result
}

pub fn bdf3<F>(
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
            let prev_x1 = result[result.len() - 1].1;
            let prev_x2 = result[result.len() - 2].1;
            let prev_x3 = result[result.len() - 3].1;

            let prev_z1 = result[result.len() - 1].2;
            let prev_z2 = result[result.len() - 2].2;
            let prev_z3 = result[result.len() - 3].2;

            x = (18.0 / 11.0) * prev_x1 - (9.0 / 11.0) * prev_x2 + (2.0 / 11.0) * prev_x3 + (6.0 / 11.0) * h * fx;
            z = (18.0 / 11.0) * prev_z1 - (9.0 / 11.0) * prev_z2 + (2.0 / 11.0) * prev_z3 + (6.0 / 11.0) * h * fz;
        }

        t += h;
        result.push((t, x, z));
    }

    result
}

pub fn bdf4<F>(
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
            let prev_x1 = result[result.len() - 1].1;
            let prev_x2 = result[result.len() - 2].1;
            let prev_x3 = result[result.len() - 3].1;
            let prev_x4 = result[result.len() - 4].1;

            let prev_z1 = result[result.len() - 1].2;
            let prev_z2 = result[result.len() - 2].2;
            let prev_z3 = result[result.len() - 3].2;
            let prev_z4 = result[result.len() - 4].2;

            x = (48.0 / 25.0) * prev_x1 - (36.0 / 25.0) * prev_x2 + (16.0 / 25.0) * prev_x3 - (3.0 / 25.0) * prev_x4 + (12.0 / 25.0) * h * fx;
            z = (48.0 / 25.0) * prev_z1 - (36.0 / 25.0) * prev_z2 + (16.0 / 25.0) * prev_z3 - (3.0 / 25.0) * prev_z4 + (12.0 / 25.0) * h * fz;
        }

        t += h;
        result.push((t, x, z));
    }

    result
}