mod parts;
use parts::classical_newton_polynomial_interpolation;
use parts::spline_interpolation::CubicSpline;
use parts::least_squares::LeastSquaresPolynomial;

fn main() {

    let target_year = 2010.0;
    const EXACT_POPULATION_VALUE: f64 = 308_745_538.0;
    println!("Истинное население США в 2010 году: {}", EXACT_POPULATION_VALUE);
    println!("___________");

    let years = vec![
        1910.0, 
        1920.0, 
        1930.0, 
        1940.0, 
        1950.0,
        1960.0, 
        1970.0, 
        1980.0, 
        1990.0, 
        2000.0,
    ];

    let populations = vec![
        92_228_496.0,
        106_021_537.0,
        123_202_624.0,
        132_164_569.0,
        151_325_798.0,
        179_323_175.0,
        203_211_926.0,
        226_545_805.0,
        248_709_873.0,
        281_421_906.0,
    ];

    let diff_table = classical_newton_polynomial_interpolation::divided_differences(&years, &populations);
    let estimated_population = classical_newton_polynomial_interpolation::newton_polynomial(&years, target_year, &diff_table);
    let deviation = (estimated_population - EXACT_POPULATION_VALUE).abs() * 100.0 / EXACT_POPULATION_VALUE;

    println!("Классическая полиномиальная интерполяция по Ньютону:");
    println!(
        "   Экстраполированное население США в {} году: {:.0}, относительное отклонение от истинного значения: {:.2}%",
        target_year, estimated_population, deviation
    );

    let spline = CubicSpline::new(&years, &populations);
    let estimated_population = spline.evaluate(target_year);
    let deviation = (estimated_population - EXACT_POPULATION_VALUE).abs() * 100.0 / EXACT_POPULATION_VALUE;

    println!("Сплайн-интерполяция:");
    println!(
        "   Экстраполированное население США в {} году: {:.0}, относительное отклонение от истинного значения: {:.2}%",
        target_year, estimated_population, deviation
    );

    let degree = 2;
    let polynomial = LeastSquaresPolynomial::new(&years, &populations, degree);
    let estimated_population = polynomial.evaluate(target_year);
    let deviation = (estimated_population - EXACT_POPULATION_VALUE).abs() * 100.0 / EXACT_POPULATION_VALUE;

    println!("Метод наименьших квадратов");
    println!(
        "   Экстраполированное население США в {} году: {:.0}, относительное отклонение от истинного значения: {:.2}%",
        target_year, estimated_population, deviation
    );
}
