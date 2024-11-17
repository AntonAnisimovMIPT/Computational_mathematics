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

    println!("Классическая полиномиальная интерполяция по Ньютону:");

    struct Result {
        window_size: usize,
        start_year: f64,
        end_year: f64,
        estimated_population: f64,
        deviation: f64,
    }
    let mut results: Vec<Result> = Vec::new();

    let min_window_size = 2;
    let max_window_size = years.len();

    for window_size in min_window_size..=max_window_size {
        println!("Размер окна: {}", window_size);

        for start_index in 0..=years.len() - window_size {
            let selected_years = &years[start_index..start_index + window_size];
            let selected_populations = &populations[start_index..start_index + window_size];

            let diff_table = classical_newton_polynomial_interpolation::divided_differences(selected_years, selected_populations);
            let estimated_population = classical_newton_polynomial_interpolation::newton_polynomial(selected_years, target_year, &diff_table);
            let deviation = (estimated_population - EXACT_POPULATION_VALUE).abs() * 100.0 / EXACT_POPULATION_VALUE;

            println!(
                "Сдвиг окна: от {} до {}. Экстраполированное население в {} году: {:.0}, отклонение: {:.2}%",
                selected_years[0], selected_years[window_size - 1], target_year, estimated_population, deviation
            );

            results.push(Result {
                window_size,
                start_year: selected_years[0],
                end_year: selected_years[window_size - 1],
                estimated_population,
                deviation,
            });
        }
        println!("...........");
    }

    results.sort_by(|a, b| a.deviation.partial_cmp(&b.deviation).unwrap());

    println!("Топ 10 вариаций с наименьшим отклонением:");
    for (i, result) in results.iter().take(10).enumerate() {
        println!(
            "{:<3} | Размер окна: {:<2} | Сдвиг окна: {} - {} | Экстраполированное население: {:.0} | Отклонение: {:.2}%",
            i + 1,
            result.window_size,
            result.start_year,
            result.end_year,
            result.estimated_population,
            result.deviation
        );
    }

    let spline = CubicSpline::new(&years, &populations);
    let estimated_population = spline.evaluate(target_year);
    let deviation = (estimated_population - EXACT_POPULATION_VALUE).abs() * 100.0 / EXACT_POPULATION_VALUE;

    println!("\n\n\nСплайн-интерполяция:");
    println!(
        "   Экстраполированное население США в {} году: {:.0}, относительное отклонение от истинного значения: {:.2}%",
        target_year, estimated_population, deviation
    );

    let degree = 2;
    let polynomial = LeastSquaresPolynomial::new(&years, &populations, degree);
    let estimated_population = polynomial.evaluate(target_year);
    let deviation = (estimated_population - EXACT_POPULATION_VALUE).abs() * 100.0 / EXACT_POPULATION_VALUE;

    println!("\n\n\nМетод наименьших квадратов");
    println!(
        "   Экстраполированное население США в {} году: {:.0}, относительное отклонение от истинного значения: {:.2}%",
        target_year, estimated_population, deviation
    );
}
