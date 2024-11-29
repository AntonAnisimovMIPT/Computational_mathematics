mod parts;
use parts::nonlinear_equation;
use parts::system_of_nonlinear_equations;
use plotters::prelude::*;

fn plot_convergence(data: Vec<f64>, method_name: &str, output_file: &str) {
    let root = BitMapBackend::new(output_file, (800, 600)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    let max_value = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut chart = ChartBuilder::on(&root)
        .caption(format!("График сходимости метода {method_name}"), ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(0..data.len(), 0.0..max_value)
        .unwrap();

    chart
        .configure_mesh()
        .x_desc("Номер итерации")
        .y_desc("Значение ошибки")
        .draw()
        .unwrap();

    chart
        .draw_series(LineSeries::new(
            data.iter().enumerate().map(|(i, &val)| (i, val)),
            &BLUE,
        ))
        .unwrap()
        .label("Сходимость")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &BLUE));

    chart.configure_series_labels().border_style(&BLACK).draw().unwrap();
}

fn main() {
    nonlinear_equation::print_func();
    let tolerance = 1e-5;
    let max_iter = 10_usize.pow(3);
    let h = 1e-5;

    println!(
    "Заданные параметры:
        точность: {tolerance},
        максимальное число итераций: {max_iter},
        h: {h}");

    let mut x0 =  0.4;
    if let Some((root, errors)) = nonlinear_equation::mpi_1_with_errors(x0, tolerance, max_iter) {
        println!("Корень (МПИ, x0 = {}): x = {}", x0, root);
        plot_convergence(
            errors, 
            &format!("МПИ_x0={}", x0), 
            &format!("./графики_сходимости/part1/MPI_x0={}_convergence.png", x0));
    } else {
        println!("МПИ не сошелся для начального приближения x0 = {}", x0);
    }

    x0 = 4.5;
    if let Some((root, errors)) = nonlinear_equation::mpi_2_with_errors(x0, tolerance, max_iter) {
        println!("Корень (МПИ, x0 = {}): x = {}", x0, root);
        plot_convergence(
            errors, 
            &format!("МПИ_x0={}", x0), 
            &format!("./графики_сходимости/part1/MPI_x0={}_convergence.png", x0));
    } else {
        println!("МПИ не сошелся для начального приближения x0 = {}", x0);
    }
    

    for x0 in [0.4, 4.5].iter().cloned() {
        if let Some((root, errors)) = nonlinear_equation::newton_with_errors(x0, tolerance, max_iter, h) {
            println!("Корень (Ньютон, x0 = {}): x = {}", x0, root);
            plot_convergence(
                errors,
                &format!("Ньютон x0 = {}", x0), 
                &format!("./графики_сходимости/part1/Newton_x0_{}_convergence.png", x0) );
        } else {
            println!("Метод Ньютона не сошелся для начального приближения x0 = {}", x0);
        }
    }

    println!("-------------------------------------------------");
        
    system_of_nonlinear_equations::print_sistem_func();

    println!(
    "Заданные параметры:
        точность: {tolerance},
        максимальное число итераций: {max_iter},
        h: {h}");

    if let Some((x, y, errors)) = system_of_nonlinear_equations::mpi_system(1.0, 1.0, tolerance, max_iter) {
        println!("Корни системы (МПИ): x = {:.6}, y = {:.6}", x, y);
        plot_convergence(
            errors, 
            "Сходимость метода МПИ для системы", 
            "./графики_сходимости/part2/mpi_system_convergence.png"
        );
    } else {
        println!("МПИ для системы не сошелся");
    }

    if let Some((x, y, errors)) = system_of_nonlinear_equations::newton_system(1.0, 1.0, tolerance, max_iter, h) {
        println!("Корни системы (Ньютон): x = {:.6}, y = {:.6}", x, y);
        plot_convergence(
            errors, 
            "Сходимость метода Ньютона для системы", 
            "./графики_сходимости/part2/newton_system_convergence.png"
        );
    } else {
        println!("Метод Ньютона для системы не сошелся");
    }


}