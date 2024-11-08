mod parts;
use parts::nonlinear_equation;
use parts::system_of_nonlinear_equations;

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
    if let Some(root) = nonlinear_equation::mpi_1(x0, tolerance, max_iter) {
        println!("Корень (МПИ, x0 = {}): x = {}", x0, root);
    } else {
        println!("МПИ не сошелся для начального приближения x0 = {}", x0);
    }

    x0 = 4.5;
    if let Some(root) = nonlinear_equation::mpi_2(x0, tolerance, max_iter) {
        println!("Корень (МПИ, x0 = {}): x = {}", x0, root);
    } else {
        println!("МПИ не сошелся для начального приближения x0 = {}", x0);
    }
    

    for x0 in [0.4, 4.5].iter().cloned() {
        if let Some(root) = nonlinear_equation::newton(x0, tolerance, max_iter, h) {
            println!("Корень (Ньютон, x0 = {}): x = {}", x0, root);
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

    if let Some((x, y)) = system_of_nonlinear_equations::mpi_system(1.0, 1.0, tolerance, max_iter) {
        println!("Корни системы (МПИ): x = {:.6}, y = {:.6}", x, y);
    } else {
        println!("МПИ для системы не сошелся");
    }

    if let Some((x, y)) = system_of_nonlinear_equations::newton_system(1.0, 1.0, tolerance, max_iter, h) {
        println!("Корни системы (Ньютон): x = {:.6}, y = {:.6}", x, y);
    } else {
        println!("Метод Ньютона для системы не сошелся");
    }





}