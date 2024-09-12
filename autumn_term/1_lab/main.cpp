#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <filesystem>

auto f1(double x) {
    return sin(x * x);
}
auto derivative_of_f1(double x) {
    return 2 * x * cos(x * x);
}

auto f2(double x) {
    return cos(sin(x));
}
auto derivative_of_f2(double x) {
    return -sin(sin(x)) * cos(x);
}

auto f3(double x) {
    return exp(sin(cos(x)));
}
auto derivative_of_f3(double x) {
    return -exp(sin(cos(x))) * sin(x) * cos(cos(x));
}

auto f4(double x) {
    return log(x + 3);
}
auto derivative_of_f4(double x) {
    return 1 / (x + 3);
}

auto f5(double x) {
    return sqrt(x + 3);
}
auto derivative_of_f5(double x) {
    return 0.5 / sqrt(x + 3);
}

auto first_method(
    double (*f)(double),
    double x,
    double h) {

    return (f(x + h) - f(x)) / h;
}

auto second_method(
    double (*f)(double),
    double x,
    double h) {

    return (f(x) - f(x - h)) / h;
}

auto third_method(
    double (*f)(double),
    double x,
    double h) {

    return (f(x + h) - f(x - h)) / (2 * h);
}

auto fourth_method(
    double (*f)(double),
    double x,
    double h) {

    return (4.0 / 3.0) * (f(x + h) - f(x - h)) / (2 * h) - (1.0 / 3.0) * (f(x + 2 * h) - f(x - 2 * h)) / (4 * h);
}

auto fifth_method(
    double (*f)(double),
    double x,
    double h) {

    return (3.0 / 2.0) * (f(x + h) - f(x - h)) / (2 * h) - (3.0 / 5.0) * (f(x + 2 * h) - f(x - 2 * h)) / (4 * h) + (1.0 / 10.0) * (f(x + 3 * h) - f(x - 3 * h)) / (6 * h);
}

auto abs_error(double exact, double approx) {
    return std::abs(exact - approx);
}

void save_data(
    const std::string& filename,
    const std::vector<double>& h_values,
    const std::vector<double>& errors) {

    std::ofstream file(filename);
    file << "h,error\n";
    for (size_t i = 0; i < h_values.size(); ++i) {
        file << h_values[i] << "," << errors[i] << "\n";
    }
    file.close();
}

int main() {

    auto x0 = 10.0;
    auto n_steps = 21;
    std::vector<double> h_values(n_steps);

    auto methods = {first_method, second_method, third_method, fourth_method, fifth_method};
    std::vector<std::string> method_names = {"method1.csv", "method2.csv", "method3.csv", "method4.csv", "method5.csv"};

    auto functions = {f1, f2, f3, f4, f5};
    std::vector<std::string> function_names = {"f1", "f2", "f3", "f4", "f5"};
    auto derivatives = {derivative_of_f1, derivative_of_f2, derivative_of_f3, derivative_of_f4, derivative_of_f5};

    for (size_t f_index = 0; f_index < functions.size(); ++f_index) {

        std::filesystem::create_directory("../calculated_data/" + function_names[f_index]);

        auto f = *std::next(functions.begin(), f_index);
        auto exact_derivative = *std::next(derivatives.begin(), f_index);

        for (size_t m_index = 0; m_index < methods.size(); ++m_index) {
            auto method = *std::next(methods.begin(), m_index);

            std::vector<double> errors(n_steps);

            for (int n = 0; n < n_steps; ++n) {
                auto h = 2.0 / pow(2, n + 1);
                h_values[n] = h;

                auto exact = exact_derivative(x0);
                auto approx = method(f, x0, h);
                errors[n] = abs_error(exact, approx);
            }

            std::string filename = "../calculated_data/" + function_names[f_index] + "/" + method_names[m_index];
            save_data(filename, h_values, errors);
        }
    }

    std::cout << "Data saved in calculated_data/ directory." << std::endl;

    return 0;
}
