#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

const double X0 = 0.0, X1 = 1.0;
const double Y0 = 0.0, Y1 = 2.0;
const double EPS = 1e-6;

void rk4_step(std::function<void(double, const std::vector<double> &,
                                 std::vector<double> &)>
                  f,
              double &x, std::vector<double> &u, double h) {
  int n = u.size();
  std::vector<double> k1(n), k2(n), k3(n), k4(n), temp(n);

  f(x, u, k1);
  for (int i = 0; i < n; ++i)
    temp[i] = u[i] + 0.5 * h * k1[i];
  f(x + 0.5 * h, temp, k2);
  for (int i = 0; i < n; ++i)
    temp[i] = u[i] + 0.5 * h * k2[i];
  f(x + 0.5 * h, temp, k3);
  for (int i = 0; i < n; ++i)
    temp[i] = u[i] + h * k3[i];
  f(x + h, temp, k4);

  for (int i = 0; i < n; ++i)
    u[i] += (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
  x += h;
}

// Функция системы ОДУ для метода стрельбы
/*
ввел замену:
u_1 = y,
u_2 = y'
тогда исходное уравнение можно свести к виду:
u_1' = u_2,
u_2' = x*sqrt(u_1)
*/
void system_ODE(double x, const std::vector<double> &u,
                std::vector<double> &du) {
  du[0] = u[1];                // u_1' = u_2
  du[1] = x * std::sqrt(u[0]); // u_2' = x * sqrt(u_1)
}

double solve_cauchy(double alpha, int N, std::vector<double> &x,
                    std::vector<double> &y) {
  double h = (X1 - X0) / N;
  x.resize(N + 1);
  y.resize(N + 1);

  std::vector<double> u = {Y0, alpha};
  x[0] = X0;
  y[0] = Y0;

  for (int i = 1; i <= N; ++i) {
    rk4_step(system_ODE, x[i - 1], u, h);
    x[i] = X0 + i * h;
    y[i] = u[0];
  }
  return y[N];
}

void shooting_method(int N, const std::string &filename) {
  double alpha1 = 1.0, alpha2 = 3.0;
  std::vector<double> x, y;

  double y1 = solve_cauchy(alpha1, N, x, y) - Y1;
  double y2 = solve_cauchy(alpha2, N, x, y) - Y1;

  while (std::abs(y2) > EPS) {
    double alpha_new = alpha2 - y2 * (alpha2 - alpha1) / (y2 - y1);
    alpha1 = alpha2;
    y1 = y2;
    alpha2 = alpha_new;
    y2 = solve_cauchy(alpha2, N, x, y) - Y1;
  }

  std::ofstream out(filename);
  for (int i = 0; i <= N; ++i)
    out << x[i] << " " << y[i] << "\n";
  out.close();
}

std::vector<double> tridiagonal_solve(const std::vector<double> &a,
                                      const std::vector<double> &b,
                                      const std::vector<double> &c,
                                      const std::vector<double> &d) {
  int n = b.size();
  std::vector<double> v(n), p(n - 1), q(n - 1);

  p[0] = c[0] / b[0];
  q[0] = d[0] / b[0];

  for (int i = 1; i < n - 1; ++i) {
    double denom = b[i] - a[i - 1] * p[i - 1];
    p[i] = c[i] / denom;
    q[i] = (d[i] - a[i - 1] * q[i - 1]) / denom;
  }

  v[n - 1] =
      (d[n - 1] - a[n - 2] * q[n - 2]) / (b[n - 1] - a[n - 2] * p[n - 2]);
  for (int i = n - 2; i >= 0; --i)
    v[i] = p[i] * v[i + 1] + q[i];

  return v;
}

/*
замена y_n+1 = y_n + v
тогда исходное уравнение становится (y_n+1)'' = x*sqrt(y_n+v)
по Тейлору sqrt(y_n+v) = sqrt(y_n) + v/(2*sqrt(y_n))
тогда получаем v'' = x*v/(2*sqrt(y_n))
*/
void quasilinearization_method(int N, const std::string &filename) {
  double h = (X1 - X0) / N;
  std::vector<double> x(N + 1), y(N + 1), y_new(N + 1);

  // Начальное приближение: y_0(x) = 2x
  for (int i = 0; i <= N; ++i) {
    x[i] = X0 + i * h;
    y[i] = 2.0 * x[i];
  }

  double error;
  do {
    // Решаем для v
    std::vector<double> a(N - 1), b(N), c(N - 1), d(N);

    // Граничные условия
    b[0] = 1.0; // v_0 = 0
    d[0] = 0.0;
    b[N - 1] = 1.0; // v_N = 2 - y_N
    d[N - 1] = Y1 - y[N];

    for (int i = 1; i < N - 1; ++i) {
      double coef = x[i] / (2.0 * std::sqrt(y[i] + 1e-10));
      a[i - 1] = 1.0 / (h * h);
      b[i] = -2.0 / (h * h) - coef;
      c[i] = 1.0 / (h * h);
      d[i] = 0.0;
    }

    std::vector<double> v = tridiagonal_solve(a, b, c, d);

    // Обновляем y
    error = 0.0;
    y_new[0] = y[0];
    y_new[N] = Y1;
    for (int i = 1; i < N; ++i) {
      y_new[i] = y[i] + v[i];
      error = std::max(error, std::abs(y_new[i] - y[i]));
    }
    y = y_new;
  } while (error > EPS);

  std::ofstream out(filename);
  for (int i = 0; i <= N; ++i)
    out << x[i] << " " << y[i] << "\n";
  out.close();
}

void check_convergence() {
  std::vector<int> Ns = {100, 200, 400, 500, 750, 800, 10000};
  std::vector<std::vector<double>> solutions;

  for (int N : Ns) {
    std::string filename = "data/shooting_" + std::to_string(N) + ".txt";
    shooting_method(N, filename);

    std::vector<double> x(N + 1), y(N + 1);
    std::ifstream in(filename);
    for (int i = 0; i <= N; ++i)
      in >> x[i] >> y[i];
    in.close();
    solutions.push_back(y);
  }

  // Вычисляем ошибки
  for (size_t k = 1; k < Ns.size(); ++k) {
    int N_coarse = Ns[k - 1];
    int N_fine = Ns[k];
    double error = 0.0;

    // Сравниваем на узлах более грубой сетки
    for (int i = 0; i <= N_coarse; ++i) {
      int idx_fine = i * (N_fine / N_coarse);
      error = std::max(error,
                       std::abs(solutions[k - 1][i] - solutions[k][idx_fine]));
    }
    std::cout << "Ошибка между N=" << N_coarse << " и N=" << N_fine << ": "
              << error << "\n";
  }
}

int main() {

  shooting_method(100, "data/shooting_100.txt");
  shooting_method(200, "data/shooting_200.txt");
  shooting_method(400, "data/shooting_400.txt");
  shooting_method(500, "data/shooting_500.txt");
  shooting_method(750, "data/shooting_750.txt");
  shooting_method(800, "data/shooting_800.txt");
  shooting_method(1000, "data/shooting_1000.txt");

  quasilinearization_method(100, "data/quasi_100.txt");
  quasilinearization_method(200, "data/quasi_200.txt");
  quasilinearization_method(400, "data/quasi_400.txt");
  quasilinearization_method(500, "data/quasi_500.txt");
  quasilinearization_method(750, "data/quasi_750.txt");
  quasilinearization_method(800, "data/quasi_800.txt");
  quasilinearization_method(1000, "data/quasi_1000.txt");

  check_convergence();

  return 0;
}