#include <iostream>
#include <autodiff.hpp>


using namespace ASC_ode;


template <typename T>
T func1 (T x, T y)
{
  return x * sin(y);
  // return 1e6 + y;
}



int main()
{
  double x = 1, y = 2;
  AutoDiff<> adx(x, 0, 2);
  AutoDiff<> ady(y, 1, 2);

  std::cout << "adx = " << adx << std::endl;
  std::cout << "ady = " << ady << std::endl;

  AutoDiff<> prod = adx * ady;
  std::cout << "prod = " << prod << std::endl;

  std::cout << "func1(adx, ady) = " << func1(adx, ady) << std::endl;

  double eps = 1e-8;
  std::cout << "numdiff df/dx = " << (func1(x + eps, y) - func1(x-eps, y)) / (2*eps) << std::endl;
  std::cout << "numdiff df/dy = " << (func1(x, y + eps) - func1(x, y-eps)) / (2*eps) << std::endl;


  {
    // we can do second derivatives:
    AutoDiff<> addx(2, 0, 1);
    std::cout << "addx = " << addx << std::endl;
    // func = x*x
    // func' = 2*x
    // func'' = 2
    std::cout << "addx*addx = " << addx * addx << std::endl;

    // std::cout << "sin(addx) = " << sin(addx) << std::endl;
  }

  // Evaluate and plot Legendre polynomials
  {
    static constexpr double STEP = 0.1;
    std::vector<AutoDiff<>> P;
    for (double x = -1.0; x <= 1.0; x += STEP) {
        AutoDiff<> dx(x, 0, 1);
        LegendrePolynomials(5, dx, P);
        std::cout << "x = " << x << ": ";
        for (size_t n = 0; n < P.size(); ++n) {
            std::cout << "P_" << n << "(x) = " << P[n] << " ";
        }
        std::cout << std::endl;
    }
  }
  return 0;
}