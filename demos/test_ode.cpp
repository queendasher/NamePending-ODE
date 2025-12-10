#include <iostream>
#include <fstream> 

#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <implicitRK.hpp>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
private:
  double mass;
  double stiffness;

public:
  MassSpring(double m, double k) : mass(m), stiffness(k) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -stiffness/mass*x(0);
  }
  
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -stiffness/mass;
  }
};

class RCCircuit : public NonlinearFunction
{
private:
  double R, C;

public:
  RCCircuit(double R_, double C_) : R(R_), C(C_) {}

  size_t dimX() const override { return 2; }   // (U_C, t)
  size_t dimF() const override { return 2; }

  void evaluate(const VectorView<double> x, VectorView<double> f) const override
  {
    double Uc = x(0);   // capacitor voltage U_C
    double t  = x(1);   // time

    double U0 = std::cos(100.0 * M_PI * t);

    f(0) = (U0 - Uc) / (R * C);    // dU_C/dt
    f(1) = 1.0;                    // dt/dt
  }

  void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override
  {
    const double t = x(1);
    const double omega = 100.0 * M_PI;
    const double invRC = 1.0 / (R * C);

    df = 0.0; // zero everything first
    df(0,0) = -invRC;             
    df(0,1) = -(omega * std::sin(omega * t)) * invRC; 
  }
};


int main()
{
  double R = 1.0;
  double C = 0.01;

  double tend = 0.1;
  int steps = 250;
  double tau = tend/steps;

<<<<<<< HEAD
  Vector<> y = { 0, 0 };  // initializer list
  auto rhs = std::make_shared<RCCircuit>(R, C);
  
  ImplicitEuler stepper(rhs);
=======
  Vector<> y = { 1, 0 };  // initializer list
  auto rhs = std::make_shared<MassSpring>(1.0, 1.0);



/*
  Vector<> Radau(3), RadauWeight(3);
  GaussRadau (Radau, RadauWeight);
  // not sure about weights, comput them via ComputeABfromC
  cout << "Radau = " << Radau << ", weight = " << RadauWeight <<  endl;
        Vector<> Gauss2c(2), Gauss3c(3);
*/
 

  // ExplicitEuler stepper(rhs);
  // ImplicitEuler stepper(rhs);
>>>>>>> upstream/RungeKutta

  // RungeKutta stepper(rhs, Gauss2a, Gauss2b, Gauss2c);

  // Gauss3c .. points tabulated, compute a,b:
  auto [Gauss3a,Gauss3b] = computeABfromC (Gauss3c);
  ImplicitRungeKutta stepper(rhs, Gauss3a, Gauss3b, Gauss3c);


  /*
  // arbitrary order Gauss-Legendre
  int stages = 5;
  Vector<> c(stages), b1(stages);
  GaussLegendre(c, b1);

  auto [a, b] = computeABfromC(c);
  ImplicitRungeKutta stepper(rhs, a, b, c);
  */

  /* 
  // arbitrary order Radau
  int stages = 5;
  Vector<> c(stages), b1(stages);
  GaussRadau(c, b1);

  auto [a, b] = computeABfromC(c);
  ImplicitRungeKutta stepper(rhs, a, b, c);
  */


  std::ofstream outfile ("output_test_ode.txt");
  std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
  outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     stepper.doStep(tau, y);

     std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
     outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
  }
}
