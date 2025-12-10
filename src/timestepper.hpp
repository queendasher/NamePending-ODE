#ifndef TIMERSTEPPER_HPP
#define TIMERSTEPPER_HPP

#include <functional>
#include <exception>

#include "Newton.hpp"


namespace ASC_ode
{
  
  class TimeStepper
  { 
  protected:
    std::shared_ptr<NonlinearFunction> m_rhs;
  public:
    TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
    virtual ~TimeStepper() = default;
    virtual void doStep(double tau, VectorView<double> y) = 0;
  };

  class ExplicitEuler : public TimeStepper
  {
    Vector<> m_vecf;
  public:
    ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}
    void doStep(double tau, VectorView<double> y) override
    {
      this->m_rhs->evaluate(y, m_vecf);
      y += tau * m_vecf;
    }
  };

  class ImplicitEuler : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
  public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) 
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      m_equ = ynew - m_yold - m_tau * m_rhs;
    }

    void doStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau->set(tau);
      NewtonSolver(m_equ, y);
    }
  };

  class ImprovedEuler : public TimeStepper 
  {
  public:
    explicit ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs_in)
      : TimeStepper(std::move(rhs_in)) {}

    void doStep(double tau, VectorView<double> y) override
    {
      // f(y_n)
      Vector<> fyn(m_rhs->dimF());
      m_rhs->evaluate(y, fyn);

      // y~ = y_n + (tau/2) * f(y_n)
      Vector<> ytilde(m_rhs->dimX());
      ytilde = y;
      for (size_t i = 0; i < ytilde.size(); ++i)
        ytilde(i) += 0.5 * tau * fyn(i);

      // f(y~)
      Vector<> fytilde(m_rhs->dimF());
      m_rhs->evaluate(ytilde, fytilde);

      // y_{n+1} = y_n + tau * f(y~)
      for (size_t i = 0; i < y.size(); ++i)
        y(i) += tau * fytilde(i);
    }
  };

  class ExplicitRungeKutta : public TimeStepper 
  {
    public:
    explicit ExplicitRungeKutta(std::shared_ptr<NonlinearFunction> rhs_in)
      : TimeStepper(std::move(rhs_in)) {}
    void doStep(double tau, VectorView<double> y) override
    {
      size_t n = m_rhs->dimX();
      Vector<> k1(n), k2(n), k3(n), k4(n);
      Vector<> ytemp(n);
      // k1 = f(y_n)
      m_rhs->evaluate(y, k1);

      // k2 = f(y_n + tau/2 * k1)
      for (size_t i = 0; i < n; ++i)
        ytemp(i) = y(i) + 0.5 * tau * k1(i);
      m_rhs->evaluate(ytemp, k2);

      // k3 = f(y_n + tau/2 * k2)
      for (size_t i = 0; i < n; ++i)
        ytemp(i) = y(i) + 0.5 * tau * k2(i);
      m_rhs->evaluate(ytemp, k3);

      // k4 = f(y_n + tau * k3)
      for (size_t i = 0; i < n; ++i)
        ytemp(i) = y(i) + tau * k3(i);
      m_rhs->evaluate(ytemp, k4);


      // y_{n+1} = y_n + tau/6 * (k1 + 2*k2 + 2*k3 + k4)
      for (size_t i = 0; i < n; ++i)
        y(i) += (tau / 6.0) * (k1(i) + 2.0 * k2(i) + 2.0 * k3(i) + k4(i));
    }
  };


  class CrankNicolson : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau_half;
    std::shared_ptr<ConstantFunction> m_yold;
    std::shared_ptr<ConstantFunction> m_rhs_2;

  public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_tau_half(std::make_shared<Parameter>(0.0)) 
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      m_rhs_2 = std::make_shared<ConstantFunction>(m_rhs->dimF());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      m_equ = ynew - m_yold - m_tau_half * (m_rhs_2 + m_rhs);
    }

    void doStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau_half->set(tau * 0.5);

      Vector<> fyn(m_rhs->dimF());
      m_rhs->evaluate(y, fyn);
      m_rhs_2->set(fyn);

      NewtonSolver(m_equ, y);
    }
  };  
}


#endif
