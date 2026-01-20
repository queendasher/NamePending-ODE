#ifndef AUTODIFF_HPP
#define AUTODIFF_HPP

#include <cstddef> 
#include <ostream> 
#include <cmath>   
#include <array>  
#include <vector>


namespace ASC_ode
{

  template <size_t N, typename T = double>
  class Variable 
  {
    private:
      T m_val;
    public:
      Variable (T v) : m_val(v) {}
      T value() const { return m_val; }
  };

  template <typename T = double>
  auto derivative (T v, size_t /*index*/) { return T(0); } 


  template <typename T = double>
  class AutoDiff
  {
  private:
    T m_val;
    std::vector<T> m_deriv;
  public: 
    AutoDiff () : m_val(0), m_deriv() {}
    AutoDiff (T v, size_t n=0) : m_val(v), m_deriv(n, T(0)) {}

    AutoDiff(T v, size_t id, size_t n) : m_val(v), m_deriv(n, T(0))
    {
        if (id < n)
            m_deriv[id] = T(1);
    }

    T value() const { return m_val; }
    std::vector<T>& deriv() { return m_deriv; }
    const std::vector<T>& deriv() const { return m_deriv; }
    size_t size() const { return m_deriv.size(); }
    void resize(size_t n) { m_deriv.resize(n, T(0)); }
  };


  template <typename T = double>
  auto derivative (AutoDiff<T> v, size_t index) 
  {
    return (index < v.deriv().size()) ? v.deriv()[index] : T(0);
  }

  template <typename T>
  std::ostream & operator<< (std::ostream& os, const AutoDiff<T>& ad)
  {
    os << "Value: " << ad.value() << ", Deriv: [";
    for (size_t i = 0; i < ad.size(); ++i)
    {
      os << ad.deriv()[i];
      if (i < ad.size() - 1) os << ", ";
    }
    os << "]";
    return os;
  }

  template <typename T>
  size_t maxSize(const AutoDiff<T>& a, const AutoDiff<T>& b)
  {
      return std::max(a.size(), b.size());
  }

  template <typename T = double>
    AutoDiff<T> operator+ (const AutoDiff<T>& a, const AutoDiff<T>& b)
    {
    size_t n = maxSize(a, b);
        AutoDiff<T> result(a.value() + b.value());
        for (size_t i = 0; i < n; ++i)
        result.deriv()[i] = derivative(a, i) + derivative(b, i);
        return result;
    }

    template <typename T = double>
    auto operator+ (T a, const AutoDiff<T>& b) { return AutoDiff<T>(a, b.size()) + b; }

    template <typename T = double>
    auto operator+ (const AutoDiff<T>& a, T b) { return a + AutoDiff<T>(b, a.size()); }

    template <typename T = double>
    auto operator+= (AutoDiff<T>& a, const AutoDiff<T>& b)
    {
        a = a + b;
        return a;
    }

    template <typename T = double>
    auto operator- (const AutoDiff<T>& a, const AutoDiff<T>& b)
    {
        AutoDiff<T> result(a.value() - b.value());
        for (size_t i = 0; i < maxSize(a, b); ++i)
            result.deriv()[i] = derivative(a, i) - derivative(b, i);
        return result;
    }

    template <typename T = double>
    auto operator- (T a, const AutoDiff<T>& b) { return AutoDiff<T>(a, b.size()) - b; }

    template <typename T = double>
    auto operator- (const AutoDiff<T>& a, T b) { return a - AutoDiff<T>(b, a.size()); }

    template <typename T = double>
    auto operator-= (AutoDiff<T>& a, const AutoDiff<T>& b)
    {
        a = a - b;
        return a;
    }

    template <typename T = double>
    AutoDiff<T> operator* (const AutoDiff<T>& a, const AutoDiff<T>& b)
    {
        size_t n = maxSize(a, b);
        AutoDiff<T> result(a.value() * b.value(), n);
        for (size_t i = 0; i < n; ++i)
            result.deriv()[i] = derivative(a, i) * b.value() + a.value() * derivative(b, i);
        return result;
    }

    template <typename T = double>
    auto operator* (T a, const AutoDiff<T>& b) { return AutoDiff<T>(a, b.size()) * b; }

    template <typename T = double>
    auto operator* (const AutoDiff<T>& a, T b) { return a * AutoDiff<T>(b, a.size()); }

    template <typename T = double>
    auto operator*= (AutoDiff<T>& a, const AutoDiff<T>& b)
    {
        a = a * b;
        return a;
    }

    template <typename T = double>
    AutoDiff<T> operator/ (const AutoDiff<T>& a, const AutoDiff<T>& b)
    {
        size_t n = maxSize(a, b);
        AutoDiff<T> result(a.value() / b.value(), n);
        for (size_t i = 0; i < n; i++)
            result.deriv()[i] = (derivative(a, i) * b.value() - a.value() * derivative(b, i)) / (b.value() * b.value());
        return result;
    }

    template <typename T = double>
    auto operator/ (T a, const AutoDiff<T>& b) { return AutoDiff<T>(a, b.size()) / b; }

    template <typename T = double>
    auto operator/ (const AutoDiff<T>& a, T b) { return a / AutoDiff<T>(b, a.size()); }

    using std::sin;
    using std::cos;
    using std::exp;
    using std::log;
    using std::sqrt;

    template <typename T = double>
    AutoDiff<T> sin(const AutoDiff<T> &a)
    {
        AutoDiff<T> result(sin(a.value()));
        for (size_t i = 0; i < a.size(); ++i)
            result.deriv()[i] = cos(a.value()) * derivative(a, i);
        return result;
    }

    template <typename T = double>
    AutoDiff<T> cos(const AutoDiff<T> &a)
    {
        AutoDiff<T> result(cos(a.value()));
        for (size_t i = 0; i < a.size(); ++i)
            result.deriv()[i] = -sin(a.value()) * derivative(a, i);
        return result;
    }

    template <typename T = double>
    AutoDiff<T> exp(const AutoDiff<T> &a)
    {
        AutoDiff<T> result(exp(a.value()));
        for (size_t i = 0; i < a.size(); ++i)
            result.deriv()[i] = exp(a.value()) * derivative(a, i);
        return result;
    }

    template <typename T = double>
    AutoDiff<T> log(const AutoDiff<T> &a)
    {
        AutoDiff<T> result(log(a.value()));
        for (size_t i = 0; i < a.size(); ++i)
            result.deriv()[i] = (1 / a.value()) * derivative(a, i);
        return result;
    }

    template <typename T = double>
    AutoDiff<T> sqrt(const AutoDiff<T>& a)
    {
        AutoDiff<T> result(sqrt(a.value()));
        for (size_t i = 0; i < a.size(); ++i)
            result.deriv()[i] = (T(0.5) / sqrt(a.value())) * derivative(a, i);
        return result;
    }

    template <typename T = double>
    AutoDiff<T> norm2(const AutoDiff<T>& a)
    {
        return a * a;
    }
    
    template <typename T>
    void LegendrePolynomials(int n, T x, std::vector<T>& P) {
        if (n < 0) {
            P.clear();
            return;
        }
        P.resize(n + 1);
        P[0] = T(1);
        if (n == 0) return;
        P[1] = x;
        for (int k = 2; k <= n; ++k) {
            P[k] = ((T(2 * k - 1) * x * P[k - 1]) - T(k - 1) * P[k - 2]) / T(k);
        }
    }

} // namespace ASC_ode

#endif
