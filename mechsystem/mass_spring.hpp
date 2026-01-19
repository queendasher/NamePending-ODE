#ifndef MASS_SPRING_HPP
#define MASS_SPRING_HPP

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;

#include <vector.hpp>
using namespace nanoblas;


template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
  Vec<D> vel = 0.0;
  Vec<D> acc = 0.0;
};


template <int D>
class Fix
{
public:
  Vec<D> pos;
};


class Connector
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};

std::ostream & operator<< (std::ostream & ost, const Connector & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  std::array<Connector,2> connectors;
};

class DistanceConstraint 
{
  public: 
    double distance;
    std::array<Connector,2> connectors;
};

template <int D>
class MassSpringSystem
{
  std::vector<Fix<D>> m_fixes;
  std::vector<Mass<D>> m_masses;
  std::vector<Spring> m_springs;
  std::vector<DistanceConstraint> m_distanceConstraints;
  Vec<D> m_gravity=0.0;
public:
  void setGravity (Vec<D> gravity) { m_gravity = gravity; }
  Vec<D> getGravity() const { return m_gravity; }

  Connector addFix (Fix<D> p)
  {
    m_fixes.push_back(p);
    return { Connector::FIX, m_fixes.size()-1 };
  }

  Connector addMass (Mass<D> m)
  {
    m_masses.push_back (m);
    return { Connector::MASS, m_masses.size()-1 };
  }
  
  size_t addSpring (Spring s) 
  {
    m_springs.push_back (s); 
    return m_springs.size()-1;
  }

  size_t addDistanceConstraint(Connector c1, Connector c2) 
  {
    Vec<D> p1 = (c1.type == Connector::FIX) ? m_fixes[c1.nr].pos : m_masses[c1.nr].pos;
    Vec<D> p2 = (c2.type == Connector::FIX) ? m_fixes[c2.nr].pos : m_masses[c2.nr].pos;
    
    double current_dist = norm(p1 - p2);
    
    DistanceConstraint dc;
    dc.distance = current_dist;
    dc.connectors = {c1, c2};
    
    m_distanceConstraints.push_back(dc);
    return m_distanceConstraints.size() - 1;
  }

  auto & fixes() { return m_fixes; } 
  auto & masses() { return m_masses; } 
  auto & springs() { return m_springs; }
  auto & distanceConstraints() { return m_distanceConstraints; }

  void getState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        valmat.row(i) = m_masses[i].pos;
        dvalmat.row(i) = m_masses[i].vel;
        ddvalmat.row(i) = m_masses[i].acc;
      }
  }

  void setState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        m_masses[i].pos = valmat.row(i);
        m_masses[i].vel = dvalmat.row(i);
        m_masses[i].acc = ddvalmat.row(i);
      }
  }
};

template <int D>
std::ostream & operator<< (std::ostream & ost, MassSpringSystem<D> & mss)
{
  ost << "fixes:" << std::endl;
  for (auto f : mss.fixes())
    ost << f.pos << std::endl;

  ost << "masses: " << std::endl;
  for (auto m : mss.masses())
    ost << "m = " << m.mass << ", pos = " << m.pos << std::endl;

  ost << "springs: " << std::endl;
  for (auto sp : mss.springs())
    ost << "length = " << sp.length << ", stiffness = " << sp.stiffness
        << ", C1 = " << sp.connectors[0] << ", C2 = " << sp.connectors[1] << std::endl;

  ost << "distance constraints: " << std::endl;
  for (auto dc : mss.distanceConstraints())
    ost << "distance = " << dc.distance
        << ", C1 = " << dc.connectors[0] << ", C2 = " << dc.connectors[1] << std::endl;

  return ost;
}


template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  virtual size_t dimX() const override { return D*mss.masses().size() + mss.distanceConstraints().size(); }
  virtual size_t dimF() const override{ return D*mss.masses().size() + mss.distanceConstraints().size(); }

  virtual void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f = 0.0;

    const auto n_masses = mss.masses().size();
    const auto n_distanceConstraints = mss.distanceConstraints().size();
    auto xmat = x.asMatrix(n_masses, D);
    auto fmat = f.asMatrix(n_masses, D);

    // Gravity
    for (size_t i = 0; i < n_masses; ++i)
      fmat.row(i) = mss.masses()[i].mass*mss.getGravity();

    // Springs
    for (auto spring : mss.springs())
    {
      auto [c1,c2] = spring.connectors;
      Vec<D> p1, p2;

      p1 = (c1.type == Connector::FIX) ? mss.fixes()[c1.nr].pos : xmat.row(c1.nr);
      p2 = (c2.type == Connector::FIX) ? mss.fixes()[c2.nr].pos : xmat.row(c2.nr);

      double force = spring.stiffness * (norm(p1-p2)-spring.length);
      Vec<D> dir12 = 1.0/norm(p1-p2) * (p2-p1);

      if (c1.type == Connector::MASS)
        fmat.row(c1.nr) += force*dir12;

      if (c2.type == Connector::MASS)
        fmat.row(c2.nr) -= force*dir12;
    }

    // Distance-Constraints
    for (size_t c_id = 0; c_id < n_distanceConstraints; ++c_id)
    {
      auto& constr = mss.distanceConstraints()[c_id];

      auto [c1, c2] = constr.connectors;
      Vec<D> p1, p2;

      p1 = (c1.type == Connector::FIX) ? mss.fixes()[c1.nr].pos : xmat.row(c1.nr);
      p2 = (c2.type == Connector::FIX) ? mss.fixes()[c2.nr].pos : xmat.row(c2.nr);
          
      double lambda = x(n_masses*D + c_id);  // Lagrange-Multiplie for this Constraint

      Vec<D> dir12 = p2 - p1;

      // Force acting on both masses
      // grad(g) = 2*(p2-p1) bei g = (p2-p1)^2 - L^2 = 0
      if (c1.type == Connector::MASS)
          fmat.row(c1.nr) += 2.0 * lambda * dir12;
      if (c2.type == Connector::MASS)
          fmat.row(c2.nr) -= 2.0 * lambda * dir12;
          
      // Equations for distance constraints (g(x) = (p2-p1)^2 - L^2 = 0)
      double dist = norm(dir12);
      f(n_masses*D + c_id) = dist*dist - constr.distance*constr.distance;
    }

    // Convert forces into accelerations (f=m*a  =>  a=f/m)
    for (size_t i = 0; i < n_masses; i++)
      fmat.row(i) *= 1.0/mss.masses()[i].mass;
  }
  
  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    // TODO: exact differentiation
    double eps = 1e-8;
    Vector<> xl(dimX()), xr(dimX()), fl(dimF()), fr(dimF());
    for (size_t i = 0; i < dimX(); i++)
      {
        xl = x;
        xl(i) -= eps;
        xr = x;
        xr(i) += eps;
        evaluate (xl, fl);
        evaluate (xr, fr);
        df.col(i) = 1/(2*eps) * (fr-fl);
      }
  }
  
};

#endif
