// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "dune/disc/operators/boundaryconditions.hh"
#include "dune/disc/groundwater/groundwater.hh"

template<class G, class RT>
class TestProblem : public Dune::GroundwaterEquationParameters<G,RT>
{
  typedef typename G::ctype DT;
  enum {n=G::dimension};
  typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
  TestProblem ()
  {
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        if (i==j)
          absperm[i][j] = 1;
        else
          absperm[i][j] = 0;
  }

  virtual const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x,
                                              const Entity& e,
                                              const Dune::FieldVector<DT,n>& xi) const
  {
    return absperm;
  }

  virtual RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e,
                  const Dune::FieldVector<DT,n>& xi) const
  {
    return 0;
  }

  virtual typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x,
                                                           const Entity& e, const Dune::FieldVector<DT,n>& xi) const
  {
    if (x[0]>1-1E-8)     // right plane
      return Dune::BoundaryConditions::dirichlet;
    return Dune::BoundaryConditions::neumann;
  }

  virtual RT g (const Dune::FieldVector<DT,n>& x, const Entity& e,
                const Dune::FieldVector<DT,n>& xi) const
  {
    return 0;
  }

  virtual RT J (const Dune::FieldVector<DT,n>& x, const Entity& e,
                const Dune::FieldVector<DT,n>& xi) const
  {
    if (x[0]<1E-8)     // left plane
      for (int i=1; i<n; i++)
        if (x[i]>0.5) return -1;
    return 0;
  }

private:
  Dune::FieldMatrix<DT,n,n> absperm;
};
