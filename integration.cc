// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

// Dune includes
#include "config.h"
#include "dune/grid/sgrid.hh"
#include "dune/quadrature/quadraturerules.hh"

#include "integrate.hh"

template<typename ct, int dim>
class Exp {
public:
  Exp () {midpoint = 0.5;}
  double operator() (const Dune::FieldVector<ct,dim>& x) {
    Dune::FieldVector<ct,dim> y(x);
    y -= midpoint;
    return exp(-3.234*sqrt(y*y));
  }
private:
  Dune::FieldVector<ct,dim> midpoint;
};

int main()
{
  // make a grid
  const int dim=3;
  typedef Dune::SGrid<dim,dim> GridType;
  Dune::FieldVector<int,dim> N(3);
  Dune::FieldVector<GridType::ctype,dim> L(0.0);
  Dune::FieldVector<GridType::ctype,dim> H(1.0);
  GridType grid(N,L,H);

  // print some information about the grid
  std::cout << "result is "
            << integrate(grid,Exp<GridType::ctype,dim>(),2)
            << std::endl;

  // done
  return 0;
}
