// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_SGRID_HH
#define UNITCUBE_SGRID_HH

#include <dune/grid/sgrid.hh>

// SGrid specialization
template<int dim>
class UnitCube<Dune::SGrid<dim,dim>,1>
{
public:
  typedef Dune::SGrid<dim,dim> GridType;

  Dune::SGrid<dim,dim>& grid ()
  {
    return grid_;
  }

private:
  Dune::SGrid<dim,dim> grid_;
};

#endif
