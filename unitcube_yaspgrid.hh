// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_YASPGRID_HH
#define UNITCUBE_YASPGRID_HH

#include "unitcube.hh"

#include <dune/grid/yaspgrid.hh>

// YaspGrid specialization
template<int dim, int size>
class UnitCube<Dune::YaspGrid<dim>,size>
{
public:
  typedef Dune::YaspGrid<dim> GridType;

  UnitCube ()
  {
    Dune::FieldVector<double,dim> length(1.0);
    Dune::array<int,dim> elements;
    std::fill(elements.begin(), elements.end(), size);
    std::bitset<dim> periodicity(0);

    grid_ = std::auto_ptr<Dune::YaspGrid<dim> >(new Dune::YaspGrid<dim>(
#if HAVE_MPI
                                                                        MPI_COMM_WORLD,
#endif
                                                                        length,elements,periodicity,1));
  }

  Dune::YaspGrid<dim>& grid ()
  {
    return *grid_;
  }

private:
  std::auto_ptr<Dune::YaspGrid<dim> > grid_;
};

#endif
