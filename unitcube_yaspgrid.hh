// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_YASPGRID_HH
#define UNITCUBE_YASPGRID_HH

#include <array>
#include <memory>

#include <dune/grid/yaspgrid.hh>

#include "unitcube.hh"

// YaspGrid specialization
template<int dim, int size>
class UnitCube<Dune::YaspGrid<dim>,size>
{
public:
  typedef Dune::YaspGrid<dim> GridType;

  UnitCube ()
  {
    Dune::FieldVector<double,dim> length(1.0);
    std::array<int,dim> elements;
    std::fill(elements.begin(), elements.end(), size);

    grid_ = std::unique_ptr<Dune::YaspGrid<dim> >(new Dune::YaspGrid<dim>(length,elements));
  }

  Dune::YaspGrid<dim>& grid ()
  {
    return *grid_;
  }

private:
  std::unique_ptr<Dune::YaspGrid<dim> > grid_;
};

#endif
