// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_YASPGRID_HH
#define UNITCUBE_YASPGRID_HH

#include <dune/grid/yaspgrid.hh>

// YaspGrid specialization
template<int dim, int size>
class UnitCube<Dune::YaspGrid<dim>,size>
{
public:
  typedef Dune::YaspGrid<dim> GridType;

  UnitCube () : Len(1.0), s(size), p(false),
#if HAVE_MPI
                grid_(MPI_COMM_WORLD,Len,s,p,1)
#else
                grid_(Len,s,p,1)
#endif
  {  }

  Dune::YaspGrid<dim>& grid ()
  {
    return grid_;
  }

private:
  Dune::FieldVector<double,dim> Len;
  Dune::FieldVector<int,dim> s;
  Dune::FieldVector<bool,dim> p;
  Dune::YaspGrid<dim> grid_;
};

#endif
