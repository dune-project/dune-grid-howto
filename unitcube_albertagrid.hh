// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_ALBERTAGRID_HH
#define UNITCUBE_ALBERTAGRID_HH

#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/gridfactory.hh>

template< int dim >
class UnitCube< Dune::AlbertaGrid< dim, dim >, 1 >
  : public BasicUnitCube< dim >
{
public:
  typedef Dune::AlbertaGrid< dim, dim > GridType;

private:
  GridType *grid_;

public:
  UnitCube ()
  {
    Dune::GridFactory< GridType > factory;
    BasicUnitCube< dim >::insertVertices( factory );
    BasicUnitCube< dim >::insertSimplices( factory );
    grid_ = factory.createGrid( "UnitCube", true );
  }

  ~UnitCube ()
  {
    Dune::GridFactory< GridType >::destroyGrid( grid_ );
  }

  GridType &grid ()
  {
    return *grid_;
  }
};

#endif // #if HAVE_ALBERTA

#endif
