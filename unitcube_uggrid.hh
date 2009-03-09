// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_UGGRID_HH
#define UNITCUBE_UGGRID_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>

template< int dim, int variant >
class UnitCube< Dune::UGGrid< dim >, variant >
  : public BasicUnitCube< dim >
{
public:
  typedef Dune::UGGrid< dim > GridType;

private:
  GridType* grid_;

public:
  UnitCube ()
  {
    Dune::GridFactory< GridType > factory;
    BasicUnitCube< dim >::insertVertices( factory );
    if( variant == 1 )
      BasicUnitCube< dim >::insertCubes( factory );
    else if( variant == 2 )
      BasicUnitCube< dim >::insertSimplices( factory );
    else
      DUNE_THROW( Dune::NotImplemented, "Variant " << variant << " of UG unit cube not implemented." );
    grid_ = factory.createGrid();
  }

  ~UnitCube ()
  {
    delete grid_;
  }

  GridType &grid ()
  {
    return *grid_;
  }
};

#endif // #if HAVE_UG

#endif
