// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "dune/grid/common/referenceelements.hh"
#include "dune/grid/common/mcmgmapper.hh"
#include "dune/disc/functions/p0function.hh"
#include "dune/io/file/vtk/vtkwriter.hh"
#if HAVE_GRAPE
#include "dune/io/visual/grapedatadisplay.hh"
#endif

//! Parameter for mapper class
template<int dim>
struct P0Layout
{
  bool contains (int codim, Dune::GeometryType gt)
  {
    if (codim==0) return true;
    return false;
  }
};

// demonstrate attaching data to elements
template<class G, class F>
void elementdata (const G& grid, const F& f)
{
  // the usual stuff
  const int dim = G::dimension;
  const int dimworld = G::dimensionworld;
  typedef typename G::ctype ct;
  typedef typename G::template Codim<0>::LeafIterator ElementLeafIterator;

  // get leaf index set type needed for mapper
  typedef typename G::template Codim<0>::LeafIndexSet IS;

  // make a mapper for codim 0 entities in the leaf grid
  Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,P0Layout>
  mapper(grid,grid.leafIndexSet());

  // allocate a vector for the data
  std::vector<double> c(mapper.size());

  // iterate through all entities of codim 0 at the leafs
  for (ElementLeafIterator it = grid.template leafbegin<0>();
       it!=grid.template leafend<0>(); ++it)
  {
    // cell geometry type
    Dune::GeometryType gt = it->geometry().type();

    // cell center in reference element
    const Dune::FieldVector<ct,dim>&
    local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);

    // get global coordinate of cell center
    Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);

    // evaluate functor and store value
    c[mapper.map(*it)] = f(global);
  }

  // initialize a P0 function with the vector, needed for VTK IO
  Dune::LeafP0Function<G,double> cc(grid,c);

  // generate a VTK file
  Dune::VTKWriter<G> vtkwriter(grid);
  vtkwriter.addCellData(cc,"data");
  vtkwriter.write("elementdata",Dune::VTKOptions::binaryappended);

  // online visualization with Grape
#if HAVE_GRAPE
  Dune::GrapeDataDisplay<G> grape(grid);
  grape.displayVector("concentration",c,grid.leafIndexSet(),0,1);
#endif
}
