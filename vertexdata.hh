// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "dune/grid/common/referenceelements.hh"
#include "dune/grid/common/mcmgmapper.hh"
#include "dune/disc/functions/p1function.hh"
#include "dune/io/file/vtk/vtkwriter.hh"
#if HAVE_GRAPE
#include "dune/io/visual/grapedatadisplay.hh"
#endif

//! Parameter for mapper class
template<int dim>
struct P1Layout
{
  bool contains (int codim, Dune::GeometryType gt)
  {
    if (codim==dim) return true;
    return false;
  }
};

// demonstrate attaching data to elements
template<class G, class F>
void vertexdata (const G& grid, const F& f)
{
  // the usual stuff
  const int dim = G::dimension;
  const int dimworld = G::dimensionworld;
  typedef typename G::ctype ct;
  typedef typename G::template Codim<dim>::LeafIterator VertexLeafIterator;

  // get leaf index set type needed for mapper
  typedef typename G::template Codim<0>::LeafIndexSet IS;

  // make a mapper for codim 0 entities in the leaf grid
  Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,P1Layout>
  mapper(grid,grid.leafIndexSet());

  // allocate a vector for the data
  std::vector<double> c(mapper.size());

  // iterate through all entities of codim 0 at the leafs
  for (VertexLeafIterator it = grid.template leafbegin<dim>();
       it!=grid.template leafend<dim>(); ++it)
  {
    // evaluate functor and store value
    c[mapper.map(*it)] = f(it->geometry()[0]);
  }

  // initialize a P0 function with the vector, needed for VTK IO
  Dune::LeafP1Function<G,double> cc(grid,c);

  // generate a VTK file
  Dune::VTKWriter<G> vtkwriter(grid);
  vtkwriter.addVertexData(cc,"data");
  vtkwriter.write("vertexdata",Dune::VTKOptions::binaryappended);

  // online visualization with Grape
#if HAVE_GRAPE
  Dune::GrapeDataDisplay<G> grape(grid);
  grape.displayVector("concentration",c,grid.leafIndexSet(),1,1);
#endif
}
