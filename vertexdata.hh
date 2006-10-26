// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
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
  // get dimension and coordinate type from Grid
  const int dim = G::dimension;
  typedef typename G::ctype ct;
  // dertermine type of LeafIterator for codimension = dimension
  typedef typename G::template Codim<dim>::LeafIterator VertexLeafIterator;

  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P1Layout>
  mapper(grid);

  // allocate a vector for the data
  std::vector<double> c(mapper.size());

  // iterate through all entities of codim 0 at the leafs
  for (VertexLeafIterator it = grid.template leafbegin<dim>();
       it!=grid.template leafend<dim>(); ++it)
  {
    // evaluate functor and store value
    c[mapper.map(*it)] = f(it->geometry()[0]);
  }

  // generate a VTK file
  //   Dune::LeafP1Function<G,double> cc(grid,c);
  Dune::VTKWriter<G> vtkwriter(grid);
  vtkwriter.addVertexData(c,"data");
  vtkwriter.write("vertexdata",Dune::VTKOptions::binaryappended);

  // online visualization with Grape
#if HAVE_GRAPE
  {
    const int polynomialOrder = 1; // we piecewise linear data
    const int dimRange = 1; // we have scalar data here
    // create instance of data display
    Dune::GrapeDataDisplay<G> grape(grid);
    // display data
    grape.displayVector("concentration", // name of data that appears in grape
                        c,  // data vector
                        grid.leafIndexSet(), // used index set
                        polynomialOrder, // polynomial order of data
                        dimRange); // dimRange of data
  }
#endif
}
