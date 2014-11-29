// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __DUNE_GRID_HOWTO_VERTEXDATA_HH__
#define __DUNE_GRID_HOWTO_VERTEXDATA_HH__

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

//! Parameter for mapper class
/** This class is only here to show what such a class looks like -- it does
    exactly the same as Dune::MCMGVertexLayout. */
template<int dimgrid>
struct P1Layout
{
  bool contains (Dune::GeometryType gt)
  {
    if (gt.dim()==0) return true;
    return false;
  }
};

// demonstrate attaching data to elements
template<class G, class F>
void vertexdata (const G& grid, const F& f)
{
  // get dimension from Grid
  const int dim = G::dimension;
  typedef typename G::LeafGridView GridView;
  // determine type of LeafIterator for codimension = dimension
  typedef typename GridView::template Codim<dim>::Iterator VertexLeafIterator;

  // get grid view on the leaf part
  GridView gridView = grid.leafGridView();

  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P1Layout>
  mapper(grid);

  // allocate a vector for the data
  std::vector<double> c(mapper.size());

  // iterate through all entities of codim 0 at the leaves
  for (VertexLeafIterator it = gridView.template begin<dim>();
       it!=gridView.template end<dim>(); ++it)
  {
    // evaluate functor and store value
    c[mapper.map(*it)] = f(it->geometry().corner(0));
  }

  // generate a VTK file
  Dune::VTKWriter<typename G::LeafGridView> vtkwriter(grid.leafGridView());
  vtkwriter.addVertexData(c,"data");
  vtkwriter.write( "vertexdata", Dune::VTK::appendedraw );

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
                        gridView.indexSet(), // used index set
                        polynomialOrder, // polynomial order of data
                        dimRange); // dimRange of data
  }
#endif
}
#endif // __DUNE_GRID_HOWTO_VERTEXDATA_HH__
