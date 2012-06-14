// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __DUNE_GRID_HOWTO_ELEMENT_DATA_HH
#define __DUNE_GRID_HOWTO_ELEMENT_DATA_HH

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

//! Parameter for mapper class
/** This class is only here to show what such a class looks like -- it does
    exactly the same as Dune::MCMGElementLayout. */
template<int dimgrid>
struct P0Layout
{
  bool contains (Dune::GeometryType gt)
  {
    if (gt.dim()==dimgrid) return true;
    return false;
  }
};

// demonstrate attaching data to elements
template<class G, class F>
void elementdata (const G& grid, const F& f)
{
  // the usual stuff
  //const int dim = G::dimension;
  const int dimworld = G::dimensionworld;
  typedef typename G::ctype ct;
  typedef typename G::LeafGridView GridView;
  typedef typename GridView::template Codim<0>::Iterator ElementLeafIterator;
  typedef typename ElementLeafIterator::Entity::Geometry LeafGeometry;

  // get grid view on leaf part
  GridView gridView = grid.leafView();

  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P0Layout>
  mapper(grid);                                        /*@\label{edh:mapper}@*/

  // allocate a vector for the data
  std::vector<double> c(mapper.size());                /*@\label{edh:c}@*/

  // iterate through all entities of codim 0 at the leaves
  for (ElementLeafIterator it = gridView.template begin<0>(); /*@\label{edh:loop0}@*/
       it!=gridView.template end<0>(); ++it)
  {
    // cell geometry
    const LeafGeometry geo = it->geometry();

    // get global coordinate of cell center
    Dune::FieldVector<ct,dimworld> global = geo.center();

    // evaluate functor and store value
    c[mapper.map(*it)] = f(global);                    /*@\label{edh:feval}@*/
  }                                                    /*@\label{edh:loop1}@*/

  // generate a VTK file
  // Dune::LeafP0Function<G,double> cc(grid,c);
  Dune::VTKWriter<typename G::LeafGridView> vtkwriter(gridView); /*@\label{edh:vtk0}@*/
  vtkwriter.addCellData(c,"data");
  vtkwriter.write( "elementdata", Dune::VTK::appendedraw ); /*@\label{edh:vtk1}@*/

  // online visualization with Grape
#if HAVE_GRAPE                                         /*@\label{edh:grape0}@*/
  {
    const int polynomialOrder = 0; // we piecewise constant data
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
#endif                                                 /*@\label{edh:grape1}@*/
}

#endif //__DUNE_GRID_HOWTO_ELEMENT_DATA_HH
