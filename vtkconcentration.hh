// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "dune/io/file/vtk/vtkwriter.hh"
#include "dune/disc/functions/p0function.hh"

template<class G, class V>
void vtkconcentration (const G& grid, const V& c, int k)
{
  Dune::VTKWriter<G> vtkwriter(grid);
  Dune::LeafP0FunctionWrapper<G,std::vector<double> > cc(grid,c);
  vtkwriter.addCellData(cc,"concentration");
  vtkwriter.write("concentration",Dune::VTKOptions::binaryappended);
}
