// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __DUNE_GRID_HOWTO_VTKOUT_HH__
#define __DUNE_GRID_HOWTO_VTKOUT_HH__

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <stdio.h>

template<class G, class V>
void vtkout (const G& grid, const V& c, const char* name, int k, double time=0.0, int rank=0)
{
  Dune::VTKWriter<typename G::LeafGridView> vtkwriter(grid.leafGridView());
  char fname[128];
  char sername[128];
  sprintf(fname,"%s-%05d",name,k);
  sprintf(sername,"%s.series",name);
  vtkwriter.addCellData(c,"celldata");
  vtkwriter.write( fname, Dune::VTK::ascii );

  if ( rank == 0)
  {
    std::ofstream serstream(sername, (k==0 ? std::ios_base::out : std::ios_base::app));
    serstream << k << " " << fname << ".vtu " << time << std::endl;
    serstream.close();
  }
}

#endif // __DUNE_GRID_HOWTO_VTKOUT_HH__
