// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"               // know what grids are present
#include <iostream>               // for input/output to shell
#include <fstream>                // for input/output to files
#include <vector>                 // STL vector class
#include <dune/grid/common/mcmgmapper.hh> // mapper class

#include "unitcube.hh"
#include "transportproblem2.hh"
#include "initialize.hh"
#include "evolve.hh"
#include "finitevolumeadapt.hh"
#include "vtkout.hh"

//===============================================================
// the time loop function working for all types of grids
//===============================================================

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

template<class G>
void gnuplot (G& grid, std::vector<double>& c)
{
  // first we extract the dimensions of the grid
  const int dim = G::dimension;
  const int dimworld = G::dimensionworld;

  // type used for coordinates in the grid
  // such a type is exported by every grid implementation
  typedef typename G::ctype ct;

  // the grid has an iterator providing the access to
  // all elements (better codim 0 entities) which are leafs
  // of the refinement tree.
  // Note the use of the typename keyword and the traits class
  typedef typename G::template Codim<0>::LeafIterator ElementLeafIterator;

  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P0Layout>
  mapper(grid);

  // iterate through all entities of codim 0 at the leafs
  int count = 0;
  for (ElementLeafIterator it = grid.template leafbegin<0>();
       it!=grid.template leafend<0>(); ++it)
  {
    Dune::GeometryType gt = it->geometry().type();
    const Dune::FieldVector<ct,dim>&
    local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);
    Dune::FieldVector<ct,dimworld>
    global = it->geometry().global(local);
    std::cout << global[0] << " " << c[mapper.map(*it)] << std::endl;
    count++;
  }
}


template<class G>
void timeloop (G& grid, double tend, int lmin, int lmax)
{
  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P0Layout>
  mapper(grid);

  // allocate a vector for the concentration
  std::vector<double> c(mapper.size());

  // initialize concentration with initial values
  initialize(grid,mapper,c);
  for (int i=grid.maxLevel(); i<lmax; i++)
  {
    if (grid.maxLevel()>=lmax) break;
    finitevolumeadapt(grid,mapper,c,lmin,lmax,0);
    initialize(grid,mapper,c);
  }
  vtkout(grid,c,"concentration",0);

  double dt, t=0;
  int k=0;
  int modulo=5;
  std::cout << "s=" << grid.size(0) << " k=" << k
            << " t=" << t << std::endl;
  while (t<tend)
  {
    k++;
    evolve(grid,mapper,c,t,dt);
    t += dt;
    if (k%modulo==0) vtkout(grid,c,"concentration",k/modulo);
    std::cout << "s=" << grid.size(0) << " k=" << k
              << " t=" << t << " dt=" << dt << std::endl;
    //            finitevolumeadapt(grid,mapper,c,lmin,lmax,k);
  }
  vtkout(grid,c,"concentration",k/modulo);
}

//===============================================================
// The main function creates objects and does the time loop
//===============================================================

int main (int argc , char ** argv)
{
#if HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  UnitCube<Dune::OneDGrid<1,1>,1> uc0;
  UnitCube<Dune::SGrid<1,1>,1> uc1;
  UnitCube<Dune::YaspGrid<2,2>,1> uc;
#if HAVE_UG
  UnitCube<Dune::UGGrid<3,3>,1> uc2;
#endif
#if HAVE_ALBERTA
#if DUNE_PROBLEM_DIM==2
  UnitCube<Dune::AlbertaGrid<2,2>,1> uc3;
#endif
#endif
  //    uc3.grid().globalRefine(8);
  //    timeloop(uc3.grid(),0.5,8,18);
  //  uc2.grid().globalRefine(3);
  //  timeloop(uc2.grid(),0.5,3,7);
  uc0.grid().globalRefine(4);
  timeloop(uc0.grid(),0.5,4,9);

#if HAVE_MPI
  MPI_Finalize();
#endif

  // done
  return 0;
}
