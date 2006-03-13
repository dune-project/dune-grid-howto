// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_HH
#define UNITCUBE_HH

#include "dune/grid/common/grid.hh"
#include "dune/grid/onedgrid.hh"
#include "dune/grid/sgrid.hh"
#include "dune/grid/yaspgrid.hh"

#if HAVE_UG
#include "grid/uggrid.hh"
#include <io/file/amirameshreader.hh>
#endif

#if HAVE_ALUGRID
#include <dune/grid/alu3dgrid.hh>
#endif

#if HAVE_ALBERTA
#include "dune/grid/albertagrid.hh"
#endif

// default implementation for any template parameter
template<typename T, int variant>
class UnitCube
{
public:

  // constructor throwing exception
  UnitCube ()
  {
    DUNE_THROW(Dune::GridError,"no specialization for this grid available");
  }

  const T& grid ()
  {
    return grid_;
  }

private:
  // the constructed grid object
  T grid_;
};


// OneDGrid specialization
template<>
class UnitCube<Dune::OneDGrid<1,1>,1>
{
public:
  UnitCube () : grid_(1,0.0,1.0)
  {}

  const Dune::OneDGrid<1,1>& grid ()
  {
    return grid_;
  }

private:
  Dune::OneDGrid<1,1> grid_;
};

// SGrid specialization
template<int dim>
class UnitCube<Dune::SGrid<dim,dim>,1>
{
public:
  const Dune::SGrid<dim,dim>& grid ()
  {
    return grid_;
  }

private:
  Dune::SGrid<dim,dim> grid_;
};

// YaspGrid specialization
template<int dim>
class UnitCube<Dune::YaspGrid<dim,dim>,1>
{
public:

  UnitCube () : Len(1.0), s(1), p(false),
#if HAVE_MPI
                grid_(MPI_COMM_WORLD,Len,s,p,0)
#else
                grid_(Len,s,p,0)
#endif
  {  }

  const Dune::YaspGrid<dim,dim>& grid ()
  {
    return grid_;
  }

private:
  Dune::FieldVector<double,dim> Len;
  Dune::FieldVector<int,dim> s;
  Dune::FieldVector<bool,dim> p;
  Dune::YaspGrid<dim,dim> grid_;
};


#if HAVE_UG
// UGGrid 3d, variant 1 (hexahedra) specialization
template<>
class UnitCube<Dune::UGGrid<3,3>,1>
{
public:
  UnitCube () : grid_(800,10)
  {
    Dune::AmiraMeshReader<Dune::UGGrid<3,3> >::read(grid_,"grids/ug3dhexagrid.am");
  }

  const Dune::UGGrid<3,3>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<3,3> grid_;
};

// UGGrid 3d, variant 2 (tetrahedra) specialization
template<>
class UnitCube<Dune::UGGrid<3,3>,2>
{
public:
  UnitCube () : grid_(800,10)
  {
    Dune::AmiraMeshReader<Dune::UGGrid<3,3> >::read(grid_,"grids/ug3dtetragrid.am");
  }

  const Dune::UGGrid<3,3>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<3,3> grid_;
};

// UGGrid 2d, variant 1 (quadrilaterals) specialization
template<>
class UnitCube<Dune::UGGrid<2,2>,1>
{
public:
  UnitCube () : grid_(800,10)
  {
    Dune::AmiraMeshReader<Dune::UGGrid<2,2> >::read(grid_,"grids/quadgrid.am");
  }

  const Dune::UGGrid<2,2>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<2,2> grid_;
};

// UGGrid 2d, variant 2 (triangles) specialization
template<>
class UnitCube<Dune::UGGrid<2,2>,2>
{
public:
  UnitCube () : grid_(800,10)
  {
    Dune::AmiraMeshReader<Dune::UGGrid<2,2> >::read(grid_,"grids/trianggrid.am");
  }

  const Dune::UGGrid<2,2>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<2,2> grid_;
};
#endif


#if HAVE_ALBERTA
// AlbertaGrid 2d, variant 1 (2 triangles) specialization
template<>
class UnitCube<Dune::AlbertaGrid<2,2>,1>
{
public:
  UnitCube () : grid_("grids/2dgrid.al")
  {}

  const Dune::AlbertaGrid<2,2>& grid ()
  {
    return grid_;
  }

private:
  Dune::AlbertaGrid<2,2> grid_;
};

// AlbertaGrid 3d, variant 1 (6 tetrahedra) specialization
template<>
class UnitCube<Dune::AlbertaGrid<3,3>,1>
{
public:
  UnitCube () : grid_("grids/3dgrid.al")
  {}

  const Dune::AlbertaGrid<3,3>& grid ()
  {
    return grid_;
  }

private:
  Dune::AlbertaGrid<3,3> grid_;
};
#endif

#if HAVE_ALUGRID
// ALU3dGrid tetrahedra specialization. Note: element type determined by type
template<>
class UnitCube<Dune::ALU3dGrid<3,3,Dune::tetra>,1>
{
public:
  UnitCube () : filename("grids/cube.tetra"), grid_(filename.c_str())
  {}

  const Dune::ALU3dGrid<3,3,Dune::tetra>& grid ()
  {
    return grid_;
  }

private:
  std::string filename;
  Dune::ALU3dGrid<3,3,Dune::tetra> grid_;
};

// ALU3dGrid hexahedra specialization. Note: element type determined by type
template<>
class UnitCube<Dune::ALU3dGrid<3,3,Dune::hexa>,1>
{
public:
  UnitCube () : filename("grids/cube.hexa"), grid_(filename.c_str())
  {}

  const Dune::ALU3dGrid<3,3,Dune::hexa>& grid ()
  {
    return grid_;
  }

private:
  std::string filename;
  Dune::ALU3dGrid<3,3,Dune::hexa> grid_;
};
#endif



#endif
