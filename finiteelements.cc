// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>
#include <iostream>
#include <vector>
#include <set>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/albertagrid.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#else
#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>
#endif // HAVE_DUNE_ISTL

#include "shapefunctions.hh"

// P1Elements:
// a P1 finite element discretization for elliptic problems Dirichlet
// boundary conditions on simplicial conforming grids
template<class GV, class F>
class P1Elements
{
public:
  static const int dim = GV::dimension;

  typedef typename GV::ctype ctype;
#if HAVE_DUNE_ISTL
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<ctype,1,1> > Matrix;
  typedef Dune::BlockVector<Dune::FieldVector<ctype,1> > ScalarField;
#else
  typedef Dune::DynamicMatrix<ctype> Matrix;
  typedef Dune::DynamicVector<ctype> ScalarField;
#endif // HAVE_DUNE_ISTL

private:
  typedef typename GV::template Codim<0>::Iterator LeafIterator;
  typedef typename GV::template Codim<0>::Geometry::JacobianInverseTransposed JacobianInverseTransposed;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename GV::IndexSet LeafIndexSet;

  const GV& gv;
  const F& f;

public:
  Matrix A;
  ScalarField b;
  ScalarField u;
  std::vector< std::set<int> > adjacencyPattern;

  P1Elements(const GV& gv_, const F& f_) : gv(gv_), f(f_) {}

  // store adjacency information in a vector of sets
  void determineAdjacencyPattern();

  // assemble stiffness matrix A and right side b
  void assemble();

  // finally solve Au = b for u
  void solve();
};

template<class GV, class F> /*@\label{fem:adjpat1}@*/
void P1Elements<GV, F>::determineAdjacencyPattern()
{
  const int N = gv.size(dim);
  adjacencyPattern.resize(N);

  const LeafIndexSet& set = gv.indexSet();
  const LeafIterator itend = gv.template end<0>();

  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    auto geo = it->geometry();
    auto ref = referenceElement(geo);

    // traverse all codim-1-entities of the current element and store all
    // pairs of vertices in adjacencyPattern
    const IntersectionIterator isend = gv.iend(*it);
    for (IntersectionIterator is = gv.ibegin(*it) ; is != isend ; ++is)
    {
      int vertexsize = ref.size(is->indexInInside(),1,dim);
      for (int i=0; i < vertexsize; i++)
      {
        int indexi = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,i,dim),dim);
        for (int j=0; j < vertexsize; j++)
        {
          int indexj = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,j,dim),dim);
          adjacencyPattern[indexi].insert(indexj);
        }
      }
    }
  }
} /*@\label{fem:adjpat2}@*/

template<class GV, class F>
void P1Elements<GV, F>::assemble()
{
  const int N = gv.size(dim);

  const LeafIndexSet& set = gv.indexSet();
  const LeafIterator itend = gv.template end<0>();

  // set sizes of A and b
#if HAVE_DUNE_ISTL
  A.setSize(N, N, N + 2*gv.size(dim-1));
  A.setBuildMode(Matrix::random);
  b.resize(N);

  for (int i = 0; i < N; i++)
    A.setrowsize(i,adjacencyPattern[i].size());
  A.endrowsizes();

  // set sparsity pattern of A with the information gained in determineAdjacencyPattern
  for (int i = 0; i < N; i++)   /*@\label{fem:setpattern}@*/
  {
    std::template set<int>::iterator setend = adjacencyPattern[i].end();
    for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();
         setit != setend; ++setit)
      A.addindex(i,*setit);
  }

  A.endindices();   /*@\label{fem:endindices}@*/
#else
  A.resize(N, N);
  b.resize(N);
#endif // HAVE_DUNE_ISTL

  // initialize A and b
  A = 0.0;
  b = 0.0;

  // get a set of P1 shape functions
  const P1ShapeFunctionSet<ctype,ctype,dim>& basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)   /*@\label{fem:loop1}@*/
  {
    // determine geometry of the current element and get the matching reference element
    auto geo = it->geometry();
    auto ref = referenceElement(geo);
    int vertexsize = ref.size(dim);

    // get a quadrature rule of order one for the given geometry type
    Dune::GeometryType gt = it->type();
    const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,1);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
         r != rule.end() ; ++r)
    {
      // compute the jacobian inverse transposed to transform the gradients
      JacobianInverseTransposed jacInvTra =
        it->geometry().jacobianInverseTransposed(r->position());

      // get the weight at the current quadrature point
      ctype weight = r->weight();

      // compute Jacobian determinant for the transformation formula
      ctype detjac = it->geometry().integrationElement(r->position());
      for (int i = 0; i < vertexsize; i++)
      {
        // compute transformed gradients
        Dune::FieldVector<ctype,dim> grad1;
        jacInvTra.mv(basis[i].evaluateGradient(r->position()),grad1);
        for (int j = 0; j < vertexsize; j++)
        {
          Dune::FieldVector<ctype,dim> grad2;
          jacInvTra.mv(basis[j].evaluateGradient(r->position()),grad2);

          // gain global inidices of vertices i and j and update associated matrix entry
          A[set.subIndex(*it,i,dim)][set.subIndex(*it,j,dim)]           /*@\label{fem:calca}@*/
            += (grad1*grad2) * weight * detjac;
        }
      }
    }

    // get a quadrature rule of order two for the given geometry type
    const Dune::QuadratureRule<ctype,dim>& rule2 = Dune::QuadratureRules<ctype,dim>::rule(gt,2);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule2.begin();
         r != rule2.end() ; ++r)
    {
      ctype weight = r->weight();
      ctype detjac = it->geometry().integrationElement(r->position());
      for (int i = 0 ; i<vertexsize; i++)
      {
        // evaluate the integrand of the right side
        ctype fval = basis[i].evaluateFunction(r->position())
                     * f(it->geometry().global(r->position())) ;
        b[set.subIndex(*it,i,dim)] += fval * weight * detjac;         /*@\label{fem:calcb}@*/
      }
    }
  }   /*@\label{fem:loop2}@*/

  // Dirichlet boundary conditions:
  // replace lines in A related to Dirichlet vertices by trivial lines
  for ( LeafIterator it = gv.template begin<0>() ; it != itend ; ++it)   /*@\label{fem:boundary1}@*/
  {
    const IntersectionIterator isend = gv.iend(*it);
    for (IntersectionIterator is = gv.ibegin(*it) ; is != isend ; ++is)
    {
      // determine geometry of the current element and get the matching reference element
      auto geo = it->geometry();
      auto ref = referenceElement(geo);

      // check whether current intersection is on the boundary
      if ( is->boundary() )
      {
        // traverse all vertices the intersection consists of
        for (int i=0; i < ref.size(is->indexInInside(),1,dim); i++)
        {
          // and replace the associated line of A and b with a trivial one
          int indexi = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,i,dim),dim);

          A[indexi] = 0.0;           /*@\label{fem:trivialline}@*/
          A[indexi][indexi] = 1.0;
          b[indexi] = 0.0;
        }
      }
    }
  }   /*@\label{fem:boundary2}@*/
}

#if HAVE_DUNE_ISTL
template<class GV, class E>
void P1Elements<GV, E>::solve()
{
  // make linear operator from A
  Dune::MatrixAdapter<Matrix,ScalarField,ScalarField> op(A);

  // initialize preconditioner
  Dune::SeqILU<Matrix,ScalarField,ScalarField> ilu1(A, 1, 0.92);

  // the inverse operator
  Dune::BiCGSTABSolver<ScalarField> bcgs(op, ilu1, 1e-15, 5000, 0);
  Dune::InverseOperatorResult r;

  // initialize u to some arbitrary value to avoid u being the exact
  // solution
  u.resize(b.N());
  u = 2.0;

  // finally solve the system
  bcgs.apply(u, b, r);
}
#endif // HAVE_DUNE_ISTL

// an example right hand side function
template<class ctype, int dim>
class Bump {
public:
  ctype operator() (Dune::FieldVector<ctype,dim> x) const
  {
    ctype result = 0;
    for (int i=0 ; i < dim ; i++)
      result += 2.0 * x[i]* (1-x[i]);
    return result;
  }
};

int main(int argc, char** argv)
{
#if HAVE_ALBERTA && ALBERTA_DIM==2
  static const int dim = 2;                             /*@\label{fem:dim}@*/
  std::stringstream gridfile;
    gridfile << DUNE_GRID_HOWTO_EXAMPLE_GRIDS_PATH
      << "2dgrid.al";                                   /*@\label{fem:file}@*/

  typedef Dune::AlbertaGrid<dim,dim> GridType;
  typedef GridType::LeafGridView GV;

  typedef GridType::ctype ctype;
  typedef Bump<ctype,dim> Func;

  GridType grid(gridfile.str());
  const GV& gv = grid.leafGridView();

  Func f;
  P1Elements<GV,Func> p1(gv, f);

#if HAVE_DUNE_ISTL
  grid.globalRefine(16);
#else
  grid.globalRefine(10);
#endif // HAVE_DUNE_ISTL

  std::cout << "-----------------------------------" << "\n";
  std::cout << "number of unknowns: " << grid.size(dim) << "\n";

  std::cout << "determine adjacency pattern..." << "\n";
  p1.determineAdjacencyPattern();

  std::cout << "assembling..." << "\n";
  p1.assemble();

#if HAVE_DUNE_ISTL
  std::cout << "solving..." << "\n";
  p1.solve();

  std::cout << "visualizing..." << "\n";
  Dune::VTKWriter<GridType::LeafGridView> vtkwriter(grid.leafGridView());
  vtkwriter.addVertexData(p1.u, "u");
  vtkwriter.write("fem2d", Dune::VTK::appendedraw);
#else
  std::cout << "for solving and visualizing dune-istl is necessary." << "\n";
#endif // HAVE_DUNE_ISTL
#else
  std::cerr << "You need Alberta in 2d for this program." << std::endl;
#endif // HAVE_ALBERTA && ALBERTA_DIM==2
}
