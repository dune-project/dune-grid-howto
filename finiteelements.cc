// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"
#include <iostream>
#include <vector>
#include <set>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/grid/albertagrid.hh>
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
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<ctype,1,1> > Matrix;
  typedef Dune::BlockVector<Dune::FieldVector<ctype,1> > ScalarField;

private:
  typedef typename GV::template Codim<0>::Iterator LeafIterator;
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

template<class GV, class F>
void P1Elements<GV, F>::determineAdjacencyPattern()
{
  const int N = gv.size(dim);
  adjacencyPattern.resize(N);

  const LeafIndexSet& set = gv.indexSet();
  const LeafIterator itend = gv.template end<0>();

  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    Dune::GeometryType gt = it->type();
    const Dune::template GenericReferenceElement<ctype,dim> &ref =
      Dune::GenericReferenceElements<ctype,dim>::general(gt);

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
}

template<class GV, class F>
void P1Elements<GV, F>::assemble()
{
  const int N = gv.size(dim);

  const LeafIndexSet& set = gv.indexSet();
  const LeafIterator itend = gv.template end<0>();

  // set sizes of A and b
  A.setSize(N, N, N + 2*gv.size(dim-1));
  A.setBuildMode(Matrix::random);
  b.resize(N, false);

  for (int i = 0; i < N; i++)
    A.setrowsize(i,adjacencyPattern[i].size());
  A.endrowsizes();

  // set sparsity pattern of A with the information gained in determineAdjacencyPattern
  for (int i = 0; i < N; i++)
  {
    std::template set<int>::iterator setend = adjacencyPattern[i].end();
    for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();
         setit != setend; ++setit)
      A.addindex(i,*setit);
  }

  A.endindices();

  // initialize A and b
  A = 0.0;
  b = 0.0;

  // get a set of P1 shape functions
  P1ShapeFunctionSet<ctype,ctype,dim> basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    // determine geometry type of the current element and get the matching reference element
    Dune::GeometryType gt = it->type();
    const Dune::template GenericReferenceElement<ctype,dim> &ref =
      Dune::GenericReferenceElements<ctype,dim>::general(gt);
    int vertexsize = ref.size(dim);

    // get a quadrature rule of order one for the given geometry type
    const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,1);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
         r != rule.end() ; ++r)
    {
      // compute the jacobian inverse transposed to transform the gradients
      Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
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
          A[set.subIndex(*it,i,dim)][set.subIndex(*it,j,dim)]
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
        ctype fval = basis[i].evaluateFunction(it->geometry().global(r->position()))
                     * f(it->geometry().global(r->position())) ;
        b[set.subIndex(*it,i,dim)] += fval * weight * detjac;
      }
    }
  }

  // Dirichlet boundary conditions:
  // replace lines in A related to Dirichlet vertices by trivial lines
  for ( LeafIterator it = gv.template begin<0>() ; it != itend ; ++it)
  {
    const IntersectionIterator isend = gv.iend(*it);
    for (IntersectionIterator is = gv.ibegin(*it) ; is != isend ; ++is)
    {
      // determine geometry type of the current element and get the matching reference element
      Dune::GeometryType gt = it->type();
      const Dune::template GenericReferenceElement<ctype,dim> &ref =
        Dune::GenericReferenceElements<ctype,dim>::general(gt);

      // check whether current intersection is on the boundary
      if ( is->boundary() )
      {
        // traverse all vertices the intersection consists of
        for (int i=0; i < ref.size(is->indexInInside(),1,dim); i++)
        {
          // and replace the associated line of A and b with a trivial one
          int indexi = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,i,dim),dim);

          A[indexi] = 0.0;
          A[indexi][indexi] = 1.0;
          b[indexi] = 0.0;
        }
      }
    }
  }
}

template<class GV, class E>
void P1Elements<GV, E>::solve()
{
  // make linear operator from A
  Dune::MatrixAdapter<Matrix,ScalarField,ScalarField> op(A);

  // initialize preconditioner
  Dune::SeqILUn<Matrix,ScalarField,ScalarField> ilu1(A, 1, 0.92);

  // the inverse operator
  Dune::BiCGSTABSolver<ScalarField> bcgs(op, ilu1, 1e-15, 5000, 0);
  Dune::InverseOperatorResult r;

  // initialue u to some arbitrary value to avoid u being the exact
  // solution
  u.resize(b.N(), false);
  u = 2.0;

  // finally solve the system
  bcgs.apply(u, b, r);
}


// an example right hand side function
template<class ctype, int dim>
class Bump {
public:
  ctype operator() (Dune::FieldVector<ctype,dim> x) const
  {
    return 2.0 * (x[0]*(1-x[0]) + x[1]*(1-x[1]));
  }
};

int main(int argc, char** argv)
{
  static const int dim = 2;
  const char* gridfile = "grids/2dgrid.al";

#if HAVE_ALBERTA
#if ALBERTA_DIM==2

  typedef Dune::AlbertaGrid<dim,dim> GridType;
  typedef GridType::LeafGridView GV;

  typedef GridType::ctype ctype;
  typedef Bump<ctype,dim> Func;

  GridType grid(gridfile);
  const GV& gv=grid.leafView();

  Func f;
  P1Elements<GV,Func> p1(gv, f);

  grid.globalRefine(1);

  std::cout << "-----------------------------------" << "\n";
  std::cout << "number of unknowns: " << grid.size(dim) << "\n";

  std::cout << "determine adjacency pattern..." << "\n";
  p1.determineAdjacencyPattern();

  std::cout << "assembling..." << "\n";
  p1.assemble();

  std::cout << "solving..." << "\n";
  p1.solve();

  std::cout << "visualizing..." << "\n";
  Dune::VTKWriter<GridType::LeafGridView> vtkwriter(grid.leafView());
  vtkwriter.addVertexData(p1.u, "u");
  vtkwriter.write("test", Dune::VTKOptions::binaryappended);

#endif
#endif
}
