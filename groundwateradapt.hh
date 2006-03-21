// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "dune/common/timer.hh"
#include "dune/grid/common/gridinfo.hh"
#include "dune/disc/functions/p1function.hh"
#include "dune/disc/operators/p1operator.hh"
#include "dune/disc/groundwater/p1groundwater.hh"
#include "dune/disc/groundwater/p1groundwaterestimator.hh"
#include "dune/io/file/vtk/vtkwriter.hh"
#include "dune/istl/bvector.hh"
#include "dune/istl/vbvector.hh"
#include "dune/istl/bcrsmatrix.hh"
#include "dune/istl/io.hh"
#include "dune/istl/gsetc.hh"
#include "dune/istl/ilu.hh"
#include "dune/istl/operators.hh"
#include "dune/istl/solvers.hh"
#include "dune/istl/preconditioners.hh"
#include "dune/istl/scalarproducts.hh"

template<class G>
void groundwateradapt (G& grid, int maxsteps, int modulo,
                       bool globalrefine=false, bool picture=false)
{
  // control the adaptive algorithm
  const double tolerance=0.01;
  const double fraction=0.15;

  // first we extract the dimensions of the grid
  const int dim = G::dimension;

  // type used for coordinates in the grid
  typedef typename G::ctype ct;
  typedef double nt;

  // Get the iterator type
  // Note the use of the typename and template keywords
  typedef typename G::template Codim<0>::LeafIterator ElementLeafIterator;

  // watches for time measurements
  Dune::Timer watch, watch2;

  // iterate through all entities of codim 0 on the given level
  grid.globalRefine(2*modulo);

  // The solution to be computed
  Dune::LeafP1Function<G,nt> u(grid);
  *u = 0;

  // adaptation loop
  for (int step=1; step<=maxsteps; step++)
  {
    std::cout << "=== [ STEP=" << step << " DOF=" <<  (*u).size() << std::endl;
    Dune::gridinfo(grid," ");

    // make  functions for solution and right hand side
    Dune::LeafP1Function<G,nt> f(grid);

    // discretize the model problem
    TestProblem<G,nt> mp;
    watch.reset();
    Dune::LeafP1OperatorAssembler<G,nt,1> A(grid);
    std::cout << "=== TIME for matrix setup is "
              <<  watch.elapsed() << " second(s)" << std::endl;
    watch.reset();
    Dune::GroundwaterEquationLocalStiffness<G,nt> lstiff(mp);
    A.assemble(lstiff,u,f);
    std::cout << "=== TIME for assemble is "
              <<  watch.elapsed() << " second(s)" << std::endl;

    // set up the solver
    watch.reset();
    typedef typename Dune::LeafP1Function<G,nt>::RepresentationType VectorType;
    typedef typename Dune::LeafP1OperatorAssembler<G,nt,1>::RepresentationType MatrixType;
    typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
    Operator op(*A);
    Dune::SeqSSOR<MatrixType,VectorType,VectorType> ssor(*A,1,1.0);
    double red=1E-4; if (step<=1) red=1E-12;
    Dune::CGSolver<VectorType> cg(op,ssor,red,10000,2);
    Dune::LoopSolver<VectorType> loop(op,ssor,red,10000,2);
    std::cout << "=== TIME for solver setup is "
              <<  watch.elapsed() << " second(s)" << std::endl;

    // solve the linear system
    Dune::InverseOperatorResult r;
    cg.apply(*u,*f,r);
    A.interpolateHangingNodes(u);

    // error estimation
    watch.reset();
    Dune::LeafP0Function<G,nt> eta2(grid);
    Dune::GroundwaterEstimator<G,nt> estimator(grid,mp);
    estimator.estimate(u,eta2);
    nt error = sqrt((*eta2).one_norm());
    std::cout << "=== TIME for error estimation is "
              <<  watch.elapsed() << " second(s)" << std::endl;
    std::cout << "=== estimated H1 error=" << error << std::endl;

    // exit after computing max number of steps
    if (step==maxsteps || error<tolerance)
      break;

    // mark elements for refinement. To be done with estimator later.
    watch.reset();
    if (globalrefine)
    {
      for (ElementLeafIterator it = grid.template leafbegin<0>();
           it!=grid.template leafend<0>(); ++it)
        grid.mark(1,it);
      std::cout << "=== " << (*eta2).size()
                << " elements out of " << (*eta2).size()
                << " refined (100%)"
                << std::endl;
    }
    else
    {
      // find threshold using bisection
      nt upper = (*eta2).infinity_norm();
      nt lower = 0.0;
      double threshold = upper;
      int refinedelements = 0;
      double refine_target = fraction*100.0;
      double refine_actual = 0;
      int count=0;
      for (int k=0; k<10; k++)
      {
        // determine new threshold and # refined elements
        threshold = 0.5*(upper+lower);
        count = 0;
        for (int i=0; i<(*eta2).size(); i++)
          if ((*eta2)[i]>threshold) count++;
        refinedelements = count;
        refine_actual = 100.0*(((double)count)/((double)(*eta2).size()));

        // is it good enough
        if (std::abs(refine_actual-refine_target)<0.05*refine_target)
          break;

        // prepare next round
        if (refine_actual<refine_target)
          upper = threshold;
        else
          lower = threshold;
      }
      std::cout << "=== " << count
                << " elements out of " << (*eta2).size()
                << " refined (" << refine_actual << "%)"
                << std::endl;
      A.preMark();
      for (ElementLeafIterator it = grid.template leafbegin<0>();
           it!=grid.template leafend<0>(); ++it)
      {
        if ((*eta2)[eta2.mapper().map(*it)]>threshold)
          A.mark(grid,it);
      }
      A.postMark(grid);
    }

    // adapt the grid
    Dune::P1FunctionManager<G,nt> manager(grid);
    grid.preAdapt();
    u.preAdapt();
    watch2.reset();
    grid.adapt();
    std::cout << "=== TIME for grid adapt only is "
              <<  watch2.elapsed() << " second(s)" << std::endl;
    u.postAdapt(manager);
    grid.postAdapt();
    std::cout << "=== TIME for grid adaptation (total) is "
              <<  watch.elapsed() << " second(s)" << std::endl;
    std::cout << "=== ]" << std::endl;
  }

  if (picture)
  {
    watch.reset();
    std::ostringstream os;
    os << "u." << transformToGridName(grid.type()) << "." << dim << "d";
    if (globalrefine)
      os << ".global";
    else
      os << ".local";
    Dune::VTKWriter<G> vtkwriter(grid);
    vtkwriter.addVertexData(u,"solution");
    std::string s(os.str());
    vtkwriter.write(s.c_str(),Dune::VTKOptions::binaryappended);
    std::cout << "=== TIME for VTK output of " << os.str()
              << " is " <<  watch.elapsed() << " second(s)" << std::endl;
  }
}
