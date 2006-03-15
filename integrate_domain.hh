// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_INTEGRATE_DOMAIN_HH
#define DUNE_INTEGRATE_DOMAIN_HH

#include "dune/quadrature/quadraturerules.hh"

template<class G, class Functor>
double integrate (G& grid, Functor f, int p)
{
  // first we extract the dimensions of the grid
  const int dim = G::dimension;

  // type used for coordinates in the grid
  typedef typename G::ctype ct;

  // iterate through all entities of codim 0 on the given level
  typedef typename G::template Codim<0>::LeafIterator LeafIterator;
  double sum = 0.0;
  LeafIterator eendit = grid.template leafend<0>();
  for (LeafIterator it = grid.template leafbegin<0>(); it!=eendit; ++it)
  {
    Dune::GeometryType gt = it->geometry().type();
    const Dune::QuadratureRule<ct,dim>&
    myrule = Dune::QuadratureRules<ct,dim>::rule(gt,p);
    for (typename Dune::QuadratureRule<ct,dim>::const_iterator i=myrule.begin();
         i!=myrule.end(); ++i)
    {
      double fval = f(it->geometry().global(i->position()));
      double weight = i->weight();
      double detjac = it->geometry().integrationElement(i->position());
      sum += fval * weight * detjac;
    }
  }

  return sum;
}
#endif
