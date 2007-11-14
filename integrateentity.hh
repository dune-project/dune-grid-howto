// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_INTEGRATE_ENTITY_HH
#define DUNE_INTEGRATE_ENTITY_HH

#include <dune/common/exceptions.hh>
#include <dune/grid/common/quadraturerules.hh>

//! compute integral of function over entity with given order
template<class Iterator, class Functor>
double integrateentity (const Iterator& it, const Functor& f, int p)
{
  // dimension of the entity
  const int dim = Iterator::Entity::dimension;

  // type used for coordinates in the grid
  typedef typename Iterator::Entity::ctype ct;

  // get geometry type
  Dune::GeometryType gt = it->type();

  // get quadrature rule of order p
  const Dune::QuadratureRule<ct,dim>&
  rule = Dune::QuadratureRules<ct,dim>::rule(gt,p);     /*@\label{ieh:qr}@*/

  // ensure that rule has at least the requested order
  if (rule.order()<p)
    DUNE_THROW(Dune::Exception,"order not available");

  // compute approximate integral
  double result=0;
  for (typename Dune::QuadratureRule<ct,dim>::const_iterator i=rule.begin(); /*@\label{ieh:for}@*/
       i!=rule.end(); ++i)
  {
    double fval = f(it->geometry().global(i->position()));       /*@\label{ieh:fval}@*/
    double weight = i->weight();                        /*@\label{ieh:weight}@*/
    double detjac = it->geometry().integrationElement(i->position());       /*@\label{ieh:detjac}@*/
    result += fval * weight * detjac;                   /*@\label{ieh:result}@*/
  }

  // return result
  return result;
}
#endif
