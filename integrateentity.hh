// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_INTEGRATE_ENTITY_HH
#define DUNE_INTEGRATE_ENTITY_HH

#include "dune/common/exceptions.hh"
#include "dune/quadrature/quadraturerules.hh"

//! compute integral of function over entity with error estimate
template<class Iterator, class Functor>
double integrate (const Iterator& it, const Functor& f, double& localerror)
{
  // dimension of the entity
  const int dim = Iterator::Entity::dimension;

  // type used for coordinates in the grid
  typedef typename Iterator::Entity::ctype ct;

  // get geometry type
  Dune::GeometryType gt = it->geometry().type();

  // get low order quadrature
  const Dune::QuadratureRule<ct,dim>&
  lowrule = Dune::QuadratureRules<ct,dim>::rule(gt,1);

  // get higher order quadrature
  const Dune::QuadratureRule<ct,dim>&
  highrule = Dune::QuadratureRules<ct,dim>::rule(gt,3);

  // ensure that highrule has higher order
  if (!(highrule.order()>lowrule.order()))
    DUNE_THROW(Dune::Exception,"high order rule not better than low order rule");

  // compute integral with low order rule
  double lowresult=0;
  for (typename Dune::QuadratureRule<ct,dim>::const_iterator i=lowrule.begin();
       i!=lowrule.end(); ++i)
  {
    double fval = f(it->geometry().global(i->position()));
    double weight = i->weight();
    double detjac = it->geometry().integrationElement(i->position());
    lowresult += fval * weight * detjac;
  }

  // compute integral with high order rule
  double highresult=0;
  for (typename Dune::QuadratureRule<ct,dim>::const_iterator i=highrule.begin();
       i!=highrule.end(); ++i)
  {
    double fval = f(it->geometry().global(i->position()));
    double weight = i->weight();
    double detjac = it->geometry().integrationElement(i->position());
    highresult += fval * weight * detjac;
  }

  // compute estimate for local error
  localerror = std::abs(lowresult-highresult);

  // return the more accurate integral
  return highresult;
}


//! compute integral of function over entity with given order
template<class Iterator, class Functor>
double integrateentity (const Iterator& it, const Functor& f, int p)
{
  // dimension of the entity
  const int dim = Iterator::Entity::dimension;

  // type used for coordinates in the grid
  typedef typename Iterator::Entity::ctype ct;

  // get geometry type
  Dune::GeometryType gt = it->geometry().type();

  // get quadrature rule of order p
  const Dune::QuadratureRule<ct,dim>&
  rule = Dune::QuadratureRules<ct,dim>::rule(gt,p);

  // ensure that rule has at least the requested order
  if (rule.order()<p)
    DUNE_THROW(Dune::Exception,"order not available");

  // compute integral with low order rule
  double result=0;
  for (typename Dune::QuadratureRule<ct,dim>::const_iterator i=rule.begin();
       i!=rule.end(); ++i)
  {
    double fval = f(it->geometry().global(i->position()));
    double weight = i->weight();
    double detjac = it->geometry().integrationElement(i->position());
    result += fval * weight * detjac;
  }

  // return the more accurate integral
  return result;
}

#endif
