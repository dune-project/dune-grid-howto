// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __DUNE_GRID_HOWTO_EVOLVE_HH__
#define __DUNE_GRID_HOWTO_EVOLVE_HH__

#include <dune/common/fvector.hh>

template<class G, class M, class V>
void evolve (const G& grid, const M& mapper, V& c, double t, double& dt)
{
  // first we extract the dimensions of the grid
  const int dimworld = G::dimensionworld;

  // type used for coordinates in the grid
  typedef typename G::ctype ct;

  // type of grid view on leaf part
  typedef typename G::LeafGridView GridView;

  // element iterator type
  typedef typename GridView::template Codim<0>::Iterator LeafIterator;

  // leaf entity geometry
  typedef typename LeafIterator::Entity::Geometry LeafGeometry;

  // intersection iterator type
  typedef typename GridView::IntersectionIterator IntersectionIterator;

  // intersection geometry
  typedef typename IntersectionIterator::Intersection::Geometry IntersectionGeometry;

  // entity type
  typedef typename G::template Codim<0>::Entity Entity;

  // get grid view on leaf part
  GridView gridView = grid.leafGridView();

  // allocate a temporary vector for the update
  V update(c.size());                                  /*@\label{evh:update}@*/
  for (typename V::size_type i=0; i<c.size(); i++) update[i] = 0;

  // initialize dt very large
  dt = 1E100;

  // compute update vector and optimum dt in one grid traversal
  LeafIterator endit = gridView.template end<0>();     /*@\label{evh:loop0}@*/
  for (LeafIterator it = gridView.template begin<0>(); it!=endit; ++it)
  {
    // cell geometry
    const LeafGeometry geo = it->geometry();


    // cell volume, assume linear map here
    double volume = geo.volume();

    // cell index
    int indexi = mapper.index(*it);

    // variable to compute sum of positive factors
    double sumfactor = 0.0;

    // run through all intersections with neighbors and boundary
    IntersectionIterator isend = gridView.iend(*it);       /*@\label{evh:flux0}@*/
    for (IntersectionIterator is = gridView.ibegin(*it); is!=isend; ++is)
    {
      // get geometry type of face
      const IntersectionGeometry igeo = is->geometry();

      // get normal vector scaled with volume
      Dune::FieldVector<ct,dimworld> integrationOuterNormal
        = is->centerUnitOuterNormal();
      integrationOuterNormal *= igeo.volume();

      // center of face in global coordinates
      Dune::FieldVector<ct,dimworld> faceglobal = igeo.center();

      // evaluate velocity at face center
      Dune::FieldVector<double,dimworld> velocity = u(faceglobal,t);

      // compute factor occuring in flux formula
      double factor = velocity*integrationOuterNormal/volume;

      // for time step calculation
      if (factor>=0) sumfactor += factor;

      // handle interior face
      if (is->neighbor())             // "correct" version /*@\label{evh:neighbor}@*/
      {
        // access neighbor
        Entity outside = is->outside();
        int indexj = mapper.index(outside);

        // compute flux from one side only
        if (indexi<indexj)
        {
          // compute factor in neighbor
          const LeafGeometry nbgeo = outside.geometry();
          double nbvolume = nbgeo.volume();
          double nbfactor = velocity*integrationOuterNormal/nbvolume;

          if (factor<0)                         // inflow
          {
            update[indexi] -= c[indexj]*factor;
            update[indexj] += c[indexj]*nbfactor;
          }
          else                         // outflow
          {
            update[indexi] -= c[indexi]*factor;
            update[indexj] += c[indexi]*nbfactor;
          }
        }
      }

      // handle boundary face
      if (is->boundary())                               /*@\label{evh:bndry}@*/
      {
        if (factor<0)                 // inflow, apply boundary condition
          update[indexi] -= b(faceglobal,t)*factor;
        else                 // outflow
          update[indexi] -= c[indexi]*factor;
      }
    }             // end all intersections             /*@\label{evh:flux1}@*/

    // compute dt restriction
    dt = std::min(dt,1.0/sumfactor);                   /*@\label{evh:dt}@*/

  }       // end grid traversal                        /*@\label{evh:loop1}@*/

  // scale dt with safety factor
  dt *= 0.99;                                          /*@\label{evh:.99}@*/

  // update the concentration vector
  for (unsigned int i=0; i<c.size(); ++i)
    c[i] += dt*update[i];                              /*@\label{evh:updc}@*/

  return;
}

#endif //__DUNE_GRID_HOWTO_EVOLVE_HH__
