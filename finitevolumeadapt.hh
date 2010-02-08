// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __DUNE_GRID_HOWTO_FINITEVOLUMEADAPT_HH__
#define __DUNE_GRID_HOWTO_FINITEVOLUMEADAPT_HH__
#include <map>
#include <cmath>

struct RestrictedValue
{
  double value;
  int count;
  RestrictedValue ()
  {
    value = 0;
    count = 0;
  }
};

template<class G, class M, class V>
bool finitevolumeadapt (G& grid, M& mapper, V& c, int lmin, int lmax, int k)
{
  // tol value for refinement strategy
  const double refinetol  = 0.05;
  const double coarsentol = 0.001;

  // type used for coordinates in the grid
  typedef typename G::ctype ct;

  // grid view types
  typedef typename G::LeafGridView LeafGridView;
  typedef typename G::LevelGridView LevelGridView;

  // iterator types
  typedef typename LeafGridView::template Codim<0>::Iterator LeafIterator;
  typedef typename LevelGridView::template Codim<0>::Iterator LevelIterator;

  // entity and entity pointer
  typedef typename G::template Codim<0>::Entity Entity;
  typedef typename G::template Codim<0>::EntityPointer EntityPointer;

  // intersection iterator type
  typedef typename LeafGridView::IntersectionIterator LeafIntersectionIterator;

  // global id set types, local means that the numbering is unique in a single process only.
  typedef typename G::template Codim<0>::LocalIdSet IdSet;
  // type for the index set, note that this is _not_ an integer
  typedef typename IdSet::IdType IdType;

  // get grid view on leaf grid
  LeafGridView leafView = grid.leafView();

  // compute cell indicators
  V indicator(c.size(),-1E100);
  double globalmax = -1E100;
  double globalmin =  1E100;
  for (LeafIterator it = leafView.template begin<0>(); /*@\label{fah:loop0}@*/
       it!=leafView.template end<0>(); ++it)
  {
    // my index
    int indexi = mapper.map(*it);

    // global min/max
    globalmax = std::max(globalmax,c[indexi]);
    globalmin = std::min(globalmin,c[indexi]);

    LeafIntersectionIterator isend = leafView.iend(*it);
    for (LeafIntersectionIterator is = leafView.ibegin(*it); is!=isend; ++is)
    {
      const typename LeafIntersectionIterator::Intersection &intersection = *is;
      if( !intersection.neighbor() )
        continue;

      // access neighbor
      const EntityPointer pOutside = intersection.outside();
      const Entity &outside = *pOutside;
      int indexj = mapper.map( outside );

      // handle face from one side only
      if ( it.level() > outside.level() ||
           (it.level() == outside.level() && indexi<indexj) )
      {
        double localdelta = std::abs(c[indexj]-c[indexi]);
        indicator[indexi] = std::max(indicator[indexi],localdelta);
        indicator[indexj] = std::max(indicator[indexj],localdelta);
      }
    }
  }                                              /*@\label{fah:loop1}@*/

  // mark cells for refinement/coarsening
  double globaldelta = globalmax-globalmin;
  int marked=0;
  for (LeafIterator it = leafView.template begin<0>(); /*@\label{fah:loop2}@*/
       it!=leafView.template end<0>(); ++it)
  {
    if (indicator[mapper.map(*it)]>refinetol*globaldelta
        && (it.level()<lmax || !it->isRegular()))
    {
      const Entity &entity = *it;
      grid.mark( 1, entity );
      ++marked;
      LeafIntersectionIterator isend = leafView.iend(entity);
      for( LeafIntersectionIterator is = leafView.ibegin(entity); is != isend; ++is )
      {
        const typename LeafIntersectionIterator::Intersection intersection = *is;
        if( !intersection.neighbor() )
          continue;

        const EntityPointer pOutside = intersection.outside();
        const Entity &outside = *pOutside;
        if( (outside.level() < lmax) || !outside.isRegular() )
          grid.mark( 1, outside );
      }
    }
    if (indicator[mapper.map(*it)]<coarsentol*globaldelta && it.level()>lmin)
    {
      grid.mark( -1, *it );
      ++marked;
    }
  }                                              /*@\label{fah:loop3}@*/
  if( marked==0 )
    return false;

  // restrict to coarse elements
  std::map<IdType,RestrictedValue> restrictionmap; // restricted concentration /*@\label{fah:loop4}@*/
  const IdSet& idset = grid.localIdSet();
  for (int level=grid.maxLevel(); level>=0; level--)
  {
    // get grid view on level grid
    LevelGridView levelView = grid.levelView(level);
    for (LevelIterator it = levelView.template begin<0>();
         it!=levelView.template end<0>(); ++it)
    {
      // get your map entry
      IdType idi = idset.id(*it);
      RestrictedValue& rv = restrictionmap[idi];

      // put your value in the map
      if (it->isLeaf())
      {
        int indexi = mapper.map(*it);
        rv.value = c[indexi];
        rv.count = 1;
      }

      // average in father
      if (it.level()>0)
      {
        EntityPointer ep = it->father();
        IdType idf = idset.id(*ep);
        RestrictedValue& rvf = restrictionmap[idf];
        rvf.value += rv.value/rv.count;
        rvf.count += 1;
      }
    }                                                  /*@\label{fah:loop5}@*/
  }
  grid.preAdapt();

  // adapt mesh and mapper
  bool rv=grid.adapt();                                /*@\label{fah:adapt}@*/
  mapper.update();                                     /*@\label{fah:update}@*/
  c.resize(mapper.size());                             /*@\label{fah:resize}@*/

  // interpolate new cells, restrict coarsened cells
  for (int level=0; level<=grid.maxLevel(); level++)   /*@\label{fah:loop6}@*/
  {
    LevelGridView levelView = grid.levelView(level);
    for (LevelIterator it = levelView.template begin<0>();
         it!=levelView.template end<0>(); ++it)
    {
      // get your id
      IdType idi = idset.id(*it);

      // check map entry
      typename std::map<IdType,RestrictedValue>::iterator rit
        = restrictionmap.find(idi);
      if (rit!=restrictionmap.end())
      {
        // entry is in map, write in leaf
        if (it->isLeaf())
        {
          int indexi = mapper.map(*it);
          c[indexi] = rit->second.value/rit->second.count;
        }
      }
      else
      {
        // value is not in map, interpolate from father element
        if (it.level()>0)
        {
          EntityPointer ep = it->father();
          IdType idf = idset.id(*ep);
          RestrictedValue& rvf = restrictionmap[idf];
          if (it->isLeaf())
          {
            int indexi = mapper.map(*it);
            c[indexi] = rvf.value/rvf.count;
          }
          else
          {
            // create new entry
            RestrictedValue& rv = restrictionmap[idi];
            rv.value = rvf.value/rvf.count;
            rv.count = 1;
          }
        }
      }
    }                                                  /*@\label{fah:loop7}@*/
  }
  grid.postAdapt();

  return rv;
}

#endif //__DUNE_GRID_HOWTO_FINITEVOLUMEADAPT_HH__
