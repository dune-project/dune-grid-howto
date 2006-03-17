// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <map>

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
bool finitevolumeadapt (G& grid, M& mapper, V& c, int lmin, int lmax)
{
  // tol value for refinement strategy
  const double refinetol = 0.02;
  const double coarsentol = 0.0002;

  // type used for coordinates in the grid
  typedef typename G::ctype ct;

  // iterator types
  typedef typename G::template Codim<0>::LeafIterator LeafIterator;
  typedef typename G::template Codim<0>::LevelIterator LevelIterator;

  // entity pointer
  typedef typename G::template Codim<0>::EntityPointer EntityPointer;

  // intersection iterator type
  typedef typename G::template Codim<0>::IntersectionIterator IntersectionIterator;

  // global id set types
  typedef typename G::template Codim<0>::LocalIdSet IdSet;
  typedef typename IdSet::IdType IdType;

  // compute cell indicators
  V indicator(c.size());
  double globaldelta = -1E100;
  for (LeafIterator it = grid.template leafbegin<0>();
       it!=grid.template leafend<0>(); ++it)
  {
    int indexi = mapper.map(*it);
    double localdelta = -1E100;
    IntersectionIterator isend = it->iend();
    for (IntersectionIterator is = it->ibegin(); is!=isend; ++is)
      if (is.neighbor())
      {
        int indexj = mapper.map(*(is.outside()));
        localdelta = std::max(localdelta,std::abs(c[indexj]-c[indexi]));
      }
    indicator[indexi] = localdelta;
    globaldelta = std::max(globaldelta,localdelta);
  }

  // mark cells for refinement/coarsening
  for (LeafIterator it = grid.template leafbegin<0>();
       it!=grid.template leafend<0>(); ++it)
  {
    if (indicator[mapper.map(*it)]>refinetol*globaldelta && it.level()<lmax)
      grid.mark(1,it);
    if (indicator[mapper.map(*it)]<coarsentol*globaldelta && it.level()>lmin)
      grid.mark(-1,it);
  }

  // restrict to coarse elements
  std::map<IdType,RestrictedValue> restrictionmap; // restricted concentration
  const IdSet& idset = grid.localIdSet();
  for (int level=grid.maxLevel(); level>=0; level--)
    for (LevelIterator it = grid.template lbegin<0>(level);
         it!=grid.template lend<0>(level); ++it)
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
    }
  grid.preAdapt();

  // adapt mesh and mapper
  bool rv=grid.adapt();
  mapper.update();
  c.resize(mapper.size());

  // interpolate new cells, restrict coarsened cells
  for (int level=0; level<=grid.maxLevel(); level++)
    for (LevelIterator it = grid.template lbegin<0>(level);
         it!=grid.template lend<0>(level); ++it)
    {
      // get your id
      IdType idi = idset.id(*it);

      // check map entry
      typename std::map<IdType,RestrictedValue>::iterator rit = restrictionmap.find(idi);
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
        // value is not in map, interpolate
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
    }
  grid.postAdapt();

  return rv;
}
