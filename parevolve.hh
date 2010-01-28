// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

template<class G, class M, class V>
void parevolve (const G& grid, const M& mapper, V& c, double t, double& dt)
{
  // check data partitioning
  assert(grid.overlapSize(0)>0 || (grid.ghostSize(0)>0)); /*@\label{peh:assert}@*/

  // first we extract the dimensions of the grid
  const int dim = G::dimension;
  const int dimworld = G::dimensionworld;

  // type used for coordinates in the grid
  typedef typename G::ctype ct;

  // type for grid view on leaf part
  typedef typename G::LeafGridView GridView;

  // iterator type
  typedef typename GridView::template Codim<0>::
  template Partition<Dune::All_Partition>::Iterator LeafIterator;       /*@\label{peh:pit}@*/

  // leaf entity geometry
  typedef typename LeafIterator::Entity::Geometry LeafGeometry;

  // intersection iterator type
  typedef typename GridView::IntersectionIterator IntersectionIterator;

  // type of intersection
  typedef typename IntersectionIterator::Intersection Intersection;

  // intersection geometry
  typedef typename Intersection::Geometry IntersectionGeometry;

  // entity pointer type
  typedef typename G::template Codim<0>::EntityPointer EntityPointer;

  // allocate a temporary vector for the update
  V update(c.size());
  for (typename V::size_type i=0; i<c.size(); i++) update[i] = 0;

  // initialize dt very large
  dt = 1E100;

  // get grid view instance on leaf grid
  GridView gridView = grid.leafView();

  // compute update vector and optimum dt in one grid traversal
  // iterate over all entities, but update is only used on interior entities
  LeafIterator endit = gridView.template end<0,Dune::All_Partition>(); /*@\label{peh:end}@*/
  for (LeafIterator it = gridView.template begin<0,Dune::All_Partition>(); it!=endit; ++it) /*@\label{peh:begin}@*/
  {
    // cell geometry
    const LeafGeometry & gt = it->geometry();

    // cell volume
    double volume = gt.volume();

    // cell index
    int indexi = mapper.map(*it);

    // variable to compute sum of positive factors
    double sumfactor = 0.0;

    // run through all intersections with neighbors and boundary
    const IntersectionIterator isend = gridView.iend(*it);
    for( IntersectionIterator is = gridView.ibegin(*it); is != isend; ++is )
    {
      const Intersection &intersection = *is;

      // get geometry type of face
      const IntersectionGeometry & gtf = intersection.geometry();

      // get normal vector scaled with volume
      Dune::FieldVector< ct, dimworld > integrationOuterNormal
        = intersection.centerUnitOuterNormal();
      integrationOuterNormal *= gtf.volume();

      // center of face in global coordinates
      Dune::FieldVector< ct, dimworld > faceglobal = gtf.center();

      // evaluate velocity at face center
      Dune::FieldVector<double,dim> velocity = u(faceglobal,t);

      // compute factor occuring in flux formula
      double factor = velocity*integrationOuterNormal/volume;

      // for time step calculation
      if (factor>=0) sumfactor += factor;

      // handle interior face
      if( intersection.neighbor() )
      {
        // access neighbor
        EntityPointer outside = intersection.outside();
        int indexj = mapper.map(*outside);

        const int insideLevel = it->level();
        const int outsideLevel = outside->level();

        // handle face from one side
        if( (insideLevel > outsideLevel)
            || ((insideLevel == outsideLevel) && (indexi < indexj)) )
        {
          // compute factor in neighbor
          const LeafGeometry & nbgt = outside->geometry();
          double nbvolume = nbgt.volume();
          double nbfactor = velocity*integrationOuterNormal/nbvolume;

          if( factor < 0 )       // inflow
          {
            update[indexi] -= c[indexj]*factor;
            update[indexj] += c[indexj]*nbfactor;
          }
          else       // outflow
          {
            update[indexi] -= c[indexi]*factor;
            update[indexj] += c[indexi]*nbfactor;
          }
        }
      }

      // handle boundary face
      if( intersection.boundary() )
      {
        if( factor < 0 )       // inflow, apply boundary condition
          update[indexi] -= b(faceglobal,t)*factor;
        else       // outflow
          update[indexi] -= c[indexi]*factor;
      }
    }       // end all intersections

    // compute dt restriction
    if (it->partitionType()==Dune::InteriorEntity)       /*@\label{peh:inter}@*/
      dt = std::min(dt,1.0/sumfactor);

  }       // end grid traversal

  // global min over all partitions
  dt = grid.comm().min(dt);                            /*@\label{peh:min}@*/
  // scale dt with safety factor
  dt *= 0.99;

  // exchange update
  VectorExchange<M,V> dh(mapper,update);               /*@\label{peh:dist0}@*/
  grid.template
  communicate<VectorExchange<M,V> >(dh,Dune::InteriorBorder_All_Interface,
                                    Dune::ForwardCommunication);                                       /*@\label{peh:dist1}@*/

  // update the concentration vector
  for (unsigned int i=0; i<c.size(); ++i)
    c[i] += dt*update[i];

  return;
}
