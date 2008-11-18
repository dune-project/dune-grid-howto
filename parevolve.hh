// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/grid/common/referenceelements.hh>

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

  // iterator type
  typedef typename G::template Codim<0>::
  template Partition<Dune::All_Partition>::LeafIterator LeafIterator;       /*@\label{peh:pit}@*/

  // intersection iterator type
  typedef typename G::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

  // type of intersection
  typedef typename IntersectionIterator::Intersection Intersection;

  // entity pointer type
  typedef typename G::template Codim<0>::EntityPointer EntityPointer;

  // allocate a temporary vector for the update
  V update(c.size());
  for (typename V::size_type i=0; i<c.size(); i++) update[i] = 0;

  // initialize dt very large
  dt = 1E100;

  // compute update vector and optimum dt in one grid traversal
  // iterate over all entities, but update is only used on interior entities
  LeafIterator endit = grid.template leafend<0,Dune::All_Partition>(); /*@\label{peh:end}@*/
  for (LeafIterator it = grid.template leafbegin<0,Dune::All_Partition>(); it!=endit; ++it) /*@\label{peh:begin}@*/
  {
    // cell geometry type
    Dune::GeometryType gt = it->type();

    // cell center in reference element
    const Dune::FieldVector<ct,dim>&
    local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);

    // cell center in global coordinates
    Dune::FieldVector<ct,dimworld>
    global = it->geometry().global(local);

    // cell volume, assume linear map here
    double volume = it->geometry().integrationElement(local)
                    *Dune::ReferenceElements<ct,dim>::general(gt).volume();

    // cell index
    int indexi = mapper.map(*it);

    // variable to compute sum of positive factors
    double sumfactor = 0.0;

    // run through all intersections with neighbors and boundary
    const IntersectionIterator isend = it->ileafend();
    for( IntersectionIterator is = it->ileafbegin(); is != isend; ++is )
    {
      const Intersection &intersection = *is;

      // get geometry type of face
      Dune::GeometryType gtf = intersection.intersectionSelfLocal().type();

      const Dune::ReferenceElement< ct, dim-1 > &refElement
        = Dune::ReferenceElements< ct, dim-1 >::general( gtf );

      // center in face's reference element
      const Dune::FieldVector< ct, dim-1 > &facelocal = refElement.position( 0, 0 );

      // get normal vector scaled with volume
      Dune::FieldVector< ct, dimworld > integrationOuterNormal
        = intersection.integrationOuterNormal( facelocal );
      integrationOuterNormal *= refElement.volume();

      // center of face in global coordinates
      Dune::FieldVector< ct, dimworld > faceglobal
        = intersection.intersectionGlobal().global( facelocal );

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
          Dune::GeometryType nbgt = outside->type();
          const Dune::FieldVector<ct,dim>&
          nblocal = Dune::ReferenceElements<ct,dim>::general(nbgt).position(0,0);
          double nbvolume = outside->geometry().integrationElement(nblocal)
                            *Dune::ReferenceElements<ct,dim>::general(nbgt).volume();
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
