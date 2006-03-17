// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
//! initialize the vector of unknowns with initial value
template<class G, class M, class V>
void initialize (const G& grid, const M& mapper, V& c)
{
  // first we extract the dimensions of the grid
  const int dim = G::dimension;
  const int dimworld = G::dimensionworld;

  // type used for coordinates in the grid
  typedef typename G::ctype ct;

  // leaf iterator type
  typedef typename G::template Codim<0>::LeafIterator LeafIterator;

  // iterate through leaf grid an evaluate c0 at cell center
  LeafIterator endit = grid.template leafend<0>();
  for (LeafIterator it = grid.template leafbegin<0>(); it!=endit; ++it)
  {
    // get geometry type
    Dune::GeometryType gt = it->geometry().type();

    // get cell center in reference element
    const Dune::FieldVector<ct,dim>&
    local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);

    // get global coordinate of cell center
    Dune::FieldVector<ct,dimworld> global =
      it->geometry().global(local);

    // initialize update with source term
    c[mapper.map(*it)] = c0(global);
  }
}
