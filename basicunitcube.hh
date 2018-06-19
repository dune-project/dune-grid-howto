// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef  BASICUNITCUBE_HH
#define  BASICUNITCUBE_HH

#include <dune/geometry/type.hh>
#include <dune/grid/common/gridfactory.hh>

// declaration of a basic unit cube that uses the GridFactory
template< int dim >
class BasicUnitCube;

// unit cube in two dimensions with 2 variants: triangle and rectangle elements
template<>
class BasicUnitCube< 2 >
{
protected:
  template< class Grid >
  static void insertVertices ( Dune::GridFactory< Grid > &factory )
  {
    Dune::FieldVector<double,2> pos;

    pos[0] = 0;  pos[1] = 0;
    factory.insertVertex(pos);                         /*@\label{uc:iv}@*/

    pos[0] = 1;  pos[1] = 0;
    factory.insertVertex(pos);

    pos[0] = 0;  pos[1] = 1;
    factory.insertVertex(pos);

    pos[0] = 1;  pos[1] = 1;
    factory.insertVertex(pos);
  }

  template< class Grid >
  static void insertSimplices ( Dune::GridFactory< Grid > &factory )
  {
    const Dune::GeometryType type( Dune::GeometryTypes::triangle );
    std::vector< unsigned int > cornerIDs( 3 );

    cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;
    factory.insertElement( type, cornerIDs );          /*@\label{uc:ie}@*/

    cornerIDs[0] = 2;  cornerIDs[1] = 1;  cornerIDs[2] = 3;
    factory.insertElement( type, cornerIDs );
  }

  template< class Grid >
  static void insertCubes ( Dune::GridFactory< Grid > &factory )
  {
    const Dune::GeometryType type( Dune::GeometryTypes::quadrilateral );
    std::vector< unsigned int > cornerIDs( 4 );
    for( int i = 0; i < 4; ++i )
      cornerIDs[ i ] = i;
    factory.insertElement( type, cornerIDs );
  }
};

// unit cube in 3 dimensions with two variants: tetraheda and hexahedra
template<>
class BasicUnitCube< 3 >
{
protected:
  template< class Grid >
  static void insertVertices ( Dune::GridFactory< Grid > &factory )
  {
    Dune::FieldVector< double, 3 > pos;

    pos[0] = 0;  pos[1] = 0;  pos[2] = 0;    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 0;  pos[2] = 0;    factory.insertVertex(pos);
    pos[0] = 0;  pos[1] = 1;  pos[2] = 0;    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 1;  pos[2] = 0;    factory.insertVertex(pos);
    pos[0] = 0;  pos[1] = 0;  pos[2] = 1;    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 0;  pos[2] = 1;    factory.insertVertex(pos);
    pos[0] = 0;  pos[1] = 1;  pos[2] = 1;    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 1;  pos[2] = 1;    factory.insertVertex(pos);
  }

  template< class Grid >
  static void insertSimplices ( Dune::GridFactory< Grid > &factory )
  {
    const Dune::GeometryType type( Dune::GeometryTypes::tetrahedron );
    std::vector< unsigned int > cornerIDs( 4 );

    cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;  cornerIDs[3] = 4;
    factory.insertElement( type, cornerIDs );

    cornerIDs[0] = 1;  cornerIDs[1] = 3;  cornerIDs[2] = 2;  cornerIDs[3] = 7;
    factory.insertElement( type, cornerIDs );

    cornerIDs[0] = 1;  cornerIDs[1] = 7;  cornerIDs[2] = 2;  cornerIDs[3] = 4;
    factory.insertElement( type, cornerIDs );

    cornerIDs[0] = 1;  cornerIDs[1] = 7;  cornerIDs[2] = 4;  cornerIDs[3] = 5;
    factory.insertElement( type, cornerIDs );

    cornerIDs[0] = 4;  cornerIDs[1] = 7;  cornerIDs[2] = 2;  cornerIDs[3] = 6;
    factory.insertElement( type, cornerIDs );
  }

  template< class Grid >
  static void insertCubes ( Dune::GridFactory< Grid > &factory )
  {
    const Dune::GeometryType type( Dune::GeometryTypes::hexahedron );
    std::vector< unsigned int > cornerIDs( 8 );
    for( int i = 0; i < 8; ++i )
      cornerIDs[ i ] = i;
    factory.insertElement( type, cornerIDs );
  }
};

#endif  /*BASICUNITCUBE_HH*/
