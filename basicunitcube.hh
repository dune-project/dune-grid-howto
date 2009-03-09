// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef BASICUNITCUBE_HH
#define BASICUNITCUBE_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/gridfactory.hh>

template< int dim >
class BasicUnitCube;


template<>
class BasicUnitCube< 2 >
{
protected:
  template< class Grid >
  static void insertVertices ( Dune::GridFactory< Grid > &factory )
  {
    Dune::FieldVector<double,2> pos;

    pos[0] = 0;  pos[1] = 0;
    factory.insertVertex(pos);

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
    const Dune::GeometryType type( Dune::GeometryType::simplex, 2 );
    std::vector< unsigned int > cornerIDs( 3 );

    cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;
    factory.insertElement( type, cornerIDs );

    cornerIDs[0] = 2;  cornerIDs[1] = 1;  cornerIDs[2] = 3;
    factory.insertElement( type, cornerIDs );
  }

  template< class Grid >
  static void insertCubes ( Dune::GridFactory< Grid > &factory )
  {
    const Dune::GeometryType type( Dune::GeometryType::cube, 2 );
    std::vector< unsigned int > cornerIDs( 4 );
    for( int i = 0; i < 4; ++i )
      cornerIDs[ i ] = i;
    factory.insertElement( type, cornerIDs );
  }
};


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
    const Dune::GeometryType type( Dune::GeometryType::simplex, 3 );
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
    const Dune::GeometryType type( Dune::GeometryType::cube, 3 );
    std::vector< unsigned int > cornerIDs( 8 );
    for( int i = 0; i < 8; ++i )
      cornerIDs[ i ] = i;
    factory.insertElement( type, cornerIDs );
  }
};

#endif
