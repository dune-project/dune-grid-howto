// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// a smooth function
template<typename ct, int dim>
class Exp {
public:
  Exp () {midpoint = 0.5;}
  double operator() (const Dune::FieldVector<ct,dim>& x) const
  {
    Dune::FieldVector<ct,dim> y(x);
    y -= midpoint;
    return exp(-3.234*(y*y));
  }
private:
  Dune::FieldVector<ct,dim> midpoint;
};

// a function with a local feature
template<typename ct, int dim>
class Needle {
public:
  Needle ()
  {
    midpoint = 0.5;
    midpoint[dim-1] = 1;
  }
  double operator() (const Dune::FieldVector<ct,dim>& x) const
  {
    Dune::FieldVector<ct,dim> y(x);
    y -= midpoint;
    return 1.0/(1E-4+y*y);
  }
private:
  Dune::FieldVector<ct,dim> midpoint;
};
