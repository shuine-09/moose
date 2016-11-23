/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef CRACKFRONTRTHETA_H
#define CRACKFRONTRTHETA_H

#include "DiscreteElementUserObject.h"

// libMesh
//#include "libmesh/point.h"
//#include "libmesh/vector_value.h"
//#include "libmesh/elem.h"

#include "SymmTensor.h"

class CrackFrontRTheta : public DiscreteElementUserObject
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  CrackFrontRTheta(const InputParameters & parameters);

  virtual ~CrackFrontRTheta() {}

protected:
  void calcRTheta(Point & p, Point & crack_front_point, Point & crack_direction);
  RealVectorValue rotateToCrackFrontCoords(const RealVectorValue vector) const;
  RealVectorValue rotateToCrackFrontCoords(const Point point) const;
  ColumnMajorMatrix rotateToCrackFrontCoords(const SymmTensor tensor) const;
  ColumnMajorMatrix rotateToCrackFrontCoords(const ColumnMajorMatrix tensor) const;
  Real getRadius(){return _r;};
  Real getTheta(){return _theta;};

private:
  Real _r;
  Real _theta;
  ColumnMajorMatrix _rot_mat;
  RealVectorValue _crack_plane_normal;
};

template<>
InputParameters validParams<CrackFrontRTheta>();

#endif //CRACKFRONTRTHETA_H
