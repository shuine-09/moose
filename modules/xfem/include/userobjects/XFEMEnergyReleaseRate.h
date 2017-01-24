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

#ifndef XFEMENERGYRELEASERATE_H
#define XFEMENERGYRELEASERATE_H

#include "ElementUserObject.h"
#include "libmesh/string_to_enum.h"

// libMesh
#include "libmesh/point.h"
#include "libmesh/vector_value.h"
#include "libmesh/elem.h"

class XFEM;

class XFEMEnergyReleaseRate : public ElementUserObject
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  XFEMEnergyReleaseRate(const InputParameters & parameters);

  virtual ~XFEMEnergyReleaseRate() {}

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject &y);
  virtual void finalize();

protected:
  virtual std::vector<Real> computeIntegrals();
  virtual std::vector<Real> computeQpIntegrals(const std::vector<std::vector<Real> > & N_shape_func, const std::vector<std::vector<RealGradient> > & dN_shape_func);
  Real calcQValue(Point & node, Point & crack_front);  
  unsigned int _crack_front_point_index;
  const MaterialProperty<ColumnMajorMatrix> & _Eshelby_tensor;
  const MaterialProperty<RealVectorValue> * _J_thermal_term_vec;
  const MaterialProperty<Real> & _weibull;

private:
  unsigned int _qp;
  std::vector<Real> _integral_values;
  Real _radius_inner;
  Real _radius_outer;
  Real _critical_energy_release_rate;
  MooseMesh & _mesh;
  MooseSharedPointer<XFEM> _xfem;
  std::map<unsigned int, const Elem* > _elem_id_crack_tip;
  std::vector<Point> _crack_front_points;
  std::vector<Point> _crack_directions;
  unsigned int _num_crack_front_points;
  bool _use_weibull;
  std::vector<Real> _weibull_at_tip;
  const PostprocessorValue * const _postprocessor;
};

template<>
InputParameters validParams<XFEMEnergyReleaseRate>();

#endif //XFEMENERGYRELEASERATE_H
