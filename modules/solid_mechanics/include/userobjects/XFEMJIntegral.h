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

#ifndef XFEMJINTEGRAL_H
#define XFEMJINTEGRAL_H

#include "ElementUserObject.h"
#include "CrackFrontDefinition.h"
#include "libmesh/string_to_enum.h"

// libMesh
#include "libmesh/point.h"
#include "libmesh/vector_value.h"
#include "libmesh/elem.h"

class XFEM;

class XFEMJIntegral : public ElementUserObject
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  XFEMJIntegral(const InputParameters & parameters);

  virtual ~XFEMJIntegral() {}

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject &y);
  virtual void finalize();

protected:
  virtual DenseVector<Real> computeIntegrals();
  virtual DenseVector<Real> computeQpIntegrals(const std::vector<std::vector<Real> > & N_shape_func, const std::vector<std::vector<RealGradient> > & dN_shape_func);
  void  projectToFrontAtPoint(Real & dist_to_front, Real & dist_along_tangent, Point p);
  Real calcQValue(Node* node);  
  const CrackFrontDefinition * const _crack_front_definition;
  bool _has_crack_front_point_index;
  unsigned int _crack_front_point_index;
  bool _treat_as_2d;
  const MaterialProperty<ColumnMajorMatrix> & _Eshelby_tensor;
  const MaterialProperty<RealVectorValue> * _J_thermal_term_vec;
  bool _convert_J_to_K;
  bool _has_symmetry_plane;
  Real _poissons_ratio;
  Real _youngs_modulus; 

private:
  unsigned int _qp;
  DenseVector<Real> _integral_values;
  Real _radius_inner;
  Real _radius_outer;
  MooseMesh & _mesh;
  XFEM *_xfem;
  std::map<unsigned int, Real> _map_jintegral;
  unsigned int _num_crack_front_points;
};

template<>
InputParameters validParams<XFEMJIntegral>();

#endif //XFEMJINTEGRAL_H
