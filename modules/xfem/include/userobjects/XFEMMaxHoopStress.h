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

#ifndef XFEMMAXHOOPSTRESS_H
#define XFEMMAXHOOPSTRESS_H

#include "ElementUserObject.h"
#include "libmesh/string_to_enum.h"

// libMesh
#include "libmesh/point.h"
#include "libmesh/vector_value.h"
#include "libmesh/elem.h"

#include "SymmTensor.h"

class XFEM;

class XFEMMaxHoopStress : public ElementUserObject
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  XFEMMaxHoopStress(const InputParameters & parameters);

  virtual ~XFEMMaxHoopStress() {}

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject &y);
  virtual void finalize();

protected:
  virtual std::vector<Real> computeIntegrals();
  virtual std::vector<Real> computeQpIntegrals(const std::vector<std::vector<Real> > & N_shape_func, const std::vector<std::vector<RealGradient> > & dN_shape_func);
  Real calcQValue(Point & node, Point & crack_front);  
  void calcRTheta(Point & p, Point & crack_front_point, Point & crack_direction);

  RealVectorValue rotateToCrackFrontCoords(const RealVectorValue vector) const;
  RealVectorValue rotateToCrackFrontCoords(const Point point) const;
  ColumnMajorMatrix rotateToCrackFrontCoords(const SymmTensor tensor) const;
  ColumnMajorMatrix rotateToCrackFrontCoords(const ColumnMajorMatrix tensor) const;

  unsigned int _crack_front_point_index;
  
  const MaterialProperty<ColumnMajorMatrix> & _Eshelby_tensor;
  const MaterialProperty<RealVectorValue> * _J_thermal_term_vec;
  const MaterialProperty<SymmTensor> & _stress;
  const MaterialProperty<SymmTensor> & _strain;
  const VariableGradient & _grad_disp_x;
  const VariableGradient & _grad_disp_y;
  const VariableGradient & _grad_disp_z;
  const VariableGradient & _temp_grad;

  ColumnMajorMatrix  _aux_stress;
  ColumnMajorMatrix  _aux_disp;
  ColumnMajorMatrix  _aux_grad_disp;
  ColumnMajorMatrix  _aux_strain;

  enum SIF_MODE
  {
    KI,
    KII
  };

  void computeAuxFields(const SIF_MODE sif_mode, ColumnMajorMatrix & stress, ColumnMajorMatrix & disp, ColumnMajorMatrix & grad_disp, ColumnMajorMatrix & strain);

private:
  unsigned int _qp;
  std::vector<Real> _integral_values;
  Real _radius_inner;
  Real _radius_outer;
  Real _thermal_expansion;
  MooseMesh & _mesh;
  MooseSharedPointer<XFEM> _xfem;
  std::map<unsigned int, const Elem* > _elem_id_crack_tip;
  std::vector<Point> _crack_front_points;
  unsigned int _num_crack_front_points;
  std::vector<Point> _crack_front_directions;
  FEProblemBase * _fe_problem;

  std::vector<SIF_MODE> _sif_mode;

  Real _poissons_ratio;
  Real _youngs_modulus;
  Real _kappa;
  Real _shear_modulus;
  Real _K_factor;
  Real _r;
  Real _theta;

  ColumnMajorMatrix _rot_mat;
  RealVectorValue _crack_plane_normal;
};

template<>
InputParameters validParams<XFEMMaxHoopStress>();

#endif //XFEMMAXHOOPSTRESS_H
