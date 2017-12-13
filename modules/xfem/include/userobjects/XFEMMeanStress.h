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

#ifndef XFEMMEANSTRESS_H
#define XFEMMEANSTRESS_H

#include "ElementUserObject.h"
#include "libmesh/string_to_enum.h"
#include "MaterialTensorCalculator.h"

// libMesh
#include "libmesh/point.h"
#include "libmesh/vector_value.h"
#include "libmesh/elem.h"

class XFEM;

class XFEMMeanStress : public ElementUserObject
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  XFEMMeanStress(const InputParameters & parameters);

  virtual ~XFEMMeanStress() {}

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject &y);
  virtual void finalize();

protected:
  virtual std::vector<Real> getStressTensor();
  unsigned int _crack_front_point_index;

private:
  MaterialTensorCalculator _material_tensor_calculator;
  const MaterialProperty<SymmTensor> & _tensor;
  std::vector<Real> _stress_tensor;
  Real _radius;
  Real _weibull_radius;
  MooseSharedPointer<XFEM> _xfem;
  std::map<unsigned int, const Elem* > _elem_id_crack_tip;
  std::vector<Point> _crack_front_points;
  std::vector<Real> _weights;
  unsigned int _num_crack_front_points;
  FEProblemBase *_fe_problem;

  Real _critical_stress;
  const MaterialProperty<Real> & _weibull;
  bool _use_weibull;
  std::vector<Real> _weibull_at_tip;
  const PostprocessorValue * const _postprocessor;
};

template<>
InputParameters validParams<XFEMMeanStress>();

#endif //XFEMMEANSTRESS_H
