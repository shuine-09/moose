#pragma once

#include "PointValueAtXFEMInterface.h"
//#include "DiscreteElementUserObject.h"
#include "Postprocessor.h"

class GradValueAtXFEMInterfacePostprocessor : public Postprocessor//,public DiscreteElementUserObject
{
public:
  static InputParameters validParams();

  GradValueAtXFEMInterfacePostprocessor(const InputParameters & parameters);

  virtual void initialize() override;

  virtual void getGradientValue(unsigned int point_id);

  virtual Real getValue() override ;

protected:
  /// Pointer to PointValueAtXFEMInterface object
  const PointValueAtXFEMInterface * _value_at_interface_uo;

  /// Value to indicate which side of the interface we want the gradient of (+1 or -1)
  const Real _side;

  Real _grad_value;
};
