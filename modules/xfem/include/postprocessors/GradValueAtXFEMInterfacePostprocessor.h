#pragma once

#include "PointValueAtXFEMInterface.h"
#include "GeneralPostprocessor.h"



/**
*
* Retrieves the value of the gradient at a specified interface, at the
* specified side of that interface.
*
*/

class GradValueAtXFEMInterfacePostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  GradValueAtXFEMInterfacePostprocessor(const InputParameters & parameters);

  virtual void initialize() override;

  virtual void execute() override {}

  virtual Real getValue() override;

protected:
  /// Pointer to PointValueAtXFEMInterface object
  const PointValueAtXFEMInterface * _value_at_interface_uo;

  /// Value to indicate which side of the interface we want the gradient of (+1 or -1)
  const Real _side;

};
