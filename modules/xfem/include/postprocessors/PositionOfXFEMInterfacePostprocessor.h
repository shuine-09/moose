#pragma once

#include "PointValueAtXFEMInterface.h"
#include "GeneralPostprocessor.h"

/**
*
* Retrieves the position of a specified interface
*
*/


class PositionOfXFEMInterfacePostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  PositionOfXFEMInterfacePostprocessor(const InputParameters & parameters);

  virtual void initialize() override;

  virtual void execute() override {}

  virtual Real getValue() override;

protected:
  /// Pointer to PointValueAtXFEMInterface object
  const PointValueAtXFEMInterface * _value_at_interface_uo;

};
