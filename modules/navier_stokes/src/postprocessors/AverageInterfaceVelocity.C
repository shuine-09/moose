#include "AverageInterfaceVelocity.h"

registerMooseObject("NavierStokesApp", AverageInterfaceVelocity);

InputParameters
AverageInterfaceVelocity::validParams()
{
  InputParameters params = ElementPostprocessor::validParams();
  params.addClassDescription("Compute average interface velocity.");
  params.addRequiredCoupledVar("velocity", "Velocity vector variable.");
  return params;
}

AverageInterfaceVelocity::AverageInterfaceVelocity(const InputParameters & parameters)
  : ElementPostprocessor(parameters),
    _delta_function(getADMaterialProperty<Real>("delta_function")),
    _velocity(adCoupledVectorValue("velocity"))
{
}

void
AverageInterfaceVelocity::initialize()
{
  // Compute maximum velocity
  _max_velocity = std::numeric_limits<Real>::min();
}

void
AverageInterfaceVelocity::execute()
{
  for (unsigned int qp = 0; qp < _q_point.size(); ++qp)
  {
    if (_delta_function[qp] > 0.1)
    {
      RealVectorValue vel = MetaPhysicL::raw_value(_velocity[qp]);
      _max_velocity = std::max(_max_velocity, std::abs(vel.norm()));
    }
  }
}

void
AverageInterfaceVelocity::finalize()
{
  gatherMax(_max_velocity);
}

void
AverageInterfaceVelocity::threadJoin(const UserObject & user_object)
{
  const AverageInterfaceVelocity & pps = static_cast<const AverageInterfaceVelocity &>(user_object);
  _max_velocity = std::max(_max_velocity, pps._max_velocity);
}

PostprocessorValue
AverageInterfaceVelocity::getValue()
{
  return _max_velocity;
}
