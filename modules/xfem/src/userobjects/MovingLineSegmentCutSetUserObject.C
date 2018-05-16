//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MovingLineSegmentCutSetUserObject.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "VectorPostprocessorInterface.h"
#include "VectorPostprocessor.h"
#include "PointValueAtXFEMInterface.h"

// MOOSE includes
#include "MooseError.h"

registerMooseObject("XFEMApp", MovingLineSegmentCutSetUserObject);

template <>
InputParameters
validParams<MovingLineSegmentCutSetUserObject>()
{
  InputParameters params = validParams<GeometricCut2DUserObject>();
  params.addRequiredParam<UserObjectName>("interface_value_uo", "XXX");
  params.addRequiredParam<std::vector<Real>>("cut_data",
                                             "Vector of Real values providing cut information");
  params.addRequiredParam<VariableName>(
      "var", "The name of solution variable used to calcuate interface velocity.");
  params.addRequiredParam<Real>("diffusivity_at_positive_level_set_side",
                                "Diffusivity at level set positive side.");
  params.addRequiredParam<Real>("diffusivity_at_negative_level_set_side",
                                "Diffusivity at level set negative side.");
  params.addClassDescription("Creates a UserObject for a line segment cut on 2D meshes for XFEM");
  return params;
}

MovingLineSegmentCutSetUserObject::MovingLineSegmentCutSetUserObject(
    const InputParameters & parameters)
  : GeometricCut2DUserObject(parameters),
    VectorPostprocessorInterface(this),
    _cut_data(getParam<std::vector<Real>>("cut_data")),
    _var_number(_subproblem.getVariable(_tid, parameters.get<VariableName>("var")).number()),
    _system(_subproblem.getSystem(getParam<VariableName>("var"))),
    _solution(_system.current_local_solution.get()),
    _diffusivity_at_positive_level_set_side(
        getParam<Real>("diffusivity_at_positive_level_set_side")),
    _diffusivity_at_negative_level_set_side(
        getParam<Real>("diffusivity_at_negative_level_set_side"))
{
  // Set up constant parameters
  const int line_cut_data_len = 6;

  // Throw error if length of cut_data is incorrect
  if (_cut_data.size() % line_cut_data_len != 0)
    mooseError("Length of MovingLineSegmentCutSetUserObject cut_data must be a multiple of 6.");

  unsigned int num_cuts = _cut_data.size() / line_cut_data_len;

  // Clear _start_times & _end_times vectors initialized from
  // time_start_cut & time_end_cut values
  _cut_time_ranges.clear();

  for (unsigned int i = 0; i < num_cuts; ++i)
  {
    Real x0 = _cut_data[i * line_cut_data_len + 0];
    Real y0 = _cut_data[i * line_cut_data_len + 1];
    Real x1 = _cut_data[i * line_cut_data_len + 2];
    Real y1 = _cut_data[i * line_cut_data_len + 3];
    _cut_line_endpoints.push_back(std::make_pair(Point(x0, y0, 0.0), Point(x1, y1, 0.0)));

    _cut_time_ranges.push_back(
        std::make_pair(_cut_data[i * line_cut_data_len + 4], _cut_data[i * line_cut_data_len + 5]));
  }

  if (_cut_line_endpoints.size() != _cut_time_ranges.size())
    mooseError("Number of start/end times must match number of cut line endpoint sets");
}

const std::vector<Point>
MovingLineSegmentCutSetUserObject::getCrackFrontPoints(
    unsigned int /*num_crack_front_points*/) const
{
  mooseError("getCrackFrontPoints() is not implemented for this object.");
}

void
MovingLineSegmentCutSetUserObject::initialize()
{
  FEProblemBase * fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);

  const UserObject * uo =
      &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("interface_value_uo")));

  if (dynamic_cast<const PointValueAtXFEMInterface *>(uo) == nullptr)
    mooseError(
        "UserObject casting to PointValueAtXFEMInterface in MovingLineSegmentCutSetUserObject");

  _interface_value_uo = dynamic_cast<const PointValueAtXFEMInterface *>(uo);
}

Real
MovingLineSegmentCutSetUserObject::calculateInterfaceVelocity(Real value_positive,
                                                              Real value_negative,
                                                              RealVectorValue grad_positive,
                                                              RealVectorValue grad_negative)
{
  if (!MooseUtils::absoluteFuzzyEqual(value_positive, value_negative))
    return std::abs((_diffusivity_at_positive_level_set_side * grad_positive(0) -
                     _diffusivity_at_negative_level_set_side * grad_negative(0)) /
                    (value_positive - value_negative));
  else
    return 0.0;
}

void
MovingLineSegmentCutSetUserObject::execute()
{
  std::vector<Real> cut_data_copy = _cut_data;

  std::map<unsigned int, Real> value_positive = _interface_value_uo->getValueAtPositiveLevelSet();
  std::map<unsigned int, Real> value_negative = _interface_value_uo->getValueAtNegativeLevelSet();
  std::map<unsigned int, RealVectorValue> grad_positive =
      _interface_value_uo->getGradientAtPositiveLevelSet();
  std::map<unsigned int, RealVectorValue> grad_negative =
      _interface_value_uo->getGradientAtNegativeLevelSet();

  // Set up constant parameters
  const int line_cut_data_len = 6;

  // Throw error if length of cut_data is incorrect
  if (_cut_data.size() % line_cut_data_len != 0)
    mooseError("Length of MovingLineSegmentCutSetUserObject cut_data must be a multiple of 6.");

  unsigned int num_cuts = _cut_data.size() / line_cut_data_len;

  for (unsigned i = 0; i < value_positive.size() - 1; ++i)
  {
    cut_data_copy[i * 6 + 0] +=
        calculateInterfaceVelocity(
            value_positive[i], value_negative[i], grad_positive[i], grad_negative[i]) *
        _dt;
    cut_data_copy[i * 6 + 2] += calculateInterfaceVelocity(value_positive[i + 1],
                                                           value_negative[i + 1],
                                                           grad_positive[i + 1],
                                                           grad_negative[i + 1]) *
                                _dt;
    ;
  }

  _cut_line_endpoints.clear();
  for (unsigned int i = 0; i < num_cuts; ++i)
  {
    Real x0 = cut_data_copy[i * line_cut_data_len + 0];
    Real y0 = cut_data_copy[i * line_cut_data_len + 1];
    Real x1 = cut_data_copy[i * line_cut_data_len + 2];
    Real y1 = cut_data_copy[i * line_cut_data_len + 3];
    _cut_line_endpoints.push_back(std::make_pair(Point(x0, y0, 0.0), Point(x1, y1, 0.0)));
  }

  GeometricCutUserObject::execute();
}

void
MovingLineSegmentCutSetUserObject::finalize()
{
  std::map<unsigned int, Real> value_positive = _interface_value_uo->getValueAtPositiveLevelSet();
  std::map<unsigned int, Real> value_negative = _interface_value_uo->getValueAtNegativeLevelSet();
  std::map<unsigned int, RealVectorValue> grad_positive =
      _interface_value_uo->getGradientAtPositiveLevelSet();
  std::map<unsigned int, RealVectorValue> grad_negative =
      _interface_value_uo->getGradientAtNegativeLevelSet();

  // for (unsigned i = 0; i < value_positive.size(); ++i)
  // {
  //   std::cout << "positive_size_value[" << i << "] = " << value_positive[i] << std::endl;
  //   std::cout << "negative_size_value[" << i << "] = " << value_negative[i] << std::endl;
  // }

  for (unsigned i = 0; i < value_positive.size() - 1; ++i)
  {
    _cut_data[i * 6 + 0] +=
        calculateInterfaceVelocity(
            value_positive[i], value_negative[i], grad_positive[i], grad_negative[i]) *
        _dt;
    _cut_data[i * 6 + 2] += calculateInterfaceVelocity(value_positive[i + 1],
                                                       value_negative[i + 1],
                                                       grad_positive[i + 1],
                                                       grad_negative[i + 1]) *
                            _dt;
  }

  GeometricCutUserObject::finalize();
}
