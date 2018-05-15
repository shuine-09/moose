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
  // Get input parameters from parent class
  InputParameters params = validParams<GeometricCut2DUserObject>();

  params.addRequiredParam<UserObjectName>("interface_value_uo", "XXX");
  // Add required parameters
  params.addRequiredParam<std::vector<Real>>("cut_data",
                                             "Vector of Real values providing cut information");
  // Add optional parameters
  params.addParam<std::vector<Real>>("cut_scale", "X,Y scale factors for geometric cuts");
  params.addParam<std::vector<Real>>("cut_translate", "X,Y translations for geometric cuts");
  params.addRequiredParam<VariableName>(
      "var", "The name of solution variable used to calcuate interface velocity.");
  // Class description
  params.addClassDescription("Creates a UserObject for a line segment cut on 2D meshes for XFEM");
  // Return the parameters
  return params;
}

MovingLineSegmentCutSetUserObject::MovingLineSegmentCutSetUserObject(
    const InputParameters & parameters)
  : GeometricCut2DUserObject(parameters),
    VectorPostprocessorInterface(this),
    _cut_data(getParam<std::vector<Real>>("cut_data")),
    _var_number(_subproblem.getVariable(_tid, parameters.get<VariableName>("var")).number()),
    _system(_subproblem.getSystem(getParam<VariableName>("var"))),
    _solution(_system.current_local_solution.get())
{
  // FEProblemBase * fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);

  // const UserObject * uo =
  //     &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("interface_value_uo")));
  //
  // if (dynamic_cast<const PointValueAtXFEMInterface *>(uo) == nullptr)
  //   mooseError(
  //       "UserObject casting to PointValueAtXFEMInterface in MovingLineSegmentCutSetUserObject");
  //
  // _interface_value_uo = dynamic_cast<const PointValueAtXFEMInterface *>(uo);

  // Set up constant parameters
  const int line_cut_data_len = 6;

  // Throw error if length of cut_data is incorrect
  if (_cut_data.size() % line_cut_data_len != 0)
    mooseError("Length of MovingLineSegmentCutSetUserObject cut_data must be a multiple of 6.");

  unsigned int num_cuts = _cut_data.size() / line_cut_data_len;

  // Assign scale and translate parameters
  std::pair<Real, Real> scale;
  if (isParamValid("cut_scale"))
  {
    auto vec_scale = getParam<std::vector<Real>>("cut_scale");
    scale = std::make_pair(vec_scale[0], vec_scale[1]);
  }
  else
  {
    scale = std::make_pair(1.0, 1.0);
  }

  std::pair<Real, Real> trans;
  if (isParamValid("cut_translate"))
  {
    auto vec_trans = getParam<std::vector<Real>>("cut_translate");
    trans = std::make_pair(vec_trans[0], vec_trans[1]);
  }
  else
  {
    trans = std::make_pair(0.0, 0.0);
  }

  // Clear _start_times & _end_times vectors initialized from
  // time_start_cut & time_end_cut values
  _cut_time_ranges.clear();

  for (unsigned int i = 0; i < num_cuts; ++i)
  {
    Real x0 = (_cut_data[i * line_cut_data_len + 0] + trans.first) * scale.first;
    Real y0 = (_cut_data[i * line_cut_data_len + 1] + trans.second) * scale.second;
    Real x1 = (_cut_data[i * line_cut_data_len + 2] + trans.first) * scale.first;
    Real y1 = (_cut_data[i * line_cut_data_len + 3] + trans.second) * scale.second;
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

void
MovingLineSegmentCutSetUserObject::execute()
{
  std::vector<Real> cut_data_copy = _cut_data;

  std::vector<Real> value_positive = _interface_value_uo->getValueAtPositiveLevelSet();
  std::vector<Real> value_negative = _interface_value_uo->getValueAtNegativeLevelSet();

  // Set up constant parameters
  const int line_cut_data_len = 6;

  // Throw error if length of cut_data is incorrect
  if (_cut_data.size() % line_cut_data_len != 0)
    mooseError("Length of MovingLineSegmentCutSetUserObject cut_data must be a multiple of 6.");

  unsigned int num_cuts = _cut_data.size() / line_cut_data_len;

  if (value_positive.size())
  {
    for (unsigned i = 0; i < value_positive.size() - 1; ++i)
    {
      // cut_data_copy[i * line_cut_data_len + 0] = _values[0] / 50 + 0.5;
      // cut_data_copy[i * line_cut_data_len + 2] = _values[0] / 50 + 0.5;
      cut_data_copy[i * 6 + 0] += value_positive[i] / 10;
      cut_data_copy[i * 6 + 2] += value_positive[i + 1] / 10;
    }
  }

  _cut_line_endpoints.clear();
  for (unsigned int i = 0; i < num_cuts; ++i)
  {
    Real x0 = cut_data_copy[i * line_cut_data_len + 0];
    Real y0 = cut_data_copy[i * line_cut_data_len + 1];
    Real x1 = cut_data_copy[i * line_cut_data_len + 2];
    Real y1 = cut_data_copy[i * line_cut_data_len + 3];
    std::cout << "x0 = " << x0 << ", y0 = " << y0 << ", x1 = " << x1 << ", y1 = " << y1
              << std::endl;
    _cut_line_endpoints.push_back(std::make_pair(Point(x0, y0, 0.0), Point(x1, y1, 0.0)));
  }

  GeometricCutUserObject::execute();
}

Real
MovingLineSegmentCutSetUserObject::getLocationX() const
{
  return _cut_data[0];
}

void
MovingLineSegmentCutSetUserObject::finalize()
{
  // {
  //   if (_values.size())
  //   {
  //     // cut_data_copy[i * line_cut_data_len + 0] = _values[0] / 50 + 0.5;
  //     // cut_data_copy[i * line_cut_data_len + 2] = _values[0] / 50 + 0.5;
  //     _cut_data[0] += _values[0] / 200;
  //     _cut_data[2] += _values[0] / 200;
  //   }
  //   else
  //   {
  //     _cut_data[0] = 0.5;
  //     _cut_data[2] = 0.5;
  //   }
  // }

  std::vector<Real> value_positive = _interface_value_uo->getValueAtPositiveLevelSet();
  std::vector<Real> value_negative = _interface_value_uo->getValueAtNegativeLevelSet();

  std::cout << "called " << std::endl;

  for (unsigned i = 0; i < value_positive.size(); ++i)
  {
    std::cout << "positive_size_value[" << i << "] = " << value_positive[i] << std::endl;
    std::cout << "negative_size_value[" << i << "] = " << value_negative[i] << std::endl;
  }

  {
    if (value_positive.size())
    {
      for (unsigned i = 0; i < value_positive.size() - 1; ++i)
      {
        // cut_data_copy[i * line_cut_data_len + 0] = _values[0] / 50 + 0.5;
        // cut_data_copy[i * line_cut_data_len + 2] = _values[0] / 50 + 0.5;
        _cut_data[i * 6 + 0] += value_positive[i] / 10;
        _cut_data[i * 6 + 2] += value_positive[i + 1] / 10;
      }
    }
    // else
    // {
    //   _cut_data[0] = 0.5;
    //   _cut_data[2] = 0.5;
    // }
  }

  GeometricCutUserObject::finalize();
}
