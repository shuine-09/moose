/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MaterialTensorIntegralXFEM.h"
#include "RankTwoScalarTools.h"
#include "XFEM.h"
#include "FEProblem.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<MaterialTensorIntegralXFEM>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredParam<MaterialPropertyName>("rank_two_tensor",
                                                "The rank two material tensor name");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i",
      "index_i >= 0 & index_i <= 2",
      "The index i of ij for the tensor to output (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j",
      "index_j >= 0 & index_j <= 2",
      "The index j of ij for the tensor to output (0, 1, 2)");
  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<bool>("normalized_error", false, "Calculate the normalized value.");
  return params;
}

MaterialTensorIntegralXFEM::MaterialTensorIntegralXFEM(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _tensor(getMaterialProperty<RankTwoTensor>("rank_two_tensor")),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j")),
    _normalized_error(getParam<bool>("normalized_error"))
{
  // FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);

  // if (fe_problem == NULL)
  //  mooseError("Problem casting _subproblem to FEProblem in MaterialTensorIntegralXFEM");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(_fe_problem.getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in MaterialTensorIntegralXFEM");
}

Real
MaterialTensorIntegralXFEM::getValue()
{
  gatherSum(_integral_value);
  return std::sqrt(_integral_value);
}

Real
MaterialTensorIntegralXFEM::computeIntegral()
{
  Real sum = 0;

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    Point crack_tip(0.5, 1.0, 0); // crack tip is at (0.5, 0.5, 0)
    Point q_pt = _q_point[_qp];

    Real x_to_tip = q_pt(0) - crack_tip(0);
    Real y_to_tip = q_pt(1) - crack_tip(1);

    Real alpha = 0.0;

    Real x_local = std::cos(alpha) * x_to_tip + std::sin(alpha) * y_to_tip;
    Real y_local = -std::sin(alpha) * x_to_tip + std::cos(alpha) * y_to_tip;

    Real r = std::sqrt(x_local * x_local + y_local * y_local);

    if (r < 0.0001)
      r = 0.0001;

    Real theta = std::atan2(y_local, x_local);

    Real KI = 3.562545435319310;
    Real sigma_xx_exact = KI / std::sqrt(2.0 * r * libMesh::pi) * std::cos(0.5 * theta) *
                          (1.0 - std::sin(0.5 * theta) * std::sin(1.5 * theta));
    Real sigma_yy_exact = KI / std::sqrt(2.0 * r * libMesh::pi) * std::cos(0.5 * theta) *
                          (1.0 + std::sin(0.5 * theta) * std::sin(1.5 * theta));
    Real sigma_xy_exact = KI / std::sqrt(2.0 * r * libMesh::pi) * std::cos(0.5 * theta) *
                          std::sin(0.5 * theta) * std::cos(1.5 * theta);

    Real error = 0.0;
    if (!_normalized_error)
    {
      if (_i == 0 && _j == 0)
      {
        error =
            std::pow((sigma_xx_exact - RankTwoScalarTools::component(_tensor[_qp], _i, _j)), 2.0);
      }
      else if ((_i == 1 && _j == 0) || (_i == 0 && _j == 1))
      {
        error =
            std::pow((sigma_xy_exact - RankTwoScalarTools::component(_tensor[_qp], _i, _j)), 2.0);
      }
      else if (_i == 1 & _j == 1)
      {
        error =
            std::pow((sigma_yy_exact - RankTwoScalarTools::component(_tensor[_qp], _i, _j)), 2.0);
      }
    }
    else
    {
      if (_i == 0 && _j == 0)
      {
        error = std::pow(sigma_xx_exact, 2.0);
      }
      else if ((_i == 1 && _j == 0) || (_i == 0 && _j == 1))
      {
        error = std::pow(sigma_xy_exact, 2.0);
      }
      else if (_i == 1 & _j == 1)
      {
        error = std::pow(sigma_yy_exact, 2.0);
      }
    }

    if (r > 0.2)
      error = 0.0;

    sum += _JxW[_qp] * _coord[_qp] * error;
  }
  return sum;
}

Real
MaterialTensorIntegralXFEM::computeQpIntegral()
{
  return RankTwoScalarTools::component(_tensor[_qp], _i, _j);
}
