/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "RankTwoAuxXFEM.h"
#include "RankTwoScalarTools.h"
#include "XFEM.h"

#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/quadrature.h"


template<>
InputParameters validParams<RankTwoAuxXFEM>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Access a component of a RankTwoTensor");
  params.addRequiredParam<MaterialPropertyName>("rank_two_tensor", "The rank two material tensor name");
  params.addRequiredRangeCheckedParam<unsigned int>("index_i", "index_i >= 0 & index_i <= 2", "The index i of ij for the tensor to output (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>("index_j", "index_j >= 0 & index_j <= 2", "The index j of ij for the tensor to output (0, 1, 2)");
  return params;
}

RankTwoAuxXFEM::RankTwoAuxXFEM(const InputParameters & parameters) :
    AuxKernel(parameters),
    _tensor(getMaterialProperty<RankTwoTensor>("rank_two_tensor")),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j"))
{
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);

  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in RankTwoAuxXFEM");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in RankTwoAuxXFEM");

}

void
RankTwoAuxXFEM::compute()
{
  Point tip_edge((_current_elem->point(0))(0), 0.0, 0);
  Point tip(0.5, 0.0, 0);
  if(_current_elem->contains_point(tip))
  {
    std::vector<Point> intersectionPoints;

    intersectionPoints.push_back(_current_elem->point(0));
    intersectionPoints.push_back(_current_elem->point(1));
    intersectionPoints.push_back(_current_elem->point(2));
    intersectionPoints.push_back(_current_elem->point(3));
    intersectionPoints.push_back(tip_edge);

    std::vector<Point> q_points;
    std::vector<Real> weights;

    _xfem->getXFEMqRuleOnSurface(intersectionPoints, tip, q_points, weights);

    FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);

    fe_problem->reinitElemPhys(_current_elem, q_points, 0);

    //fe_problem->prepareShapes(_var.number(), 0);

    fe_problem->reinitMaterials(_current_elem->subdomain_id(), 0, false);

    _n_local_dofs = _var.numberOfDofs();

    if (_n_local_dofs==1)  /* p0 */
    {
      Real value = 0;
      for (_qp=0; _qp < q_points.size(); _qp++)
        value += weights[_qp] * computeValue();
      value /= (_bnd ? _current_side_volume : _current_elem_volume);
      // update the variable data refernced by other kernels.
      // Note that this will update the values at the quadrature points too
      // (because this is an Elemental variable)
      _var.setNodalValue(value);
    }
    else                   /* high-order */
    {
      _local_re.resize(_n_local_dofs);
      _local_re.zero();
      _local_ke.resize(_n_local_dofs, _n_local_dofs);
      _local_ke.zero();

      // assemble the local mass matrix and the load
      for (unsigned int i = 0; i < _test.size(); i++)
        for (_qp = 0; _qp < q_points.size(); _qp++)
        {
          Real t = weights[_qp] * _test[i][_qp];
          _local_re(i) += t * computeValue();
          for (unsigned int j = 0; j < _test.size(); j++)
            _local_ke(i, j) += t * _test[j][_qp];
        }

      // mass matrix is always SPD
      _local_sol.resize(_n_local_dofs);
      _local_ke.cholesky_solve(_local_re, _local_sol);

      _var.setNodalValue(_local_sol);
    }

  }
  else
    AuxKernel::compute();
}

Real
RankTwoAuxXFEM::computeValue()
{
  return RankTwoScalarTools::component(_tensor[_qp], _i, _j);
}
