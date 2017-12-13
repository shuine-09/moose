/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This userobject calculates the configurational force
//
#include "XFEMMeanStress.h"
#include "libmesh/fe_interface.h"
#include "DisplacedProblem.h"
#include "XFEM.h"
#include "libmesh/quadrature.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<XFEMMeanStress>()
{
  InputParameters params = validParams<ElementUserObject>();
  params += validParams<MaterialTensorCalculator>();
  params.addRequiredParam<std::string>("tensor", "The material tensor name.");
  params.addParam<Real>("radius", "Radius for average out stress");
  params.addParam<Real>("weibull_radius", "Radius for weibull distrbution");
  params.addParam<Real>("critical_stress", 0.0, "Critical stress.");
  params.addParam<bool>("use_weibull", false,"Use weibull distribution to propagate crack?");
  params.addParam<PostprocessorName>("average_h", "Postprocessor that gives average element size");
  return params;
}

XFEMMeanStress::XFEMMeanStress(const InputParameters & parameters):
    ElementUserObject(parameters),
    _material_tensor_calculator(parameters),
    _tensor(getMaterialProperty<SymmTensor>(getParam<std::string>("tensor"))),
    _critical_stress(getParam<Real>("critical_stress")),
    _weibull(getMaterialProperty<Real>("weibull")),
    _use_weibull(getParam<bool>("use_weibull")),
    _postprocessor( isParamValid("average_h") ? &getPostprocessorValue("average_h") : NULL )
{

  _fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);
  if (_fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblemBase in XFEMMarkerUserObject");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(_fe_problem->getXFEM());

  if (isParamValid("radius"))
  {
    _radius = getParam<Real>("radius");
  }
  else
    mooseError("MeanStress error: must set radius.");

  if (isParamValid("weibull_radius"))
  {
    _weibull_radius = getParam<Real>("weibull_radius");
  }
  else
    mooseError("MeanStress error: must set weibull radius.");

}

void
XFEMMeanStress::initialize()
{
  if (_postprocessor)
  {
    _radius = 3.0 * *_postprocessor;
    _weibull_radius = 3.0 * *_postprocessor;
  }

  _crack_front_points.clear();
  _elem_id_crack_tip.clear();
  _stress_tensor.clear();

  _num_crack_front_points = _xfem->numberCrackTips();

  _stress_tensor.resize(_num_crack_front_points*9);

  _weibull_at_tip.clear();
  _weibull_at_tip.resize(_num_crack_front_points);
  
  _weights.clear();
  _weights.resize(_num_crack_front_points);
  
  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {  
    _weibull_at_tip[i] = 9999e15;
    _weights[i] = 0;
  }

  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    _stress_tensor[i] = 0.0;

  _xfem->getCrackTipOrigin(_elem_id_crack_tip, _crack_front_points);

}

std::vector<Real>
XFEMMeanStress::getStressTensor()
{
  std::vector<Real> StressTensor(_num_crack_front_points*9);
  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    StressTensor[i] = 0.0;


  unsigned int numqp = _qrule->n_points();

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    Point crack_front = _crack_front_points[i];

    const Elem * undisplaced_elem  = NULL;
    if(_fe_problem->getDisplacedProblem() != NULL)
      undisplaced_elem = _fe_problem->getDisplacedProblem()->refMesh().elemPtr(_current_elem->id());
    else
      undisplaced_elem = _current_elem;
    
    for ( unsigned int qp = 0; qp < numqp; ++qp )
    {
      Real flag = _xfem->flagQpoint(undisplaced_elem, _q_point[qp]); //qp inside (flag = 1) or ouside (flag = 0) real domain
      Point dist_to_crack_front_vector = _q_point[qp] - crack_front;
      Real dist = std::pow(dist_to_crack_front_vector.size_sq(),0.5);
      Real fact = 1.0/(pow(2*libMesh::pi, 1.5) * pow(_radius, 3.0)) * std::exp(-0.5 * pow(dist/_radius,2.0));
      //Real fact = (1.0-dist/_radius);
      if (dist < 2.0 * _radius && flag > 0.5)
      {
        StressTensor[i*9+0] += _tensor[qp](0,0) * fact;
        StressTensor[i*9+1] += _tensor[qp](0,1) * fact;
        StressTensor[i*9+2] += _tensor[qp](0,2) * fact;
        StressTensor[i*9+3] += _tensor[qp](1,0) * fact;
        StressTensor[i*9+4] += _tensor[qp](1,1) * fact;
        StressTensor[i*9+5] += _tensor[qp](1,2) * fact;
        StressTensor[i*9+6] += _tensor[qp](2,0) * fact;
        StressTensor[i*9+7] += _tensor[qp](2,1) * fact;
        StressTensor[i*9+8] += _tensor[qp](2,2) * fact;
        _weights[i] += fact;
        
        if (dist < _weibull_radius && flag > 1.5)
        {
          if (_weibull[qp] < _weibull_at_tip[i])
          {
            //std::cout << "weillbull at tip = " << _weibull_at_tip[i] << ", weibull_eta[qp] = " << _weibull[qp] << std::endl;
            _weibull_at_tip[i] = _weibull[qp];
          }
        }
      }
    }
  }
  return StressTensor;
}

void 
XFEMMeanStress::execute()
{
  std::vector<Real> StressTensor = getStressTensor();
  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    _stress_tensor[i] += StressTensor[i];

//  for (unsigned int i = 0; i < _num_crack_front_points; i++)
//  {
//    if (_current_elem == _elem_id_crack_tip[i])
//    {
//      _weibull_at_tip[i] = _weibull[0];
//      break;
//    }
//  }
}

void
XFEMMeanStress::threadJoin(const UserObject & y)
{
  const XFEMMeanStress & pps = static_cast<const XFEMMeanStress &>(y);
  for (unsigned int i = 0; i < _num_crack_front_points*9; i++)
    _stress_tensor[i] += pps._stress_tensor[i];
}

void 
XFEMMeanStress::finalize()
{
  //_xfem->clearCrackGrowthDirection();
  _xfem->clearDoesCrackGrowth();

  gatherSum(_stress_tensor);
  gatherMin(_weibull_at_tip);
  gatherSum(_weights);

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    RealVectorValue direction;
    SymmTensor average_tensor;
    //std::cout << "total weight = " << _weights[i] << std::endl;
    average_tensor(0,0) = _stress_tensor[i*9+0] / _weights[i];
    average_tensor(0,1) = _stress_tensor[i*9+1] / _weights[i];
    average_tensor(0,2) = _stress_tensor[i*9+2] / _weights[i];
    average_tensor(1,0) = _stress_tensor[i*9+3] / _weights[i];
    average_tensor(1,1) = _stress_tensor[i*9+4] / _weights[i];
    average_tensor(1,2) = _stress_tensor[i*9+5] / _weights[i];
    average_tensor(2,0) = _stress_tensor[i*9+6] / _weights[i];
    average_tensor(2,1) = _stress_tensor[i*9+7] / _weights[i];
    average_tensor(2,2) = _stress_tensor[i*9+8] / _weights[i];
    Real tensor_quantity = _material_tensor_calculator.getTensorQuantity(average_tensor,_q_point[0],direction);
    direction /= pow(direction.size_sq(),0.5);
    Point normal(0.0,0.0,0.0);
    normal(0) = -direction(1);
    normal(1) = direction(0);
    bool does_elem_crack = false;
    if (_use_weibull)
      does_elem_crack = (tensor_quantity > _critical_stress * _weibull_at_tip[i]);
    else
      does_elem_crack = (tensor_quantity > _critical_stress);

    //std::cout << "stress = " << tensor_quantity << ", weibull_tip = " << _weibull_at_tip[i] << " does elem crack = " << does_elem_crack << std::endl;
    _xfem->updateDoesCrackGrowth(_elem_id_crack_tip[i], does_elem_crack);

    //_xfem->updateCrackGrowthDirection(_elem_id_crack_tip[i], normal);
    //std::cout << "crack front index (" << i << ") : average stress direction  = " << normal << std::endl; 
  } 
}

