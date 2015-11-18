/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This userobject calculates the configurational force
//
#include "XFEMMeanDirection.h"
#include "libmesh/fe_interface.h"
#include "DisplacedProblem.h"

template<>
InputParameters validParams<XFEMMeanDirection>()
{
  InputParameters params = validParams<ElementUserObject>();
  params += validParams<MaterialTensorCalculator>();
  params.addRequiredParam<std::string>("tensor", "The material tensor name.");
  params.addParam<Real>("radius", "Radius for average out stress");
  return params;
}

XFEMMeanDirection::XFEMMeanDirection(const InputParameters & parameters):
    ElementUserObject(parameters),
    _material_tensor_calculator(parameters),
    _tensor(getMaterialProperty<SymmTensor>(getParam<std::string>("tensor")))
{

  _fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (_fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMMarkerUserObject");
  _xfem = _fe_problem->get_xfem();

  if (isParamValid("radius"))
  {
    _radius = getParam<Real>("radius");
  }
  else
    mooseError("MeanStress error: must set radius.");
}

void
XFEMMeanDirection::initialize()
{ 
  _crack_front_points.clear();
  _elem_id_crack_tip.clear();
  _directions.clear();

  _num_crack_front_points = _xfem->num_crack_tips();

  _directions.resize(_num_crack_front_points*3);
  
  for (unsigned int i = 0; i < _num_crack_front_points*3; i++)
    _directions[i] = 0.0;

  _xfem->get_crack_tip_origin(_elem_id_crack_tip, _crack_front_points);

}

std::vector<Real>
XFEMMeanDirection::getDirection()
{
  std::vector<Real> Direction(_num_crack_front_points*3);
  for (unsigned int i = 0; i < _num_crack_front_points*3; i++)
    Direction[i] = 0.0;

  unsigned int numqp = _qrule->n_points();

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    Point crack_front = _crack_front_points[i];

    const Elem * undisplaced_elem  = NULL;
    if(_fe_problem->getDisplacedProblem() != NULL)
      undisplaced_elem = _fe_problem->getDisplacedProblem()->refMesh().elem(_current_elem->id());
    else
      undisplaced_elem = _current_elem;
    
    unsigned int count_numqp = 0; 
    SymmTensor average_tensor;
    Real average_dist = 0.0;
    RealVectorValue direction(0.0,0.0,0.0);
    
    for ( unsigned int qp = 0; qp < numqp; ++qp )
    {
      Real flag = _xfem->flag_qp_inside(undisplaced_elem, _q_point[qp]); //qp inside (flag = 1) or ouside (flag = 0) real domain
      Point dist_to_crack_front_vector = _q_point[qp] - crack_front;
      Real dist = std::pow(dist_to_crack_front_vector.size_sq(),0.5);
      if (dist < _radius && flag > 0.5)
      {
        average_tensor += _tensor[qp];
        average_dist += dist;
        count_numqp++;
      }
    }
    if (count_numqp)
    {
      average_tensor *= 1.0/(Real)count_numqp;
      average_dist *= 1.0/(Real)count_numqp;
      Real tensor_quantity = _material_tensor_calculator.getTensorQuantity(average_tensor,&_q_point[0],direction);
      direction *= (1.0 - average_dist/_radius);
    }

    Direction[i*3+0] = direction(0);
    Direction[i*3+1] = direction(1);
    Direction[i*3+2] = direction(2);
  }
  return Direction;
}

void 
XFEMMeanDirection::execute()
{
  std::vector<Real> Direction = getDirection();
  for (unsigned int i = 0; i < _num_crack_front_points*3; i++)
    _directions[i] += Direction[i];
}

void
XFEMMeanDirection::threadJoin(const UserObject & y)
{
  const XFEMMeanDirection & pps = static_cast<const XFEMMeanDirection &>(y);
  for (unsigned int i = 0; i < _num_crack_front_points*3; i++)
    _directions[i] += pps._directions[i];
}

void 
XFEMMeanDirection::finalize()
{
  _xfem->clear_crack_propagation_direction();
  gatherSum(_directions);

  for (unsigned int i = 0; i < _num_crack_front_points; i++)
  {
    RealVectorValue direction;
    direction(0) = _directions[i*3+0];
    direction(1) = _directions[i*3+1];
    direction(2) = _directions[i*3+2];
    std::cout << "crack tip = " << _crack_front_points[i] << std::endl;
    std::cout << "direction = " << direction << std::endl;
    if (direction.size_sq() > 1.0e-10)
      direction /= pow(direction.size_sq(),0.5);
    Point normal(0.0,0.0,0.0);
    normal(0) = direction(1);
    normal(1) = -direction(0);
    _xfem->update_crack_propagation_direction(_elem_id_crack_tip[i], normal);
    std::cout << "crack front index (" << i << ") : average direction  = " << normal << std::endl; 
  } 
}

