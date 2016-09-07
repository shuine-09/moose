/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "XFEMMeshCutUserObject.h"

template<>
InputParameters validParams<XFEMMeshCutUserObject>()
{
  InputParameters params = validParams<DiscreteElementUserObject>();
  params.addClassDescription("XFEM mesh cut base class. Override the virtual functions in your class");
  return params;
}

XFEMMeshCutUserObject::XFEMMeshCutUserObject(const InputParameters & parameters) :
    DiscreteElementUserObject(parameters)
{
}
