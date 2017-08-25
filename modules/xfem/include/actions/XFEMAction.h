/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMACTION_H
#define XFEMACTION_H

#include "Action.h"
#include "UserObjectInterface.h"

class XFEMAction;

template <>
InputParameters validParams<XFEMAction>();

class XFEMAction : public Action
{
public:
  XFEMAction(InputParameters params);

  virtual void act();

protected:
  std::vector<UserObjectName> _geom_cut_userobjects;
  std::string _xfem_qrule;
  std::string _order;
  std::string _family;
  bool _xfem_cut_plane;
  bool _xfem_use_crack_growth_increment;
  Real _xfem_crack_growth_increment;
  bool _use_crack_tip_enrichment;
  Real _crack_tip_enrichment_cutoff_radius;
  UserObjectName _crack_front_definition;
  std::vector<NonlinearVariableName> _enrich_displacement;
};

#endif // XFEMACTION_H
