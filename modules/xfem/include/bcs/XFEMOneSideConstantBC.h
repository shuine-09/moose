//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMONESIDECONSTANTBC_H
#define XFEMONESIDECONSTANTBC_H

#include "PresetBC.h"
#include "CrackFrontDefinition.h"

class XFEMOneSideConstantBC;

template <>
InputParameters validParams<XFEMOneSideConstantBC>();

/**
 * XFEMOneSideConstantBC is used in XFEM Crack Tip Enrichment to fix DOFs to zero for those
 * nodes with basis function supports that are far away from any crack tip.
 */
class XFEMOneSideConstantBC : public PresetBC
{
public:
  XFEMOneSideConstantBC(const InputParameters & parameters);

protected:
  virtual bool shouldApply() override;

  const Real _left_x;
};

#endif /* XFEMONESIDECONSTANTBC_H */
