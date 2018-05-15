//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMTWOSIDEDIRICHLETNITSCHE_H
#define XFEMTWOSIDEDIRICHLETNITSCHE_H

// MOOSE includes
#include "XFEMTwoSideDirichlet.h"

// Forward Declarations
class XFEMTwoSideDirichletNitsche;

class XFEM;

template <>
InputParameters validParams<XFEMTwoSideDirichletNitsche>();

class XFEMTwoSideDirichletNitsche : public XFEMTwoSideDirichlet
{
public:
  XFEMTwoSideDirichletNitsche(const InputParameters & parameters);
  virtual ~XFEMTwoSideDirichletNitsche();

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;

  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  /// Diffusivity in positive level set side
  Real _diffusivity_at_positive_level_set_side;

  /// Diffusivity in negative level set side
  Real _diffusivity_at_negative_level_set_side;

  /// Pointer to the XFEM controller object
  std::shared_ptr<XFEM> _xfem;
};

#endif /* XFEMTWOSIDEDIRICHLETNITSCHE_H_ */
