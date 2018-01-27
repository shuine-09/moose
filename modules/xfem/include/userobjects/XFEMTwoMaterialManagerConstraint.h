/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMTWOMATERIALMANAGERCONSTRAINT_H
#define XFEMTWOMATERIALMANAGERCONSTRAINT_H

#include "ElemElemConstraint.h"
#include "XFEMElemPairMaterialManager.h"

/**
 * ElemElemConstraint that utilizes the XFEMMaterialManager
 */
class XFEMTwoMaterialManagerConstraint : public ElemElemConstraint
{
public:
  XFEMTwoMaterialManagerConstraint(const InputParameters & parameters);

  virtual void computeResidual();
  virtual void computeJacobian();

  template <typename T>
  const MaterialProperty<T> * getMaterialProperty(const std::string & name) const;
  template <typename T>
  const MaterialProperty<T> * getMaterialPropertyOld(const std::string & name) const;
  template <typename T>
  const MaterialProperty<T> * getMaterialPropertyOlder(const std::string & name) const;

  template <typename T>
  const MaterialProperty<T> * getNeighborMaterialProperty(const std::string & name) const;
  template <typename T>
  const MaterialProperty<T> * getNeighborMaterialPropertyOld(const std::string & name) const;
  template <typename T>
  const MaterialProperty<T> * getNeighborMaterialPropertyOlder(const std::string & name) const;

protected:
  const XFEMElemPairMaterialManager & _manager;
  const XFEMElemPairMaterialManager & _manager_neighbor;
};

template <>
InputParameters validParams<XFEMTwoMaterialManagerConstraint>();

template <typename T>
const MaterialProperty<T> *
XFEMTwoMaterialManagerConstraint::getMaterialProperty(const std::string & name) const
{
  return &_manager.getMaterialProperty<T>(name);
}

template <typename T>
const MaterialProperty<T> *
XFEMTwoMaterialManagerConstraint::getMaterialPropertyOld(const std::string & name) const
{
  return &_manager.getMaterialPropertyOld<T>(name);
}

template <typename T>
const MaterialProperty<T> *
XFEMTwoMaterialManagerConstraint::getMaterialPropertyOlder(const std::string & name) const
{
  return &_manager.getMaterialPropertyOlder<T>(name);
}

template <typename T>
const MaterialProperty<T> *
XFEMTwoMaterialManagerConstraint::getNeighborMaterialProperty(const std::string & name) const
{
  return &_manager_neighbor.getMaterialProperty<T>("neighbor_" + name);
}

template <typename T>
const MaterialProperty<T> *
XFEMTwoMaterialManagerConstraint::getNeighborMaterialPropertyOld(const std::string & name) const
{
  return &_manager_neighbor.getMaterialPropertyOld<T>("neighbor_" + name);
}

template <typename T>
const MaterialProperty<T> *
XFEMTwoMaterialManagerConstraint::getNeighborMaterialPropertyOlder(const std::string & name) const
{
  return &_manager_neighbor.getMaterialPropertyOlder<T>("neighbor_" + name);
}

#endif // XFEMTWOMATERIALMANAGERCONSTRAINT_H
