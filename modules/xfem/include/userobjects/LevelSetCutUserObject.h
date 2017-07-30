/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef LEVELSETCUTUSEROBJECT_H
#define LEVELSETCUTUSEROBJECT_H

#include "GeometricCut2DUserObject.h"
#include "AuxiliarySystem.h"

// Forward declarations
class LevelSetCutUserObject;
class AuxiliarySystem;

template <>
InputParameters validParams<LevelSetCutUserObject>();

class LevelSetCutUserObject : public GeometricCut2DUserObject
{
public:
  LevelSetCutUserObject(const InputParameters & parameters);

  virtual void initialize() override{};
  virtual void execute() override{};
  virtual void finalize() override{};
  virtual bool active(Real time) const override { return true; };
  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<CutEdge> & cut_edges,
                                    Real time) const override;

protected:
  /// name of level set variable
  AuxVariableName _ls_var_name;

  AuxiliarySystem & _aux_system;
  /// level set variable
  MooseVariable & _ls_var;
  const NumericVector<Number> * _aux_solution;
};

#endif // LEVELSETCUTUSEROBJECT_H
