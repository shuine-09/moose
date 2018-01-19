/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef LEVELSETCUTUSEROBJECT_H
#define LEVELSETCUTUSEROBJECT_H

#include "GeometricCutUserObject.h"

class LevelSetCutUserObject;

template <>
InputParameters validParams<LevelSetCutUserObject>();

class LevelSetCutUserObject : public GeometricCutUserObject
{
public:
  LevelSetCutUserObject(const InputParameters & parameters);

  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutEdge> & cut_edges,
                                    std::vector<Xfem::CutNode> & cut_nodes,
                                    Real time) const override;
  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutFace> & cut_faces,
                                    Real time) const override;

  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_edges,
                                     std::vector<Xfem::CutEdge> & cut_edges,
                                     Real time) const override;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_faces,
                                     std::vector<Xfem::CutFace> & cut_faces,
                                     Real time) const override;

  virtual const std::vector<Point>
  getCrackFrontPoints(unsigned int num_crack_front_points) const override;

protected:
  /// name of level set variable
  NonlinearVariableName _ls_var_name;

  /// level set variable
  MooseVariable & _ls_var;

  AuxiliarySystem & _aux_system;
  const NumericVector<Number> * _aux_solution;
  std::unique_ptr<NumericVector<Number>> _serialized_aux_solution;
};

#endif // LEVELSETCUTUSEROBJECT_H
