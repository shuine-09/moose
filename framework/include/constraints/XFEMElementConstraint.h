/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef XFEMELEMENTCONSTRAINT_H
#define XFEMELEMENTCONSTRAINT_H

// MOOSE includes
#include "Constraint.h"
#include "NeighborCoupleableMooseVariableDependencyIntermediateInterface.h"
#include "MooseMesh.h"

// Forward Declarations
class XFEMElementConstraint;
class FEProblem;

template<>
InputParameters validParams<XFEMElementConstraint>();

class XFEMElementConstraint :
  public Constraint,
  public NeighborCoupleableMooseVariableDependencyIntermediateInterface
{
public:
  XFEMElementConstraint(const InputParameters & parameters);
  virtual ~XFEMElementConstraint();

  /**
   * Set the quadrature points and normal direction
   */
  virtual void setqRuleNormal(std::vector<Point> & quadrature_pts,
                              std::vector<Real> & quadrature_wts,
                              Point & normal);

  /**
   * Computes the residual for this element or the neighbor
   */
  virtual void computeElemNeighResidual(Moose::DGResidualType type);

  /**
   * Computes the residual for the current side.
   */
  virtual void computeResidual();

  /**
   * Computes the element/neighbor-element/neighbor Jacobian
   */
  virtual void computeElemNeighJacobian(Moose::DGJacobianType type);

  /**                                  
   * Computes the jacobian for the current side.
   */                         
  virtual void computeJacobian();

protected:
  FEProblem & _fe_problem;
  unsigned int _dim;

  const Elem * & _current_elem;

  /// The volume (or length) of the current element
  const Real & _current_elem_volume;

  /// The neighboring element
  const Elem * & _neighbor_elem;

  /// The volume (or length) of the current neighbor
  const Real & _neighbor_elem_volume;  

  /// Current side            
  unsigned int & _current_side;
  /// Current side element    
  const Elem * & _current_side_elem;

  /// The volume (or length) of the current side
  const Real & _current_side_volume;

  /// Coordinate system
  const Moose::CoordinateSystemType & _coord_sys;
  const MooseArray< Point > & _q_point;
  QBase * & _qrule;
  const MooseArray<Real> & _JxW;
  const MooseArray<Real> & _coord;

  Point _xfem_normal;
  std::vector<Point> _xfem_quad_pts;
  std::vector<Real> _xfem_quad_wts;

  unsigned int _i, _j;

  BoundaryID _boundary_id;

  /// Holds the current solution at the current quadrature point on the face.
  const VariableValue & _u;

  /// Holds the current solution gradient at the current quadrature point on the face.
  const VariableGradient & _grad_u;
  // shape functions
  const VariablePhiValue & _phi;
  const VariablePhiGradient & _grad_phi;
  // test functions

  /// Side shape function.
  const VariableTestValue & _test;
  /// Gradient of side shape function
  const VariableTestGradient & _grad_test;
  /// Normal vectors at the quadrature points
  const MooseArray<Point>& _normals;

  /// Side shape function.
  const VariablePhiValue & _phi_neighbor;
  /// Gradient of side shape function
  const VariablePhiGradient & _grad_phi_neighbor;

  /// Side test function
  const VariableTestValue & _test_neighbor;
  /// Gradient of side shape function
  const VariableTestGradient & _grad_test_neighbor;

  /// Holds the current solution at the current quadrature point
  const VariableValue & _u_neighbor;
  /// Holds the current solution gradient at the current quadrature point
  const VariableGradient & _grad_u_neighbor;

  /** 
   * This is the virtual that derived classes should override for computing the residual on neighboring element. 
   */ 
  virtual Real computeQpResidual(Moose::DGResidualType type) = 0;

  /**
   * This is the virtual that derived classes should override for computing the Jacobian on neighboring element.
   */
  virtual Real computeQpJacobian(Moose::DGJacobianType type) = 0;
};

#endif /* XFEMELEMENTCONSTRAINT_H_ */
