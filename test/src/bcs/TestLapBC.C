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

#include "TestLapBC.h"
#include "Function.h"

#include "libmesh/numeric_vector.h"
#include "libmesh/string_to_enum.h"

#include <cmath>
#include "Assembly.h"
#include "NonlinearSystem.h"

template<>
InputParameters validParams<TestLapBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  return params;
}

TestLapBC::TestLapBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _second_u_from_interface(second()),
    _second_phi(_assembly.secondPhiFace()),
    _second_test(_var.secondPhiFace()),
    _second_u(_is_implicit ?_var.secondSln() : _var.secondSlnOld()),
    //_nodal_soln(_var.nodalSln())
    // value() should call nodalSln() on the _var.  This gives me a
    // vector with 2 entries on the face, but shouldn't it be 3?
    // Also, note that we need the nodal solution from the full
    // element in order to compute second derivatives on the face,
    // since the 2nd derivative depends on all the basis functions.
    // Perhaps it is actually returning _var.sln(), since this variable
    // is not actually nodal?  This is the same thing as _u in that case.
    //_nodal_soln(value())
    _second_test_from_interface(secondTest())
{
}

Real
TestLapBC::computeQpResidual()
{
  Real r = 0;
  // Moose::out << "_q_point[" << _qp << "]=" << _q_point[_qp] << std::endl;

  // In this test, _second_u is basically always zero, as the true solution is u==const.
  // WJ: I cannot compile it as norm() is not defined in my libmesh? 
  //if (_second_u[_qp].norm() > 1.e-6)
  //  Moose::out << "_second_u[" << _qp << "]=" << _second_u[_qp] << std::endl;

  // Print the diff between the result of getting _second_u in two
  // different ways if it is large enough.  I don't currently see any
  // printouts from these, so this seems to be working...
  
  // WJ: I cannot compile it as norm() is not defined in my libmesh? 
  //Real diff_second = (_second_u_from_interface[_qp] - _second_u[_qp]).norm();
  //if (diff_second > 1.e-6)
  //  Moose::out << "diff_second[" << _qp << "]=" << diff_second << std::endl;

  // Second derivatives of the test function.  These values match what
  // is being computed during Assembly::reinitFEFace(), unlike the values you get
  // by initializing a variable by calling secondTest().
  // Moose::out << "_second_test[" << _i << "][" << _qp << "]=" << _second_test[_i][_qp] << std::endl;

  // Are "_second_test" and "_second_phi" the same thing? (One comes
  // from _assembly, the other comes from _var...)  This line of code
  // actually segfaults, so possibly _second_phi is not even ready to
  // be accessed in computeQpResidual().
  // Moose::out << "phi_xx - test_xx = " << (_second_test[_i][_qp] - _second_phi[_i][_qp]).norm() << std::endl;

  // // Compute sum_i (phi_xx_i + phi_yy_i + phi_zz_i).  This should sum
  // // to zero at every quadrature point.
  // Real lap_phi = 0.;
  // for (unsigned int i=0; i<_second_test.size(); ++i)
  //   lap_phi += _second_test[i][_qp].tr();
  // if (std::abs(lap_phi) > 1.e-12)
  //   mooseError("lap_phi=" << lap_phi << ", should be zero.");

  // For some reason, the size of _nodal_soln is 1 if we call
  // _var.nodalSln() to get it?  It seems like there should be 3
  // values on a QUAD9.  Ah, the reason is that these are the
  // *quadrature* point values (there are two qp's if we use
  // SECOND-order Gauss quadrature, and that's what this is returning
  // us.)

  // How to compute sum_i u_i (phi_xx_i + phi_yy_i + phi_zz_i)?  Note
  // that _u refers to solution values at qp's, not nodal values.

  // // These two things do the same thing -- give back a vector of values at qps.
  // for (unsigned int i=0; i<_nodal_soln.size(); ++i)
  //   Moose::out << "_nodal_soln[" << i << "]=" << _nodal_soln[i] << std::endl;
  //
  // for (unsigned int i=0; i<_u.size(); ++i)
  //   Moose::out << "_u[" << i << "]=" << _u[i] << std::endl;

  // What kind of element are we on?  Even though we're in a BC, this returns QUAD9.
  // Moose::out << "On a " << Utility::enum_to_string(_current_elem->type()) << " element." << std::endl;

  // IntegratedBC provides _current_side_elem for the side!
  // Moose::out << "The side is a " << Utility::enum_to_string(_current_side_elem->type()) << " element." << std::endl;

  // // We can access the DofMap through the var and ask it for the dof indices for this side's nodes.
  // std::vector<dof_id_type> di;
  // _var.dofMap().dof_indices(_current_side_elem, di);
  // for (unsigned int i=0; i<di.size(); ++i)
  //   Moose::out << di[i] << " ";
  // Moose::out << std::endl;
  //
  // // From this, we can access solution values from the System.
  // NumericVector<Number> & soln = _fe_problem.getNonlinearSystem().solution();
  // for (unsigned int i=0; i<di.size(); ++i)
  //   Moose::out << soln(i) << " ";
  // Moose::out << std::endl;

  // // We can access the DofMap through the var and ask it for the dof indices for this entire element.
  // std::vector<dof_id_type> di;
  // _var.dofMap().dof_indices(_current_elem, di);
  // for (unsigned int i=0; i<di.size(); ++i)
  //   Moose::out << di[i] << " ";
  // Moose::out << std::endl;
  //
  // // From this, we can access solution values from the System.  It
  // // seems like we're seeing finite differencing (triggered by
  // // -snes_check_jacobian) when we do this.  I get printouts like:
  // //
  // // The side is a EDGE3 element.
  // // 0 1 2 3 4 5 6 7 8
  // // 0 0 0 0 0 1.49012e-08 0 0 0
  // const NumericVector<Number> * & current_soln = _fe_problem.getNonlinearSystem().currentSolution();
  // for (unsigned int i=0; i<di.size(); ++i)
  //   Moose::out << (*current_soln)(i) << " ";
  // Moose::out << std::endl;

  // For this simple example, the second derivative is identically 0,
  // so we should be able to get roughly the same result (including
  // incorrect Jacobians) just by returning 0.  We can simulate this
  // by simply commenting out the next line:
  // return 0: 0.51484 = ||J - Jfd||/||J||
  //           3.82971 = ||J - Jfd||
  // return _second_u.tr(): 0.580818 = ||J - Jfd||/||J||
  //                        4.32049  = ||J - Jfd||
  // They are slightly different, but probably not different enough
  // for the _second_u values to make a difference?

  // Verify that the problem is not with _test[_i]
  // r += 1.0 * _test[_i][_qp];

  // What we actually want to return for the residual.
  //r += _second_u[_qp].tr() * _test[_i][_qp];
    r += 10 * (_second_u[_qp].tr() - 2) * _test[_i][_qp]; // 10 is the penalty parameter which penalize the \Lap U = 2 on the boundary.
  return r;
}

Real
TestLapBC::computeQpJacobian()
{
  Real r = 0;
  // Moose::out << "_q_point[" << _qp << "]=" << _q_point[_qp] << std::endl;
  // Moose::out << "_second_phi[" << _j << "][" << _qp << "]=" << _second_phi[_j][_qp] << std::endl;

  // Compare second derivative values obtained in different ways.  It
  // appears to be safe to access _second_phi from the
  // computeQpJacobian() routine...  These two arrays appear to be identical.
  {
    //Real normdiff = (_second_test[_i][_qp] - _second_phi[_i][_qp]).norm();
    //if (normdiff > 1.e-12)
    //  Moose::out << "phi_xx - test_xx = " << normdiff << std::endl;
  }

  // Compare test function second derivatives obtained from _var with
  // those obtained through the interface. -- OK, these values are
  // definitely different, and this is probably a bug!  But it's not the only bug, because
  // we aren't even using _second_test_from_interface and there are still issues
  // with the Jacobian.
//  {
//    Real normdiff = (_second_test[_i][_qp] - _second_test_from_interface[_i][_qp]).norm();
//    if (normdiff > 1.e-12)
//      Moose::out << "test_xx - test_interface_xx = " << normdiff << std::endl;
//  }

  // What if we used _second_test again.  Does that make the Jacobian
  // tester work?  No, it does not seem to make any difference...
  // r += _second_test[_j][_qp].tr() * _test[_i][_qp];

  // What if we used _second_test_from_interface instead?  This actually gives a different
  // result that is still not correct, but at least it's different, which is consistent
  // with using different shape function values.
  // r += _second_test_from_interface[_j][_qp].tr() * _test[_i][_qp];

  // Verify that the problem is not with _test[_i] (this works fine).
  // r += 0.;

  // Since we believe this example to have a zero second derivative,
  // it should also "work" to simply return 0 for the Jacobian.
  // Unfortunately, it does not work, so irrespective of what
  // _second_phi actually is, it appears that _second_u may be wrong?
  // r += 0.;

  // What I think we want to use:
  //r += _second_phi[_j][_qp].tr() * _test[_i][_qp];
  r += 10 * _second_phi[_j][_qp].tr() * _test[_i][_qp]; // 10 is the penalty parameter which penalize the \Lap U = 2 on the boundary.
  return r;
}
