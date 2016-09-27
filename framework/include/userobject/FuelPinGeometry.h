/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     BISON                     */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     Under Contract No. DE-AC07-05ID14517      */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/

#ifndef FUELPINGEOMETRY_H
#define FUELPINGEOMETRY_H

#include "GeneralUserObject.h"

class FuelPinGeometry : public GeneralUserObject
{
public:
  FuelPinGeometry(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  Real pellet_OD() const;
  Real pellet_ID() const;
  Real top_of_stack() const;
  Real bottom_of_stack() const;

  Real clad_OD() const;
  Real clad_ID() const;
  Real top_of_rod() const;
  Real bottom_of_rod() const;

private:
  MooseMesh & _mesh;
  Real _fuel_outer_radius;
  Real _fuel_inner_radius;
  Real _top_of_stack;
  Real _bottom_of_stack;
  Real _clad_outer_radius;
  Real _clad_inner_radius;
  Real _top_of_rod;
  Real _bottom_of_rod;
};

template<>
InputParameters validParams<FuelPinGeometry>();

#endif // FUELPINGEOMETRY_H
