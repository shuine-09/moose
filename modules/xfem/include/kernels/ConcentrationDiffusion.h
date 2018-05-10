/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CONCENTRATIONDIFFUSION_H
#define CONCENTRATIONDIFFUSION_H

#include "Diffusion.h"
#include "Material.h"

// Forward Declarations
class ConcentrationDiffusion;

class ConcentrationDiffusion : public Diffusion
{
public:
  static InputParameters validParams();
  ConcentrationDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

private:
  /// String for the base name
  std::string _base_name;

  /// Material property name of diffusion coefficient
  std::string _diffusion_coefficient_name;

  /// Reference to diffusion coefficient property
  const MaterialProperty<Real> & _diffusion_coefficient;

  /// Pointer to material property that computes the derivatie of diffusion coefficient with repsect to temperature
  const MaterialProperty<Real> * const _diffusion_coefficient_dT;
};

#endif // CONCENTRATIONDIFFUSION_H
