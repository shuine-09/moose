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

#ifndef THERMALMAMOX_H
#define THERMALMAMOX_H

#include "Material.h"
#include "Function.h"


//Forward Declarations
class ThermalMAMOX;

template<>
InputParameters validParams<ThermalMAMOX>();

/**
 * Temperature and burnup dependent thermal properties of Uranium Dioxide
 */
class ThermalMAMOX : public Material
{
  public:
    ThermalMAMOX(const InputParameters & parameters);

  protected:
    virtual void computeQpProperties();


  private:
    bool _has_temp;
    const VariableValue  & _temp;
    bool _has_porosity;
    const VariableValue & _porosity;
    Function * _fp_function;


    //bool _has_porosity;
    //VariableValue  & _porosity;

    //  MaterialProperty<Real> & _density;
    MaterialProperty<Real> & _thermal_conductivity;
    MaterialProperty<Real> & _thermal_conductivity_dT;
    MaterialProperty<Real> & _specific_heat;

    //  Real _initial_porosity;
    //  Real _initial_density;
    Real _stoech_dev;
    Real _Am_content;
    Real _Np_content;

};

#endif //THERMALMAMOX_H
