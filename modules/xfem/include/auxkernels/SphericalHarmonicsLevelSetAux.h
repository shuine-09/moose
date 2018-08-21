//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPHERICALHARMONICSLEVELSETAUX_H
#define SPHERICALHARMONICSLEVELSETAUX_H

#include "AuxKernel.h"

class SphericalHarmonicsLevelSetAux;

template <>
InputParameters validParams<SphericalHarmonicsLevelSetAux>();

class SphericalHarmonicsLevelSetAux : public AuxKernel
{
public:
  SphericalHarmonicsLevelSetAux(const InputParameters & parameters);

  virtual ~SphericalHarmonicsLevelSetAux() {}

  virtual void initialSetup() override;
  virtual Real computeValue() override;

protected:
  Real minimumDistanceForObject();

  Real computAngle(unsigned int object_id);

  Real ytotalfunction(int n, int m, Real x, Real y);

  Real computeSHdistance(unsigned int object_id);

  FileName _point_data_file_name;
  FileName _center_rotation_data_file_name;
  std::string _object_id_header;
  std::string _SH_m;
  std::string _SH_n;
  std::string _SH_real;
  std::string _SH_imaginary;
  std::string _x_center;
  std::string _y_center;
  std::string _x_rotation;
  std::string _y_rotation;

  std::map<int, std::vector<std::array<Real, 4>>> _point_data;
  std::map<int, std::array<Real, 4>> _object_centers_rotation;

private:
  inline Real factorial(int n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n; }
};
#endif // SPHERICALHARMONICSLEVELSETAUX_H
