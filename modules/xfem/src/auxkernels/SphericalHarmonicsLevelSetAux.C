//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SphericalHarmonicsLevelSetAux.h"
#include "Conversion.h"
#include "DelimitedFileReader.h"
#include <boost/math/special_functions/legendre.hpp>

registerMooseObject("XFEMApp", SphericalHarmonicsLevelSetAux);

template <>
InputParameters validParams<SphericalHarmonicsLevelSetAux>() // Change to SH distance
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<FileName>(
      "center_rotation_data_file",
      "Name of the CSV file containing the object center coordinate and "
      "rotation data");
  params.addRequiredParam<FileName>(
      "point_data_file", "Name of the CSV file containing the  spherical harmonics coefficients.");

  // The above file has 5 columns: object_id, x_center, y_center, x_rotation,
  // y_rotation
  // The above file has 5 columns: object_id, SH_n, SH_m, SH_real, SH_imaginary
  params.addRequiredParam<std::string>(
      "object_id_header", "Header name for the column in the CSV files with the object IDs");
  params.addRequiredParam<std::string>(
      "x_center", "Header name for the column in the CSV file with the x coordinate");
  params.addRequiredParam<std::string>(
      "y_center", "Header name for the column in the CSV file with the y coordinate");
  params.addRequiredParam<std::string>(
      "x_rotation", "Header name for the column in the CSV file with the x rotation degree");
  params.addRequiredParam<std::string>(
      "y_rotation", "Header name for the column in the CSV file with the y rotation degree");
  params.addRequiredParam<std::string>(
      "SH_m", "Header name for the column in the CSV file with the m value");
  params.addRequiredParam<std::string>(
      "SH_n", "Header name for the column in the CSV file with the n value");
  params.addRequiredParam<std::string>(
      "SH_real", "Header name for the column in the CSV file with the real coefficient");
  params.addRequiredParam<std::string>("SH_imaginary",
                                       "Header name for the column in the CSV "
                                       "file with the imaginary coefficient");
  // params.addRequiredParam<std::string>("z_rotation", "Header name for the
  // column in the CSV file with the z rotation degree");
  return params;
}

SphericalHarmonicsLevelSetAux::SphericalHarmonicsLevelSetAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _point_data_file_name(getParam<FileName>("point_data_file")),
    _center_rotation_data_file_name(getParam<FileName>("center_rotation_data_file")),
    _object_id_header(getParam<std::string>("object_id_header")),
    _SH_m(getParam<std::string>("SH_m")),
    _SH_n(getParam<std::string>("SH_n")),
    _SH_real(getParam<std::string>("SH_real")),
    _SH_imaginary(getParam<std::string>("SH_imaginary")),
    _x_center(getParam<std::string>("x_center")),
    _y_center(getParam<std::string>("y_center")),
    _x_rotation(getParam<std::string>("x_rotation")),
    _y_rotation(getParam<std::string>("y_rotation"))
{
  if (!isNodal())
    mooseError("SphericalHarmonicsLevelSetAux can only be run on a nodal variable");
}

void
SphericalHarmonicsLevelSetAux::initialSetup()
{
  // Read and store the spherical harmonics data for each object
  MooseUtils::DelimitedFileReader pd_csv_reader(_point_data_file_name);
  pd_csv_reader.setDelimiter(",");
  pd_csv_reader.setHeaderFlag(MooseUtils::DelimitedFileReader::HeaderFlag::ON);
  pd_csv_reader.read();

  std::vector<Real> object_id_real = pd_csv_reader.getData(_object_id_header);
  std::vector<int> object_id(object_id_real.size());
  for (unsigned int i = 0; i < object_id_real.size(); ++i)
    object_id[i] = static_cast<int>(object_id_real[i]);

  std::vector<Real> m = pd_csv_reader.getData(_SH_m);
  std::vector<Real> n = pd_csv_reader.getData(_SH_n);
  std::vector<Real> SHreal = pd_csv_reader.getData(_SH_real);
  std::vector<Real> SHimag = pd_csv_reader.getData(_SH_imaginary);

  if ((object_id.size() != m.size()) || (object_id.size() != n.size()))
    mooseError("Columns of data in point data csv file must all be the same size");

  for (unsigned int i = 0; i < m.size(); ++i)
  {
    _point_data[object_id[i]].push_back(std::array<Real, 4>{{n[i], m[i], SHreal[i], SHimag[i]}});
  }

  // for (unsigned int i = 0; i < m.size(); ++i)
  // {
  //   for (unsigned int j = 0; j < _point_data[object_id[i]].size(); ++j)
  //   {
  //     std::cout << "data[" << object_id[i] << "] = " << (_point_data[object_id[i]])[j][0] << ", "
  //               << (_point_data[object_id[i]])[j][1] << ", " << (_point_data[object_id[i]])[j][2]
  //               << ", " << (_point_data[object_id[i]])[j][3] << std::endl;
  //   }
  // }

  // Read and store the center coordinate data for each object
  MooseUtils::DelimitedFileReader oc_csv_reader(_center_rotation_data_file_name);
  oc_csv_reader.setDelimiter(",");
  oc_csv_reader.setHeaderFlag(MooseUtils::DelimitedFileReader::HeaderFlag::ON);
  oc_csv_reader.read();
  object_id_real = oc_csv_reader.getData(_object_id_header);
  object_id.resize(object_id_real.size());

  // if (_point_data.size() != object_id.size())
  //   mooseError("Number of objects for point data and center data must match");

  for (unsigned int i = 0; i < object_id_real.size(); ++i)
  {
    object_id[i] = static_cast<int>(object_id_real[i]);

    if (_point_data.find(object_id[i]) == _point_data.end())
      mooseError("No point data for object " + Moose::stringify(object_id[i]) +
                 " in center data file");
  }
  std::vector<Real> x = oc_csv_reader.getData(_x_center);
  std::vector<Real> y = oc_csv_reader.getData(_y_center);
  std::vector<Real> x_rotation = oc_csv_reader.getData(_x_rotation);
  std::vector<Real> y_rotation = oc_csv_reader.getData(_y_rotation);

  if ((object_id.size() != x.size()) || (object_id.size() != y.size()))
    mooseError("Columns of data in object center csv file must all be the "
               "same size");

  for (unsigned int i = 0; i < object_id.size(); ++i)
    _object_centers_rotation[object_id[i]] = {{x[i], y[i], x_rotation[i], y_rotation[i]}};

  // for (unsigned int i = 0; i < object_id.size(); ++i)
  // {
  //   std::cout << "data[" << object_id[i] << "] = " << (_object_centers_rotation[object_id[i]])[0]
  //             << ", " << (_object_centers_rotation[object_id[i]])[1] << ", "
  //             << (_object_centers_rotation[object_id[i]])[2] << ", "
  //             << (_object_centers_rotation[object_id[i]])[3] << std::endl;
  // }
}

Real
SphericalHarmonicsLevelSetAux::computeValue()
{
  return minimumDistanceForObject();
}

Real
SphericalHarmonicsLevelSetAux::computAngle(unsigned int object_id)
{
  Point center(_object_centers_rotation[object_id][0], _object_centers_rotation[object_id][1], 0);
  Point center_to_node = *_current_node - center;

  Real atan2_v = atan2(center_to_node(1), center_to_node(0));
  return atan2_v > 0 ? atan2_v : 2 * libMesh::pi + atan2_v;
}

Real
SphericalHarmonicsLevelSetAux::ytotalfunction_real(int n, int m, Real x, Real y)
{
  // Calculations all from Garboczi 2002 spherical harmonics paper i = sqrt(-1);
  Real factor = sqrt(((2 * n + 1) * factorial(n - m)) / (4 * libMesh::pi * factorial(n + m)));
  Real yout = boost::math::legendre_p(n, m, x);
  Real exp_value_real = cos(m * y);
  return factor * yout * exp_value_real;
}

Real
SphericalHarmonicsLevelSetAux::ytotalfunction_imaginary(int n, int m, Real x, Real y)
{
  Real factor = sqrt(((2 * n + 1) * factorial(n - m)) / (4 * libMesh::pi * factorial(n + m)));
  Real yout = boost::math::legendre_p(n, m, x);
  Real exp_value_imaginary = sin(m * y);
  return factor * yout * exp_value_imaginary;
}

Real
SphericalHarmonicsLevelSetAux::computeSHdistance(unsigned int object_id)
{
  Real radius = 0.0;

  std::vector<std::array<Real, 4>> coef = _point_data[object_id];

  Real y = computAngle(object_id);

  for (auto & c : coef)
  {
    Real n = c[0];
    Real m = c[1];
    Real real = c[2];
    Real imaginary = c[3];
    Real y_function_value_real = ytotalfunction_real(n, m, 0, y);
    Real y_function_value_imaginary = ytotalfunction_imaginary(n, m, 0, y);
    radius = radius + (real * y_function_value_real - imaginary * y_function_value_imaginary);
  }

  return radius;
}

Real
SphericalHarmonicsLevelSetAux::minimumDistanceForObject()
{
  Real min_dist = std::numeric_limits<Real>::max();

  for (std::map<int, std::array<Real, 4>>::iterator iter = _object_centers_rotation.begin();
       iter != _object_centers_rotation.end();
       ++iter)
  {
    int object_id = iter->first;
    Point center((iter->second)[0], (iter->second)[1], 0);
    Real dist_from_node_to_center = (*_current_node - center).norm();
    Real radius = computeSHdistance(object_id);

    Real dist = dist_from_node_to_center - radius;

    if (dist < min_dist)
      min_dist = dist;
  }

  return min_dist;
}
