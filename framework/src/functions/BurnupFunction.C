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

#include "BurnupFunction.h"
#include "FuelPinGeometry.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <functional>

namespace {
const Real invalid_geometry = -99999999;
}

template<>
InputParameters validParams<BurnupFunction>()
{
  InputParameters params = validParams<Function>();

  std::vector<Real> i_enrich_default(BurnupFunction::NUM_HM_ISOTOPES);
  i_enrich_default[0] = 0.05;
  i_enrich_default[1] = 0.95;
  i_enrich_default[2] = 0.0;
  i_enrich_default[3] = 0.0;
  i_enrich_default[4] = 0.0;
  i_enrich_default[5] = 0.0;

  params.addParam<Real>("value", 1.0, "Default/scaling value.");
  params.addParam<bool>("rpf_active", true, "Flag for turning radial power factor on.");
  params.addParam<bool>("include_gadolinia", false, "Set to true for fuel containing burnable absorber Gd2O3.");
  params.addRangeCheckedParam<unsigned>("num_radial", 80, "num_radial>1", "Number of radial grid points.");
  params.addRangeCheckedParam<unsigned>("num_axial", 20, "num_axial>1", "Number of axial grid points.");
  params.addParam<unsigned>("axial_axis",1,"Coordinate axis of the axial direction of the fuel stack (0, 1, or 2 for x, y, or z");
  params.addParam<Real>("a_lower", "The lower axial coordinate.");
  params.addParam<Real>("a_upper", "The upper axial coordinate.");
  params.addParam<Real>("fuel_inner_radius", 0, "The inner radius.");
  params.addParam<Real>("fuel_outer_radius",  0.0041, "The outer radius.");
  params.addRangeCheckedParam<Real>("bias", 1.0, "bias>=0.5 & bias<=2.0", "Bias for radial point spacing.  Must be between 0.5 and 2.0");
  params.addRangeCheckedParam<Real>("fuel_volume_ratio", 1.0, "fuel_volume_ratio>0 & fuel_volume_ratio<=1", "Reduction factor for deviation from right circular cylinder fuel.  The ratio of actual volume to right circular cylinder volume.");
  params.addRequiredRangeCheckedParam<Real>("density", "density>0", "The density.");
  // 205 MeV/fission * 1.6022e-13 J/MeV = 3.28451e-11 J/fission
  params.addRangeCheckedParam<Real>("energy_per_fission", 3.28451e-11, "energy_per_fission>0", "Energy released per fission (J/fission)");
  params.addParam<std::vector<Real> >("i_enrich", i_enrich_default, "The initial enrichment for the considered isotopes.");
  params.addParam<Real>("wtfract_gadolinia", 0., "Initial weight fraction of Gd2O3");
  params.addParam<std::vector<Real> >("sigma_c", "The capture cross sections for the considered isotopes.");
  params.addParam<std::vector<Real> >("sigma_f", "The fission cross sections for the considered isotopes.");
  params.addParam<std::vector<Real> >("sigma_a_thermal", "The thermal absorption cross sections for the considered isotopes.");
  params.addParam<FunctionName>("rod_ave_lin_pow", "", "The power function");
  params.addParam<FunctionName>("axial_power_profile", "", "The axial power profile function");
  params.addParam<FunctionName>("rpf_input", "", "The radial power profile function. Used to specify the rpf from input");
  params.addParam<Real>("p1", "Distribution function coefficient p1.  If not given, will take default value based on reactor_type");
  params.addParam<UserObjectName>("fuel_pin_geometry", "", "Name of the UserObject that reads the pin geometry from the mesh." );
  params.addRangeCheckedParam<Real>("initial_burnup", 0, "initial_burnup>=0", "Initial burnup to be applied in units of MWd/kgU");
  params.addPrivateParam<Real>("heavy_metal_molar_mass");
  params.addPrivateParam<Real>("fuel_molar_mass");

  MooseEnum reactor_type("LWR HWR", "LWR");
  MooseEnum fuel_type("UO2 U3Si2", "UO2");
  params.addParam<MooseEnum>("reactor_type", reactor_type, "Reactor type.  Choices are " + reactor_type.getRawNames());
  params.addParam<MooseEnum>("fuel_type", fuel_type, "Fuel type. Choices are " + fuel_type.getRawNames());
  // Hide BurnupFunction from the input file dump
  params.addPrivateParam<std::string>("built_by_action", "");
  return params;
}

/* The BISON burnup model is based on:
   K. Lassmann et al., The radial distribution of plutonium in high burnup UO2 fuels, JNM 208 (1994) 223-231.
   K. Lassmann et al., Extension of the TRANSURANUS burnup model to heavy water reactor conditions, JNM 255 (1998) 222-233
   Y. Rashid et al., Fuel analysis and licensing code: FALCON MOD01, Tech Rep EPRI 1011308, 2004.
*/

BurnupFunction::BurnupFunction(const InputParameters & params) :
  Function(params),
  FunctionInterface(this),
  _ralp(getParam<FunctionName>("rod_ave_lin_pow") != "" ?
        &getFunction("rod_ave_lin_pow") :
        NULL),
  _axial_profile(getParam<FunctionName>("axial_power_profile") != "" ?
                 &getFunction("axial_power_profile") :
                 NULL),
  _rpf_input(getParam<FunctionName>("rpf_input") != "" ?
                 &getFunction("rpf_input") :
                 NULL),
  _pin_geometry(getParam<UserObjectName>("fuel_pin_geometry")),
  _value(getParam<Real>( "value" )),
  _rpf_active(getParam<bool>( "rpf_active" )),
  _include_gadolinia(getParam<bool>( "include_gadolinia" )),
  _num_isotopes(_include_gadolinia?NUM_HM_ISOTOPES+NUM_GD_ISOTOPES:NUM_HM_ISOTOPES),
  _num_radial(getParam<unsigned>("num_radial")),
  _num_axial(getParam<unsigned>("num_axial")),
  _axial_axis(getParam<unsigned>("axial_axis")),
  _a_lower(isParamValid("a_lower") ? getParam<Real>("a_lower") : invalid_geometry),
  _a_upper(isParamValid("a_upper") ? getParam<Real>("a_upper") : invalid_geometry),
  _height(invalid_geometry),
  _r_inner(getParam<Real>("fuel_inner_radius")),
  _r_outer(getParam<Real>("fuel_outer_radius")),
  _bias(getParam<Real>("bias")),
  _fuel_volume_ratio(getParam<Real>("fuel_volume_ratio")),
  _time_last_requested(declareRestartableData<Real>("time_last_requested", -std::numeric_limits<Real>::max())),
  _density( getParam<Real>("density") ),
  _initial_burnup( getParam<Real>("initial_burnup") ),
  _initial_burnup_complete(declareRestartableData<bool>("initial_burnup_complete", false)),
  _porosity_fraction(0.),
  _i_enrich(getParam<std::vector<Real> >("i_enrich")),
  _wtfract_gadolinia(getParam<Real>("wtfract_gadolinia")),
  // MeV/fission * 1.6022e-13 J/MeV = J/fission
  _energy_per_fission( getParam<Real>("energy_per_fission" ) ), // J/fission
  _r(_num_radial),
  _h(_num_axial),
  _volume_frac(_num_radial),
  _p1(p1()),
  _sigma_f(sigma_f()),
  _sigma_c(sigma_c()),
  _sigma_a(sigma_a()),
  _sigma_a_thermal(sigma_a_thermal()),
  // Concentrations are atoms/m^3 of UO2
  // Index [i][j][k] -> [isotope][height][radius]
  _concentration(
    declareRestartableData<std::vector<std::vector<std::vector<Real> > > >("concentration",
      std::vector<std::vector<std::vector<Real> > >(_num_isotopes,
        std::vector<std::vector<Real> >(_num_axial, std::vector<Real>(_num_radial, 0.0))))),
  _concentration_old(_num_isotopes, std::vector<std::vector<Real> >(_num_axial, std::vector<Real>(_num_radial, 0.0))),
  // Index [i][j] -> [isotope][height]
  _average_concentration(
    declareRestartableData<std::vector<std::vector<Real> > >("average_concentration",
      std::vector<std::vector<Real> >(_num_isotopes, std::vector<Real>(_num_axial, 0.0)))),
  _average_concentration_old(
    declareRestartableData<std::vector<std::vector<Real> > >("average_concentration_old",
      std::vector<std::vector<Real> >(_num_isotopes, std::vector<Real>(_num_axial, 0.0)))),
  // Index [j][k] -> [height][radius]
  _burnup(
    declareRestartableData<std::vector<std::vector<Real> > >("burnup",
      std::vector<std::vector<Real> >(_num_axial, std::vector<Real>(_num_radial, 0.0)))), // MWd/kgU
  _burnup_old(_num_axial, std::vector<Real>(_num_radial, 0.0)),
  _average_burnup( declareRestartableData<std::vector<Real> >("average_burnup", std::vector<Real>(_num_axial, 0.0)) ),
  _fission_rate(
    declareRestartableData<std::vector<std::vector<Real> > >("fission_rate",
      std::vector<std::vector<Real> >(_num_axial, std::vector<Real>(_num_radial, 0.0)))),
  _fission_rate_old(
    declareRestartableData<std::vector<std::vector<Real> > >("fission_rate_old",
      std::vector<std::vector<Real> >(_num_axial, std::vector<Real>(_num_radial, 0.0)))),
  _rpf(
    declareRestartableData<std::vector<std::vector<Real> > >("RPF",
      std::vector<std::vector<Real> >(_num_axial, std::vector<Real>(_num_radial, 1.0)))),
  // Index [k] -> [radius]
  _f_normalized(_num_radial, 0.0),
  // 205 MeV | 1.6022e-13 J | W*s |   day   |   MW
  // ------------------------------------------------ = 3.8e-22 MWd/fission
  // fission |     MeV      |  J  | 86400 s |  1e6 W
  _MWdperFission( _energy_per_fission / 86400e6 ),
  _is_initialized(false)
{
  if (getParam<MooseEnum>("fuel_type") == "U3Si2")
  {
    _hm_molar = 0.714;
    _fuel_molar = 0.770;
  }
  else if (getParam<MooseEnum>("fuel_type") == "UO2") // UO2
  {
    _hm_molar = 0.238;
    _fuel_molar = 0.270;
  }
  if (isParamValid("heavy_metal_molar_mass"))
  {
    _hm_molar = getParam<Real>("heavy_metal_molar_mass");
  }
  if (isParamValid("fuel_molar_mass"))
  {
    _fuel_molar = getParam<Real>("fuel_molar_mass");
  }

  _fisskgU_per_MWdaU = 86400e6 * (_hm_molar) / (_energy_per_fission * 6.0221415e23);

  // Number of atoms per m^3 of fuel
  // _fuel_molar = kg per mole for fuel compound
  // 6.0221415e23 = Avagadro's number
  _ntot_hm = (1.-_porosity_fraction ) * _density / (_fuel_molar) * 6.0221415e23 * (1.-_wtfract_gadolinia);
  // a_lower and a_upper must be specified as there is no default value, but r_inner and r_outer
  // have default values.
  if ( _pin_geometry.empty() )
  {
    if ( _a_lower == invalid_geometry )
      mooseError( "Burnup: If fuel_pin_geometry is not specified, then a_lower must be specified." );
    if ( _a_upper == invalid_geometry )
      mooseError( "Burnup: If fuel_pin_geometry is not specified, then a_upper must be specified." );
  }

  if (_bias < 0.5 || _bias > 2.0)
  {
    std::stringstream msg;
    msg << "Bias was given as " << _bias << ". "
        << "Allowed range: 0.5 <= bias <= 2.0";
    mooseError(msg.str());
  }
  if (_fuel_volume_ratio > 1)
  {
    mooseError("Burnup: fuel_volume_ratio must be <= 1");
  }

  if(_include_gadolinia && _wtfract_gadolinia == 0.)
  {
    std::stringstream msg;
    msg << "With include_gadolinia = true, please specify a value > 0 for wtfract_gadolinia.";
    mooseError(msg.str());
  }

  if(!_include_gadolinia && _wtfract_gadolinia > 0.)
  {
    std::stringstream msg;
    msg << "Need to set include_gadolinia = true if _wtfract_gadolinia > 0.";
    mooseError(msg.str());
  }

  if(getParam<MooseEnum>("fuel_type") == "U3Si2" && _include_gadolinia)
  {
    std::stringstream msg;
    msg << "Cannot include gadolinia when modeling U3Si2 fuel.";
    mooseError(msg.str());
  }
  if((getParam<MooseEnum>("fuel_type") == "U3Si2") && (getParam<MooseEnum>("reactor_type") == "HWR"))
  {
    std::stringstream msg;
    msg << "There are no heavy water reactor cross sections for U3Si2. Use reactor_type = LWR";
    mooseError(msg.str());
  }

  if (_sigma_f.size() != _sigma_c.size() ||
      _sigma_c.size() != _sigma_a_thermal.size() ||
      _sigma_a_thermal.size() != _num_isotopes)
  {
    std::stringstream msg;
    msg << "All cross section lists must have "
        << _num_isotopes
        << " entries.";
    mooseError(msg.str());
  }

  mooseAssert(_i_enrich.size() == NUM_HM_ISOTOPES, "Initial enrichment vector has wrong size");

  if (!_rpf_input)
  {
    std::vector<Real> i_concen(_i_enrich);
    std::transform(i_concen.begin(), i_concen.end(), i_concen.begin(),
                   std::bind1st(std::multiplies<Real>(), _ntot_hm));
    if(_include_gadolinia)
    {
      const double ntot_gd = ( 1.-_porosity_fraction ) * _density / 0.3625 * 2. * 6.0221415e23 * _wtfract_gadolinia;
      const double isotopic_fract_155gd = 0.1480;
      const double isotopic_fract_157gd = 0.1565;
      i_concen.push_back(ntot_gd * isotopic_fract_155gd);
      i_concen.push_back(ntot_gd * isotopic_fract_157gd);
    }

    for (unsigned int i(0); i < _num_isotopes; ++i)
    {
      for (unsigned int j(0); j < _num_axial; ++j)
      {
        _average_concentration[i][j] = _average_concentration_old[i][j] = i_concen[i];
        for (unsigned int k(0); k < _num_radial; ++k)
          _concentration[i][j][k] = _concentration_old[i][j][k] = i_concen[i];
      }
    }
  }
  // The rest of the setup has to be done after the geometry is known, which is after the
  // FuelPinGeometry UserObject has been constructed.
}

void
BurnupFunction::initialize( Real time )
{
  if ( ! _is_initialized )
  {
    if ( ! _pin_geometry.empty() )
    {
      const FuelPinGeometry & geometry = getUserObject<FuelPinGeometry>("fuel_pin_geometry");
      _a_lower = geometry.bottom_of_stack();
      _a_upper = geometry.top_of_stack();
      _r_inner = geometry.pellet_ID() / 2;
      _r_outer = geometry.pellet_OD() / 2;
    }

    _height = _a_upper-_a_lower;

    const Real h( _height / (_num_axial-1) );
    _h[0] = _a_lower;
    for (unsigned j(1); j < _num_axial; ++j)
    {
      _h[j] = _h[j-1] + h;
    }

    computeRadialLocations( _r_inner, _r_outer, _bias, _r );

    computeVolumeFractions( h, _num_axial, _r, _volume_frac );

    // Conversion from power to fission rate
    // p_to_fr => fission/J/m^2
    _p_to_fr = 1.0 / ( _energy_per_fission * _fuel_volume_ratio * M_PI * (_r[_num_radial-1]*_r[_num_radial-1] - _r[0]*_r[0]) );

    bool initialized_burnup = initializeBurnup();
    if ( initialized_burnup )
    {
      _time_last_requested = time;
    }

    if ( ! _rpf_input )
    {
      // p1 is 3.45 for LWR, 2.21 for Halden HWR (Handbook of Nuclear Engr, Vol. 2, also FRAPCON manual)
      // p2 is 3    for LWR and HWR
      // p3 is 0.45 for LWR and HWR
      const Real p2(3);
      const Real p3(0.45);
      Real average_f(0);
      for (unsigned int k(0); k < _num_radial; ++k)
      {
        // The equation for f requires length in millimeters
        // Multiply by 1000 to get millimeters from meters
        _f_normalized[k] = 1 + _p1 * std::exp(-p2*std::pow((_r[_num_radial-1] - _r[k])*1000, p3));
        average_f += _f_normalized[k] * _volume_frac[k];
      }
      for (unsigned int k(0); k < _num_radial; ++k)
      {
        _f_normalized[k] /= average_f;
      }
    }
  }
  _is_initialized = true;
}

bool
BurnupFunction::initializeBurnup()
{
  bool initialized_burnup = false;
  if (_initial_burnup > 0 && !_initial_burnup_complete)
  {
    _initial_burnup_complete = true;
    Real time = 0;
    Real delta_t = 1e5;
    const Real fission_rate = 1e19;

    Real slope = fission_rate * _MWdperFission / _density * _fuel_molar / _hm_molar;

    bool not_done = true;
    while (not_done)
    {
      if (_average_burnup[0] + delta_t * slope > _initial_burnup)
      {
        delta_t = (_initial_burnup - _average_burnup[0])/slope;
        not_done = false;
      }
      time += delta_t;
      updateStateVariables( time );
      takeStep( delta_t, time, fission_rate, NULL, NULL, NULL );
    }
    initialized_burnup = true;
  }
  return initialized_burnup;
}

Real
BurnupFunction::p1() const
{
  Real p1 = 3.45;
  if (isParamValid("p1"))
  {
    p1 = getParam<Real>("p1");
  }
  else if (getParam<MooseEnum>("reactor_type") == "LWR")
  {
    // Default value 3.45
  }
  else if (getParam<MooseEnum>("reactor_type") == "HWR")
  {
    p1 = 2.21;
  }
  return p1;
}

std::vector<Real>
BurnupFunction::sigma_f() const
{
  std::vector<Real> sigma_f(_num_isotopes); // fission
  sigma_f[0] = 41.5;
  sigma_f[1] = 0.;
  sigma_f[2] = 105.0;
  sigma_f[3] = 0.584;
  sigma_f[4] = 120.0;
  sigma_f[5] = 0.458;
  if(_include_gadolinia){
    sigma_f[6] = 0.;
    sigma_f[7] = 0.;
  }
  if (isParamValid("sigma_f"))
  {
    sigma_f = getParam<std::vector<Real> >("sigma_f");
  }
  else if (getParam<MooseEnum>("reactor_type") == "LWR")
  {
    if(getParam<MooseEnum>("fuel_type") == "U3Si2")
    {
      sigma_f[0] = 23.9101;
      sigma_f[1] = 0.1060;
      sigma_f[2] = 56.8242;
      sigma_f[3] = 0.5939;
      sigma_f[4] = 60.2543;
      sigma_f[5] = 0.4651;
    }
  }
  else if (getParam<MooseEnum>("reactor_type") == "HWR")
  {
    sigma_f[0] = 107.95;
    sigma_f[1] = 0;
    sigma_f[2] = 239.18;
    sigma_f[3] = 0.304;
    sigma_f[4] = 296.95;
    sigma_f[5] = 0.191;
  }
  return sigma_f;
}

std::vector<Real>
BurnupFunction::sigma_c() const
{
  std::vector<Real> sigma_c(_num_isotopes);
  sigma_c[0] = 9.7;
  sigma_c[1] = 0.78;
  sigma_c[2] = 58.6;
  sigma_c[3] = 100.0;
  sigma_c[4] = 50.0;
  sigma_c[5] = 80.0;
  if(_include_gadolinia){
    sigma_c[6] = 490.0;
    sigma_c[7] = 1267.0;
  }
  if (isParamValid("sigma_c"))
  {
    sigma_c = getParam<std::vector<Real> >("sigma_c");
  }
  else if (getParam<MooseEnum>("reactor_type") == "LWR")
  {
    if(getParam<MooseEnum>("fuel_type") == "U3Si2")
    {
      sigma_c[0] = 6.4194;
      sigma_c[1] = 0.7540;
      sigma_c[2] = 32.1773;
      sigma_c[3] = 80.0460;
      sigma_c[4] = 21.6194;
      sigma_c[5] = 26.7065;
    }
  }
  else if (getParam<MooseEnum>("reactor_type") == "HWR")
  {
    sigma_c[0] = 22.3;
    sigma_c[1] = 1.16;
    sigma_c[2] = 125.36;
    sigma_c[3] = 127.26;
    sigma_c[4] = 122.41;
    sigma_c[5] = 91.3;
    if(_include_gadolinia){
      sigma_c[6] = 1471.0;
      sigma_c[7] = 3800.0;
    }
  }
  return sigma_c;
}

std::vector<Real>
BurnupFunction::sigma_a() const
{
  std::vector<Real> sigma_a(_num_isotopes);
  for (unsigned i(0); i < _num_isotopes; ++i)
    sigma_a[i] = _sigma_f[i] + _sigma_c[i];

  return sigma_a;
}

std::vector<Real>
BurnupFunction::sigma_a_thermal() const
{
  std::vector<Real> sigma_a_thermal(_num_isotopes);
  sigma_a_thermal[0] = 359.68;
  sigma_a_thermal[1] = 1.56;
  sigma_a_thermal[2] = 1207.5;
  sigma_a_thermal[3] = 193.5;
  sigma_a_thermal[4] = 1095.24;
  sigma_a_thermal[5] = 11.11;
  if(_include_gadolinia){
    sigma_a_thermal[6] = 19800.;
    sigma_a_thermal[7] = 85000.;
  }
  if (isParamValid("sigma_a_thermal"))
  {
    sigma_a_thermal = getParam<std::vector<Real> >("sigma_a_thermal");
  }
  else if (getParam<MooseEnum>("reactor_type") == "LWR")
  {
    if(getParam<MooseEnum>("fuel_type") == "U3Si2")
    {
      sigma_a_thermal[0] = 30.3248;
      sigma_a_thermal[1] = 0.8539;
      sigma_a_thermal[2] = 89.0001;
      sigma_a_thermal[3] = 80.6386;
      sigma_a_thermal[4] = 81.8657;
      sigma_a_thermal[5] = 27.1692;
    }
  }
  else if (getParam<MooseEnum>("reactor_type") == "HWR")
  {
    sigma_a_thermal[0] = 395.59;
    sigma_a_thermal[1] = 1.7;
    sigma_a_thermal[2] = 1095.7;
    sigma_a_thermal[3] = 202.2;
    sigma_a_thermal[4] = 1113.7;
    sigma_a_thermal[5] = 11.98;
    if(_include_gadolinia){
      sigma_a_thermal[6] = 23924;
      sigma_a_thermal[7] = 102477.;
    }
  }
  return sigma_a_thermal;
}

void
BurnupFunction::setup( const Real delta_t,
                       const Real time )
{
  initialize( time );

  if (time == _time_last_requested)
  {
    // The data has already been computed and is ready for interpolation.
    return;
  }

  updateStateVariables( time );

  takeStep( delta_t, time, _value, _ralp, _axial_profile, _rpf_input );
}

void
BurnupFunction::updateStateVariables( Real time )
{
  // Check if we need to update state variables.  If the solve has been cut back,
  // we must not update them.
  if (time > _time_last_requested)
  {
    // Update state variables (save new into old)
    _burnup.swap( _burnup_old );
    _fission_rate.swap( _fission_rate_old );
    if( ! _rpf_input )
    {
      _concentration.swap( _concentration_old );
      _average_concentration.swap( _average_concentration_old );
    }
  }
  _time_last_requested = time;
}

void
BurnupFunction::takeStep( Real delta_t,
                          Real time,
                          Real value,
                          Function * ralp,
                          Function * axial_profile,
                          Function * rpf_input )
{
  computeFissionRate( time, value, ralp, axial_profile );

  computeBurnup( delta_t );

  if (_rpf_active)
    computeRadialPowerFactor( time, rpf_input );

  // update fission rate according to current radial power profile
  computeFissionRate( time, value, ralp, axial_profile );

  // update burnup
  computeBurnup( delta_t );


  // Dump data to screen
//   Moose::out << std::setprecision(8)
//              << std::setw(11) << time
//              << "\t" << std::setw(11) << _concentration[0][0][_num_radial-1]
//              << "\t" << std::setw(11) << _concentration[1][0][_num_radial-1]
//              << "\t" << std::setw(11) << _concentration[2][0][_num_radial-1]
//              << "\t" << std::setw(11) << _concentration[3][0][_num_radial-1]
//              << "\t" << std::setw(11) << _concentration[4][0][_num_radial-1]
//              << "\t" << std::setw(11) << _concentration[5][0][_num_radial-1];
//   if(_has_gadolinia)
//   {
//   Moose::out << std::setprecision(8)
//              << "\t" << std::setw(11) << _concentration[6][0][_num_radial-1]
//              << "\t" << std::setw(11) << _concentration[7][0][_num_radial-1];
//   }
//   Moose::out << std::setprecision(8)
//              << "\t" << std::setw(11) << _fission_rate[0][_num_radial-1]
//              << "\t" << std::setw(11) << _burnup[0][_num_radial-1]
//              << std::endl;
}

void
BurnupFunction::computeBurnup( const Real delta_t )
{
  for (unsigned int j(0); j < _num_axial; ++j)
  {
    _average_burnup[j] = 0;

    for (unsigned int k(0); k < _num_radial; ++k)
    {
      const Real fission_rate_av = ( _fission_rate[j][k] + _fission_rate_old[j][k] ) * 0.5;
      const Real burnup_delta = delta_t * fission_rate_av * _MWdperFission / _density * _fuel_molar/_hm_molar;

      _burnup[j][k] = _burnup_old[j][k] + burnup_delta;

      _average_burnup[j] += _burnup[j][k] * _volume_frac[k];
    }
  }
}

void
BurnupFunction::computeFissionRate( const Real time, Real value, Function * ralp, Function * axial_profile )
{
  if (!ralp)
  {
    for (unsigned int j(0); j < _num_axial; ++j)
    {

      if ( axial_profile )
      {
        Point p(0, _h[j], 0);
        value *= axial_profile->value( time, p );
      }
      if (value < 0)
      {
        mooseError("Negative fission rate in BurnupFunction");
      }

      for (unsigned int k(0); k < _num_radial; ++k)
      {
        _fission_rate[j][k] = value * _rpf[j][k];
      }
    }

  }
  else
  {
    for (unsigned int j(0); j < _num_axial; ++j)
    {

      Point p(0, _h[j], 0);
      Real power = value * ralp->value( time, p );
      if ( axial_profile )
      {
        power *= axial_profile->value( time, p );
      }
      if (power < 0)
      {
        mooseError("Negative fission rate in BurnupFunction");
      }

      for (unsigned int k(0); k < _num_radial; ++k)
      {
        // _p_to_fr => fission/J/m^2
        _fission_rate[j][k] = power * _p_to_fr * _rpf[j][k];
      }
    }
  }
}

void
BurnupFunction::computeRadialPowerFactor( const Real time, Function * rpf_input )
{
  for (unsigned int j(0); j < _num_axial; ++j)
  {
    Real average_rpf(0.);

    if ( ! rpf_input )
    {
      Real f2(_hm_molar/_fuel_molar);

      std::vector<Real> fission_prob_eff(_num_radial,0.);

      for (unsigned int k(0); k < _num_radial; ++k)
      {
        for (unsigned int i(0); i < _num_isotopes; ++i)
          fission_prob_eff[k] += _sigma_f[i] * _concentration_old[i][j][k];
        const Real phi = _density / ( _MWdperFission * fission_prob_eff[k] ) * f2;

        // Burnup is in MWd/kgU
        const Real conversion_factor( phi * ( _burnup[j][k] - _burnup_old[j][k] ) );

        _concentration[N235][j][k] = _concentration_old[N235][j][k] + conversion_factor * ( - _sigma_a[N235] * _concentration_old[N235][j][k] );
        _concentration[N238][j][k] = _concentration_old[N238][j][k] + conversion_factor * ( - _sigma_a[N238] * _average_concentration_old[N238][j] * _f_normalized[k] );
        _concentration[N239][j][k] = _concentration_old[N239][j][k] + conversion_factor * ( - _sigma_a[N239] * _concentration_old[N239][j][k] + _sigma_c[N238] * _average_concentration_old[N238][j] * _f_normalized[k] );
        _concentration[N240][j][k] = _concentration_old[N240][j][k] + conversion_factor * ( - _sigma_a[N240] *_concentration_old[N240][j][k] + _sigma_c[N239] * _concentration_old[N239][j][k] );
        _concentration[N241][j][k] = _concentration_old[N241][j][k] + conversion_factor * ( - _sigma_a[N241] *_concentration_old[N241][j][k] + _sigma_c[N240] * _concentration_old[N240][j][k] );
        _concentration[N242][j][k] = _concentration_old[N242][j][k] + conversion_factor * ( - _sigma_a[N242] *_concentration_old[N242][j][k] + _sigma_c[N241] * _concentration_old[N241][j][k] );

        _concentration[N235][j][k] = std::max(_concentration[N235][j][k], 0.);
        _concentration[N238][j][k] = std::max(_concentration[N238][j][k], 0.);
        _concentration[N239][j][k] = std::max(_concentration[N239][j][k], 0.);
        _concentration[N240][j][k] = std::max(_concentration[N240][j][k], 0.);
        _concentration[N241][j][k] = std::max(_concentration[N241][j][k], 0.);
        _concentration[N242][j][k] = std::max(_concentration[N242][j][k], 0.);

        if(_include_gadolinia)
        {
          _concentration[N155][j][k] = _concentration_old[N155][j][k] + conversion_factor * ( - _sigma_a[N155] * _concentration_old[N155][j][k] );
          _concentration[N157][j][k] = _concentration_old[N157][j][k] + conversion_factor * ( - _sigma_a[N157] * _concentration_old[N157][j][k] );

          _concentration[N155][j][k] = std::max(_concentration[N155][j][k], 0.);
          _concentration[N157][j][k] = std::max(_concentration[N157][j][k], 0.);
        }
      }

      // Compute average concentrations
      for (unsigned int i(0); i < _num_isotopes; ++i)
      {
        _average_concentration[i][j] = 0;
        for (unsigned int k(0); k < _num_radial; ++k)
          _average_concentration[i][j] += _concentration[i][j][k] * _volume_frac[k];
      }

      // Compute radial power factor (RPF)
      const Real oneem28(1.e-28); // m^2/barn
      const Real param_scatt(100.);
      const Real scatt_prob_tot( param_scatt * _ntot_hm * oneem28 );
      const Real diff( 1./(3*scatt_prob_tot) );
      std::vector<Real> nflux_unnorm(_num_radial,0.);

      if( ! _include_gadolinia )
      {
        Real av_thabsor_prob_eff(0.);
        for (unsigned int i(0); i < _num_isotopes; ++i)
          av_thabsor_prob_eff += _sigma_a_thermal[i] * _average_concentration[i][j];
        av_thabsor_prob_eff *= oneem28;
        const Real av_inv_diff_length( std::sqrt(av_thabsor_prob_eff/diff) );

        analyticNeutronDiffusion( av_inv_diff_length, _r , nflux_unnorm );
      }
      else // if( _include_gadolinia )
      {
        std::vector<Real> inv_diff_length(_num_radial,0.);

        for (unsigned int k(0); k < _num_radial; ++k)
        {
          Real thabsor_prob_eff(0.);
          for (unsigned int i(0); i < _num_isotopes; ++i)
            thabsor_prob_eff += _sigma_a_thermal[i] * _concentration[i][j][k];
          thabsor_prob_eff *= oneem28;
          inv_diff_length[k] = std::sqrt(thabsor_prob_eff/diff);
        }

        finiteDifferenceNeutronDiffusion( inv_diff_length, _r , nflux_unnorm );
      }

      std::fill(fission_prob_eff.begin(),fission_prob_eff.end(),0.);

      for (unsigned int k(0); k < _num_radial; ++k)
      {
        for (unsigned int i(0); i < _num_isotopes; ++i)
          fission_prob_eff[k] += _sigma_f[i] * _concentration[i][j][k];

        _rpf[j][k] = nflux_unnorm[k] * fission_prob_eff[k];
        average_rpf += _rpf[j][k] * _volume_frac[k];
      }

      for (unsigned int k(0); k < _num_radial; ++k)
        _rpf[j][k] /= average_rpf;
    }

    else  // RPF is prescribed by input
    {
      for (unsigned int k(0); k < _num_radial; ++k)
      {
        Point p(_r[k], _h[j], 0);
        _rpf[j][k] = rpf_input->value(time,p);
        average_rpf += _rpf[j][k] * _volume_frac[k];
      }

      for (unsigned int k(0); k < _num_radial; ++k)
        _rpf[j][k] /= average_rpf;
    }
  }
}


Real
BurnupFunction::sample( const Point & pt,
                        const std::vector<Real> & radius,
                        const std::vector<Real> & height,
                        const std::vector<std::vector<Real> > & data ) // [height][radius]
{
  unsigned r1 = 0;
  unsigned r2 = 2;
  if (_axial_axis == 0)
  {
    r1 = 1;
  }
  else if (_axial_axis == 2)
  {
    r2 = 1;
  }
  else if (_axial_axis > 2)
  {
    mooseError("axial_axis must be 0, 1, or 2");
  }

  Real h = pt(_axial_axis);
  Real r = std::sqrt( pt(r1)*pt(r1) + pt(r2)*pt(r2) );

  unsigned int ax_lower, ax_upper;
  bounds( h, height, ax_lower, ax_upper );
  unsigned int r_lower, r_upper;
  bounds( r, radius, r_lower, r_upper );

  mooseAssert( ax_upper < height.size(), "Out-of-bounds index for height." );
  mooseAssert( r_upper < radius.size(), "Out-of-bounds index for radius." );

  if (ax_lower == ax_upper && r_lower == r_upper)
  {
    return data[ax_lower][r_lower];
  }
  const Real rl(radius[r_lower]);
  const Real ru(radius[r_upper]);
  const Real dll(data[ax_lower][r_lower]);
  const Real dlu(data[ax_lower][r_upper]);
  if (ax_lower == ax_upper)
  {
    return dll + (dlu-dll)*(r-rl)/(ru-rl);
  }
  const Real hl(height[ax_lower]);
  const Real hu(height[ax_upper]);
  const Real dul(data[ax_upper][r_lower]);
  if (r_lower == r_upper)
  {
    return dll + (dul-dll)*(h-hl)/(hu-hl);
  }
  const Real duu(data[ax_upper][r_upper]);
  return (dll*(hu-h)*(ru-r)+dlu*(hu-h)*(r-rl)+duu*(h-hl)*(r-rl)+dul*(h-hl)*(ru-r))/((hu-hl)*(ru-rl));
}

Real
BurnupFunction::sample( const Point & pt,
                        const std::vector<Real> & height,
                        const std::vector<Real> & data ) // [height]
{
  Real h = pt(_axial_axis);

  unsigned int ax_lower, ax_upper;
  bounds( h, height, ax_lower, ax_upper );

  mooseAssert( ax_upper < height.size(), "Out-of-bounds index for height." );

  if (ax_lower == ax_upper)
  {
    return data[ax_lower];
  }
  const Real dl(data[ax_lower]);
  const Real du(data[ax_upper]);
  const Real hl(height[ax_lower]);
  const Real hu(height[ax_upper]);
  return dl + (du-dl)*(h-hl)/(hu-hl);
}

Real
BurnupFunction::value(Real, const Point &)
{
  mooseError("BurnupFunction::value is not implemented.  Call 'burnup', 'fissionRate', 'rpf', or one of the concentration retrieval methods.");
  return 0;
}
