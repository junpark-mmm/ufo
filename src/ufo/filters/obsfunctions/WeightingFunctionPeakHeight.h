/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_WEIGHTINGFUNCTIONPEAKHEIGHT_H_
#define UFO_FILTERS_OBSFUNCTIONS_WEIGHTINGFUNCTIONPEAKHEIGHT_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

/// \brief Height conversion method for WeightingFunctionPeakHeight
enum class HeightMethod {
  HYPSOMETRIC,
  STANDARD_ATMOSPHERE
};

/// \brief Parameter traits helper for the HeightMethod enum
struct HeightMethodParameterTraitsHelper {
  typedef HeightMethod EnumType;
  static constexpr char enumTypeName[] = "HeightMethod";
  static constexpr util::NamedEnumerator<HeightMethod> namedValues[] = {
      {HeightMethod::HYPSOMETRIC, "hypsometric"},
      {HeightMethod::STANDARD_ATMOSPHERE, "standard_atmosphere"}};
};

}  // namespace ufo

namespace oops {

/// \brief Register the HeightMethod enum with the oops parameter system
template <>
struct ParameterTraits<ufo::HeightMethod>
    : public EnumParameterTraits<ufo::HeightMethodParameterTraitsHelper> {};

}  // namespace oops

namespace ufo {

// -----------------------------------------------------------------------------

/// \brief Options for WeightingFunctionPeakHeight
class WeightingFunctionPeakHeightParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(WeightingFunctionPeakHeightParameters, Parameters)

 public:
  /// Channel list, e.g. "1-15" or "1,3,5"
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Height conversion method: "hypsometric" (default) or "standard_atmosphere"
  oops::Parameter<HeightMethod> method{"method", HeightMethod::HYPSOMETRIC, this};

  /// GeoVaLs variable providing per-location surface height (m), used as the
  /// integration base in hypsometric mode.  Required when method is
  /// hypsometric; unused (and may be omitted) for standard_atmosphere.
  /// Give the GeoVaLs name without the "GeoVaLs/" prefix,
  /// e.g., "height_above_mean_sea_level_at_surface".
  oops::OptionalParameter<std::string> surfaceHeightVariable{
      "surface height variable", this};

  /// Convert the hypsometric geopotential height to geometric height (m) using
  /// latitude.  Keep true (default) for MPAS-JEDI, whose GeometryIterator
  /// Point3[2] is geometric height.  Set to false when the model-to-JEDI
  /// interface uses geopotential height as its vertical coordinate, in which
  /// case the geopotential height is written out unchanged.  Has no effect in
  /// standard_atmosphere mode.
  oops::Parameter<bool> convertToGeometricHeight{
      "convert to geometric height", true, this};
};

// -----------------------------------------------------------------------------

/// \brief Converts the CRTM diagnostic pressure_level_at_peak_of_weightingfunction
/// (a layer index stored in ObsDiag) to geometric height in metres.
///
/// Background
/// ----------
/// The CRTM ObsDiag variable "pressure_level_at_peak_of_weightingfunction" stores
/// the Fortran layer index (1-based, 1=TOA) of the weighting-function peak, not
/// the pressure itself.  This ObsFunction looks up the pressure at that level from
/// GeoVaLs/air_pressure and then converts it to height.
///
/// The output is suitable for use as iodaVerticalCoordinate in ObsVertLocalization
/// with MPAS-JEDI, whose GeometryIterator provides Point3[2] in geometric height (m).
///
/// Methods
/// -------
/// **hypsometric** (default):
///   Integrates upward from the model surface using GeoVaLs/air_pressure and
///   GeoVaLs/air_temperature:
///     z = z_sfc + (Rd/g) * Σ T_mean_k * ln(p_k / p_{k+1})
///   Most accurate; uses the actual model temperature profile.
///
/// **standard_atmosphere**:
///   Applies the ICAO ISA formula to the pressure at the peak level.
///   Troposphere  (p > 22632 Pa):  z = (T0/L)*[1-(p/p0)^(Rd*L/g)]
///   Stratosphere (p ≤ 22632 Pa):  z = 11000 + (Rd/g)*216.65*ln(22632/p)
///   Does not require GeoVaLs/air_temperature.
///
/// Example YAML
/// ------------
///   obs function:
///     name: ObsFunction/WeightingFunctionPeakHeight
///     channels: &channels 1-15
///     options:
///       channels: *channels
///       method: hypsometric     # or: standard_atmosphere
///       surface height variable: height_above_mean_sea_level_at_surface  # from GeoVaLs

class WeightingFunctionPeakHeight : public ObsFunctionBase<float> {
 public:
  explicit WeightingFunctionPeakHeight(const eckit::LocalConfiguration &);
  ~WeightingFunctionPeakHeight();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;

 private:
  WeightingFunctionPeakHeightParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_WEIGHTINGFUNCTIONPEAKHEIGHT_H_
