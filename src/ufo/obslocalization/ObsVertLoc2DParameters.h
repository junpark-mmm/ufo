/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSVERTLOC2DPARAMETERS_H_
#define UFO_OBSLOCALIZATION_OBSVERTLOC2DPARAMETERS_H_

#include <cmath>
#include <string>
#include <vector>

#include "oops/util/missingValues.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/obslocalization/ObsLocalizationParametersBase.h"

namespace ufo {

/// \brief Options controlling 2D vertical localization, where the vertical
/// coordinate is stored per (location, ObsVector-variable) pair.
///
/// The coordinate is read from suffixed ioda variables `<group>/<name>_<idx>`
/// for each idx in `channels`. This matches what the Variable Assignment
/// filter produces when configured with a `channels:` list, e.g. when copying
/// ObsDiag/pressure_level_at_peak_of_weightingfunction into MetaData.
class ObsVertLoc2DParameters : public ObsLocalizationParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsVertLoc2DParameters, ObsLocalizationParametersBase)

 public:
  /// Scalar fallback lengthscale. Used for every channel/variable unless
  /// `vertical lengthscale per channel` is provided. Same units as the coord.
  oops::Parameter<double> lengthscale{"vertical lengthscale",
                        "lengthscale for localization beyond which localization is set to zero",
                        0.0, this};

  /// Optional per-channel lengthscales. If set, length must equal the number
  /// of entries in `channels`; this vector overrides the scalar lengthscale.
  oops::OptionalParameter<std::vector<double>> lengthscalePerChannel{
                        "vertical lengthscale per channel", this};

  oops::Parameter<bool> logTransform{"apply log transformation",
                        "apply localization function to the logarithm of the distance",
                        false, this};

  oops::Parameter<std::string> iodaVerticalCoordinateGroup{"ioda vertical coordinate group",
                       "group in the ioda file that stores vertical coordinate (default MetaData)",
                       "MetaData", this};

  /// Base variable name. Per-channel suffixes (`_1`, `_2`, ...) are appended
  /// internally to form the actual ioda variable names to read.
  oops::RequiredParameter<std::string> iodaVerticalCoordinate{"ioda vertical coordinate",
                       "base field name; the class appends `_<channel>` to read per-channel"
                       " variables", this};

  /// Channel / per-variable index list, e.g. "1-15" or "1,3,5". Parsed with
  /// oops::parseIntSet. The list length must match ObsVector::nvars().
  oops::RequiredParameter<std::string> channels{"channels",
                       "channel list whose length matches ObsVector::nvars()", this};

  oops::Parameter<std::string> localizationFunction{"localization function",
                       "localization function (Box Car, Gaspari Cohn, or SOAR)",
                       "Gaspari Cohn", this};

  oops::OptionalParameter<int> maxnobs{"max nobs",
                       "maximum number of obs to include in localization", this};

  oops::Parameter<double> SOARexpDecayH{"soar decay",
                       "soar decay",
                       util::missingValue<double>(), this};

  /// returns vertical distance between two coordinate values
  double distance(const double & vCoord1, const double & vCoord2) const {
    return std::abs(vCoord1 - vCoord2);
  }
};

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSVERTLOC2DPARAMETERS_H_
