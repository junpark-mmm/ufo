/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/WeightingFunctionPeakHeight.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ufo {

constexpr char HeightMethodParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<HeightMethod>
    HeightMethodParameterTraitsHelper::namedValues[];

static ObsFunctionMaker<WeightingFunctionPeakHeight>
       makerWeightingFunctionPeakHeight_("WeightingFunctionPeakHeight");

// -----------------------------------------------------------------------------

WeightingFunctionPeakHeight::WeightingFunctionPeakHeight(
    const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::trace() << "WeightingFunctionPeakHeight constructor start" << std::endl;
  options_.deserialize(conf);

  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // CRTM diagnostic: Fortran layer index of peak weighting function (1-based, 1=TOA)
  invars_ += Variable("ObsDiag/pressure_level_at_peak_of_weightingfunction", channels_);
  // Full pressure column
  invars_ += Variable("GeoVaLs/air_pressure");
  // Full temperature column and surface height (hypsometric mode only)
  if (options_.method.value() == HeightMethod::HYPSOMETRIC) {
    // The surface height GeoVaLs is the integration base, required in this mode.
    ASSERT(options_.surfaceHeightVariable.value() != boost::none);
    invars_ += Variable("GeoVaLs/air_temperature");
    // Per-location surface height from GeoVaLs (integration base)
    invars_ += Variable("GeoVaLs/" + options_.surfaceHeightVariable.value().get());
    // Latitude is only needed for the geopotential-to-geometric correction
    if (options_.convertToGeometricHeight.value())
      invars_ += Variable("MetaData/latitude");
  }

  oops::Log::trace() << "WeightingFunctionPeakHeight constructor complete" << std::endl;
}

// -----------------------------------------------------------------------------

WeightingFunctionPeakHeight::~WeightingFunctionPeakHeight() {
  oops::Log::trace() << "WeightingFunctionPeakHeight destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void WeightingFunctionPeakHeight::compute(const ObsFilterData & in,
                                          ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "WeightingFunctionPeakHeight compute start" << std::endl;

  const size_t nlocs  = in.nlocs();
  const size_t nchans = channels_.size();
  const size_t nlevs  = in.nlevs(Variable("GeoVaLs/air_pressure"));
  const float  miss   = util::missingValue<float>();
  const bool   hyps   = (options_.method.value() == HeightMethod::HYPSOMETRIC);
  const bool   convert = options_.convertToGeometricHeight.value();

  // ── pressure profile ──────────────────────────────────────────────────────
  // UFO GeoVaLs convention: level 0 = TOA (lowest p), level nlevs-1 = surface
  std::vector<std::vector<float>> prs(nlevs, std::vector<float>(nlocs));
  for (size_t ilev = 0; ilev < nlevs; ++ilev)
    in.get(Variable("GeoVaLs/air_pressure"), ilev, prs[ilev]);

  // ── hypsometric-only inputs: temperature profile, surface height, latitude ─
  std::vector<std::vector<float>> tmp;
  std::vector<float> z_sfc_vec;
  std::vector<float> lats;
  if (hyps) {
    tmp.resize(nlevs, std::vector<float>(nlocs));
    for (size_t ilev = 0; ilev < nlevs; ++ilev)
      in.get(Variable("GeoVaLs/air_temperature"), ilev, tmp[ilev]);
    // Per-location surface height from GeoVaLs (integration base, assumed geometric)
    z_sfc_vec.resize(nlocs);
    in.get(Variable("GeoVaLs/" + options_.surfaceHeightVariable.value().get()),
           0, z_sfc_vec);
    // Latitude only needed for the geopotential-to-geometric height correction
    if (convert) {
      lats.resize(nlocs);
      in.get(Variable("MetaData/latitude"), lats);
    }
  }

  // ── per-channel loop ──────────────────────────────────────────────────────
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    // ObsDiag stores Fortran jlevel (1-based, 1=TOA) as a scalar float per location
    std::vector<float> peak_idx(nlocs, miss);
    in.get(Variable("ObsDiag/pressure_level_at_peak_of_weightingfunction",
                    channels_)[ichan], peak_idx);

    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      out[ichan][iloc] = miss;
      if (peak_idx[iloc] == miss) continue;

      // Convert: Fortran 1-based TOA-first → C++ 0-based TOA-first
      //   jlevel=1 (TOA)        → cpp_peak = 0
      //   jlevel=nlevs (surface) → cpp_peak = nlevs-1
      const int cpp_peak = static_cast<int>(peak_idx[iloc]) - 1;
      if (cpp_peak < 0 || cpp_peak >= static_cast<int>(nlevs)) continue;

      if (hyps) {
        // ── hypsometric: integrate from surface (k=nlevs-1) upward to cpp_peak
        // prs[k] increases with k (level 0=TOA, nlevs-1=surface).  The sum
        // (Rd/g)·T̄·ln(p_lo/p_hi) yields geopotential height above the surface.
        float z = z_sfc_vec[iloc];
        bool valid = true;
        for (int k = static_cast<int>(nlevs) - 1; k > cpp_peak; --k) {
          const float p_lo = prs[k][iloc];
          const float p_hi = prs[k - 1][iloc];
          const float T_lo = tmp[k][iloc];
          const float T_hi = tmp[k - 1][iloc];
          if (p_lo == miss || p_hi == miss || T_lo == miss || T_hi == miss
              || p_lo <= 0.0f || p_hi <= 0.0f) {
            valid = false;
            break;
          }
          z += static_cast<float>(Constants::rd_over_g)
               * 0.5f * (T_lo + T_hi)
               * std::log(p_lo / p_hi);
        }
        // Optionally convert geopotential to geometric height (m) so the output
        // matches the MPAS-JEDI GeometryIterator Point3[2] used by
        // ObsVertLocalization.  When the interface's vertical coordinate is
        // geopotential height, write the geopotential height unchanged.
        if (valid)
          out[ichan][iloc] = convert
              ? formulas::Geopotential_to_Geometric_Height(lats[iloc], z)
              : z;

      } else {
        // ── ICAO Standard Atmosphere ─────────────────────────────────────────
        // Reuse the shared variable-transform formula (input pressure in Pa,
        // returns geometric height in m and handles missing/non-positive p).
        out[ichan][iloc] = formulas::Pressure_To_Height(
            prs[cpp_peak][iloc], formulas::Formulation::ICAO);
      }
    }
  }
  oops::Log::trace() << "WeightingFunctionPeakHeight compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & WeightingFunctionPeakHeight::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
