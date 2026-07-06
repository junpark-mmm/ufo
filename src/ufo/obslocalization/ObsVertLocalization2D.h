/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSVERTLOCALIZATION2D_H_
#define UFO_OBSLOCALIZATION_OBSVERTLOCALIZATION2D_H_

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <memory>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point3.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/generic/gc99.h"
#include "oops/generic/soar.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/obslocalization/ObsLocalizationBase.h"
#include "ufo/obslocalization/ObsVertLoc2DParameters.h"

namespace ufo {

/// 2D vertical observation-space localization.
///
/// Stores a (nlocs x nvars) vertical coordinate read from suffixed ioda
/// variables `<group>/<base>_<idx>` for each idx in the configured `channels`
/// list. Per-variable distances and (optionally per-variable) lengthscales
/// produce a per-(location, variable) weight written into `locvector`.
///
/// Useful for radiances (each channel has its own peak-of-weighting-function
/// pressure) and any other obs type whose ObsVector variables represent
/// distinct vertical positions.
template<class ITERATOR>
class ObsVertLocalization2D: public ObsLocalizationBase<ITERATOR> {
 public:
  ObsVertLocalization2D(const eckit::Configuration &, const ioda::ObsSpace &);

  void computeLocalization(const ITERATOR &,
                           ioda::ObsVector & locvector) const override;

 private:
  ObsVertLoc2DParameters options_;
  std::vector<int> channels_;                // size = nvars_
  std::vector<double> lengthscale_;          // size = nvars_, units same as vCoord_
  std::vector<float> vCoord_;                // size = nlocs_ * nvars_, row-major
  std::size_t nlocs_;
  std::size_t nvars_;
  double maxLengthscale_;                    // for early-out culling
  std::string distName_;

  float vCoordAt(std::size_t iloc, std::size_t jvar) const {
    return vCoord_[iloc * nvars_ + jvar];
  }

  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

template<typename ITERATOR>
ObsVertLocalization2D<ITERATOR>::ObsVertLocalization2D(const eckit::Configuration & config,
                                                       const ioda::ObsSpace & obsspace)
  : options_()
{
  options_.validateAndDeserialize(config);
  distName_ = obsspace.distribution()->name();

  // Parse channel list -> channels_.
  const std::set<int> channelset = oops::parseIntSet(options_.channels.value());
  channels_.assign(channelset.begin(), channelset.end());
  nvars_  = channels_.size();
  nlocs_  = obsspace.nlocs();

  // Build per-variable lengthscale vector.
  lengthscale_.assign(nvars_, options_.lengthscale.value());
  if (options_.lengthscalePerChannel.value() != boost::none) {
    const std::vector<double> & lpc = *options_.lengthscalePerChannel.value();
    if (lpc.size() != nvars_) {
      throw eckit::BadParameter("vertical lengthscale per channel length must equal "
                                "the number of channels");
    }
    lengthscale_ = lpc;
  }
  maxLengthscale_ = 0.0;
  for (std::size_t jvar = 0; jvar < nvars_; ++jvar) {
    if (lengthscale_[jvar] <= 0.0) {
      throw eckit::BadParameter("vertical lengthscale entries must be > 0");
    }
    maxLengthscale_ = std::max(maxLengthscale_, lengthscale_[jvar]);
  }

  // Read suffixed per-channel coordinate variables into the flat 2D buffer.
  vCoord_.assign(nlocs_ * nvars_, util::missingValue<float>());
  std::vector<float> col(nlocs_);
  for (std::size_t jvar = 0; jvar < nvars_; ++jvar) {
    const std::string varName = options_.iodaVerticalCoordinate.value()
                                + "_" + std::to_string(channels_[jvar]);
    obsspace.get_db(options_.iodaVerticalCoordinateGroup.value(), varName, col);
    for (std::size_t iloc = 0; iloc < nlocs_; ++iloc) {
      float v = col[iloc];
      if (options_.logTransform.value()) {
        if (v == 0.0f) { v = FLT_EPSILON; }
        v = std::log(v);
      }
      vCoord_[iloc * nvars_ + jvar] = v;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
void ObsVertLocalization2D<ITERATOR>::computeLocalization(const ITERATOR & i,
                                                  ioda::ObsVector & locvector) const {
  oops::Log::trace() << "ObsVertLocalization2D::computeLocalization start" << std::endl;

  if (distName_ != "Halo" && distName_ != "InefficientDistribution") {
    throw eckit::BadParameter("Can not use ObsVertLocalization2D with distribution="
                              + distName_);
  }

  // Reference vertical coordinate at the analysis point.
  eckit::geometry::Point3 refPoint = *i;
  double vCoordAtIterator = refPoint[2];
  if (options_.logTransform.value()) {
    if (vCoordAtIterator == 0.0) { vCoordAtIterator = FLT_EPSILON; }
    vCoordAtIterator = std::log(vCoordAtIterator);
  }

  const double missing = util::missingValue<double>();
  const std::size_t nvarsLoc = locvector.nvars();
  if (nvarsLoc != nvars_) {
    throw eckit::BadParameter("ObsVertLocalization2D: locvector.nvars() ("
                              + std::to_string(nvarsLoc) + ") does not match the configured "
                              "channels list size (" + std::to_string(nvars_) + ")");
  }

  // Save the incoming values; we will multiply them by the localization weight
  // for kept entries and reset to missing for excluded entries.
  ioda::ObsVector locvectorTmp(locvector);
  for (std::size_t jj = 0; jj < locvector.size(); ++jj) {
    locvector[jj] = missing;
  }

  const std::string & locFunc = options_.localizationFunction.value();
  const bool isBoxCar = (locFunc == "Box Car");
  const bool isGC     = (locFunc == "Gaspari Cohn");
  const bool isSOAR   = (locFunc == "SOAR");
  if (!isBoxCar && !isGC && !isSOAR) {
    throw eckit::BadParameter("Vertical correlation function not recognized: " + locFunc);
  }
  const double soarDecayH = options_.SOARexpDecayH.value();
  if (isSOAR && soarDecayH == util::missingValue<double>()) {
    throw eckit::BadParameter("soar decay parameter is not specified");
  }

  // Per-variable distance and weight, computed independently for each (iloc, jvar).
  for (std::size_t iloc = 0; iloc < nlocs_; ++iloc) {
    for (std::size_t jvar = 0; jvar < nvars_; ++jvar) {
      const std::size_t idx = jvar + iloc * nvars_;
      const double tmpVal = locvectorTmp[idx];
      if (tmpVal == missing) { continue; }
      const float vc = vCoordAt(iloc, jvar);
      if (vc == util::missingValue<float>()) { continue; }
      const double d = options_.distance(vCoordAtIterator, static_cast<double>(vc));
      const double L = lengthscale_[jvar];
      if (d >= L) { continue; }
      double locFactor;
      if (isBoxCar) {
        locFactor = 1.0;
      } else if (isGC) {
        locFactor = oops::gc99(d / L);
      } else {
        locFactor = oops::soar(d * soarDecayH);
      }
      locvector[idx] = locFactor * tmpVal;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename ITERATOR>
void ObsVertLocalization2D<ITERATOR>::print(std::ostream & os) const {
  os << "ObsVertLocalization2D with nvars=" << nvars_
     << ", maxLengthscale=" << maxLengthscale_
     << ", function=" << options_.localizationFunction.value() << std::endl;
}

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSVERTLOCALIZATION2D_H_
