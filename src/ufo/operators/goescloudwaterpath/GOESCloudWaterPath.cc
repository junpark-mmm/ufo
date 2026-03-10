/*
 * (C) Copyright
 *
 * Licensed under the terms of the Apache Licence Version 2.0
 * http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/goescloudwaterpath/GOESCloudWaterPath.h"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"
#include "oops/base/Variable.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// GeoVaL variables from 'ufo_variables_mod.F90'
static const std::string var_prsi  = "air_pressure_levels"; // interface level (+1)
static const std::string var_ts   = "air_temperature";
static const std::string var_mixr = "water_vapor_mixing_ratio_wrt_dry_air";
static const std::string var_qc   = "cloud_liquid_water";
static const std::string var_qi   = "cloud_liquid_ice";
static const std::string var_qr   = "rain_water";
static const std::string var_qs   = "snow_water";
static const std::string var_qg   = "graupel";
static const std::string var_qh   = "hail";

// Physical constants from DART 'assimilation_code/modules/utilities/types_mod.f90'
static constexpr double L_over_Rv = 5418.12;     // L/Rv

// Thresholds used for layer_rh and layer_tk
static constexpr double rhThreshold = 0.975;      // near-saturation RH threshold
static constexpr double freezingTk  = 273.0;      // freezing temperature (K)

// Pressure bounds for valid integration (Pa)
static constexpr double validMinPrs  = 10000.0;   // 100 hPa
static constexpr double validMaxPrs  = 110000.0;  // 1100 hPa

// -----------------------------------------------------------------------------
static ObsOperatorMaker<GOESCloudWaterPath> maker_("GOESCloudWaterPath");

// -----------------------------------------------------------------------------
GOESCloudWaterPath::GOESCloudWaterPath(const ioda::ObsSpace & odb,
                                       const Parameters_ & params)
  : ObsOperatorBase(odb, VariableNameMap(params.AliasFile.value())),
    odb_(odb),
    GOESOBSType_(params.GOESOBSType.value()),
    addHailMxr_(params.addHailMxr.value()),
    pressureTop_(params.pressureTop.value()),
    computeLayerTk_(params.computeLayerTk.value()),
    maximumHofx_(params.maximumHofx.value()),
    reverseColumnOrder_(params.reverseColumnOrder.value())
{
  // Always require pressure levels, temperature, and qv
  requiredVars_ += oops::Variables(std::vector<oops::Variable>{
                   oops::Variable(var_prsi),
                   oops::Variable(var_mixr),
                   oops::Variable(var_ts)});

  // Hydrometeors based on goes obs type
  if (GOESOBSType_ == "CWP") {
    requiredVars_ += oops::Variables(std::vector<oops::Variable>{
                     oops::Variable(var_qc), oops::Variable(var_qi),
                     oops::Variable(var_qr), oops::Variable(var_qs),
                     oops::Variable(var_qg)});
  } else if (GOESOBSType_ == "LWP") {
    requiredVars_ += oops::Variables(std::vector<oops::Variable>{
                     oops::Variable(var_qc), oops::Variable(var_qr)});
    // LWP needs qi, qs, qg for CWP computation (if check freezing)
    if (computeLayerTk_) {
      requiredVars_ += oops::Variables(std::vector<oops::Variable>{
                       oops::Variable(var_qi), oops::Variable(var_qs),
                       oops::Variable(var_qg)});
    }
  } else if (GOESOBSType_ == "IWP") {
    requiredVars_ += oops::Variables(std::vector<oops::Variable>{
                     oops::Variable(var_qi), oops::Variable(var_qs),
                     oops::Variable(var_qg)});
  } else {
    throw eckit::BadParameter(
      "GOESCloudWaterPath: Unknown goes obs type: '" + GOESOBSType_ +
      "'. Must be CWP, LWP, or IWP.", Here());
  }

  // Optional hail mixing ratio
  if (addHailMxr_) {
    requiredVars_ += oops::Variables(std::vector<oops::Variable>{
                     oops::Variable(var_qh)});
  }

  oops::Log::trace() << "GOESCloudWaterPath constructor done. "
                     << "Type=" << GOESOBSType_
                     << ", addHail=" << addHailMxr_
                     << ", pressureTop=" << pressureTop_
                     << ", layerTk=" << computeLayerTk_
                     << ", maxHofx=" << maximumHofx_
                     << ", reverseCol=" << reverseColumnOrder_ << std::endl;
}

// -----------------------------------------------------------------------------
GOESCloudWaterPath::~GOESCloudWaterPath() {
  oops::Log::trace() << "GOESCloudWaterPath destructor done" << std::endl;
}

// -----------------------------------------------------------------------------
void GOESCloudWaterPath::simulateObs(const GeoVaLs & geovals,
                                     ioda::ObsVector & hofx,
                                     ObsDiagnostics &,
                                     const QCFlags_t &) const {
  oops::Log::trace() << "GOESCloudWaterPath::simulateObs start" << std::endl;

  const std::size_t nlocs = geovals.nlocs();
  const std::size_t nlevsi = geovals.nlevs(oops::Variable{var_prsi}); // interface
  const std::size_t nlevs = nlevsi - 1; // mid-layer
  ASSERT(nlocs == hofx.nlocs());
  ASSERT(nlocs == odb_.nlocs());

  const double missing = util::missingValue<double>();
  const double pTop = static_cast<double>(pressureTop_);

  // Read satellite cloud base/top pressure from MetaData
  std::vector<float> satcbpVec(nlocs, 0.0), satctpVec(nlocs, 0.0);
  odb_.get_db("MetaData", "satcbp", satcbpVec);
  odb_.get_db("MetaData", "satctp", satctpVec);

  // Determine which hydrometeor profiles are needed for goes obs type
  const bool needQc = (GOESOBSType_ == "CWP" || GOESOBSType_ == "LWP");
  const bool needQi = (GOESOBSType_ == "CWP" || GOESOBSType_ == "IWP" ||
                       (GOESOBSType_ == "LWP" && computeLayerTk_));
  const bool needQr = (GOESOBSType_ == "CWP" || GOESOBSType_ == "LWP");
  const bool needQs = (GOESOBSType_ == "CWP" || GOESOBSType_ == "IWP" ||
                       (GOESOBSType_ == "LWP" && computeLayerTk_));
  const bool needQg = (GOESOBSType_ == "CWP" || GOESOBSType_ == "IWP" ||
                       (GOESOBSType_ == "LWP" && computeLayerTk_));
  const bool needQh = addHailMxr_;

  // Loop over observation locations
  for (std::size_t loc = 0; loc < nlocs; ++loc) {
    hofx[loc] = missing;

    const float satcbp = static_cast<float>(satcbpVec[loc]);
    const float satctp = static_cast<float>(satctpVec[loc]);

    // Read pressure profile at this location
    std::vector<double> press(nlevsi);
    geovals.getAtLocation(press, oops::Variable{var_prsi}, loc);
    // after reverse column order: index 0 = surface (high P), nlevsi-1 = top (low P)
    if (reverseColumnOrder_) std::reverse(press.begin(), press.end());

    // Read hydrometeor profiles (initialized to zero)
    std::vector<double> qc(nlevs, 0.0), qi(nlevs, 0.0), qr(nlevs, 0.0);
    std::vector<double> qs(nlevs, 0.0), qg(nlevs, 0.0), qh(nlevs, 0.0);

    if (needQc) {
        geovals.getAtLocation(qc, oops::Variable{var_qc}, loc);
        if (reverseColumnOrder_) std::reverse(qc.begin(), qc.end());
    }
    if (needQi) {
        geovals.getAtLocation(qi, oops::Variable{var_qi}, loc);
        if (reverseColumnOrder_) std::reverse(qi.begin(), qi.end());
    }
    if (needQr) {
        geovals.getAtLocation(qr, oops::Variable{var_qr}, loc);
        if (reverseColumnOrder_) std::reverse(qr.begin(), qr.end());
    }
    if (needQs) {
        geovals.getAtLocation(qs, oops::Variable{var_qs}, loc);
        if (reverseColumnOrder_) std::reverse(qs.begin(), qs.end());
    }
    if (needQg) {
        geovals.getAtLocation(qg, oops::Variable{var_qg}, loc);
        if (reverseColumnOrder_) std::reverse(qg.begin(), qg.end());
    }
    if (needQh) {
        geovals.getAtLocation(qh, oops::Variable{var_qh}, loc);
        if (reverseColumnOrder_) std::reverse(qh.begin(), qh.end());
    }

    // Read T and qv
    std::vector<double> qv(nlevs, 0.0), tmpk(nlevs, 0.0);
    geovals.getAtLocation(qv, oops::Variable{var_mixr}, loc);
    if (reverseColumnOrder_) std::reverse(qv.begin(), qv.end());
    // Convert qv from g/kg to kg/kg for MPAS
    for (std::size_t k = 0; k < nlevs; ++k) {
      qv[k] *= 0.001;
    }
    geovals.getAtLocation(tmpk, oops::Variable{var_ts}, loc);
    if (reverseColumnOrder_) std::reverse(tmpk.begin(), tmpk.end());

    // Clip negative hydrometeors and zero out levels above pressureTop
    // hydrometeors: nlevs | prsi: nlevs+1
    for (std::size_t k = 0; k < nlevs; ++k) {
      if (qc[k] < 0.0 || press[k+1] <= pTop) qc[k] = 0.0;
      if (qi[k] < 0.0 || press[k+1] <= pTop) qi[k] = 0.0;
      if (qr[k] < 0.0 || press[k+1] <= pTop) qr[k] = 0.0;
      if (qs[k] < 0.0 || press[k+1] <= pTop) qs[k] = 0.0;
      if (qg[k] < 0.0 || press[k+1] <= pTop) qg[k] = 0.0;
      if (qh[k] < 0.0 || press[k+1] <= pTop) qh[k] = 0.0;
    }

    // ---- Determine bounds from satellite CBP/CTP ----
    // Pressure profile is assumed as bottom-to-top (index 0 = surface, nlevs-1 = top)
    std::size_t mink = 0;
    std::size_t maxk = nlevs - 1;

    // in the case of CWP, it doesn't compute the below.
    if (satcbp > 0.0 && satctp > 0.0) {
      bool foundBase = false;
      bool foundTop  = false;

      for (std::size_t k = 0; k < nlevs; ++k) {
        if (press[k] <= validMinPrs) continue;

        // Cloud base: first level where pressure < satcbp (scanning from high to low P)
        if (!foundBase && satcbp > press[k]) {
          mink = k;
          foundBase = true;
        }
        // Cloud top: first level where pressure < satctp
        if (!foundTop && satctp > press[k]) {
          maxk = k;
          foundTop = true;
        }
      }
    }

    // Need at least two levels to integrate
    if (maxk <= mink) {
      hofx[loc] = missing;
      continue;
    }

    // ---- Compute RH at each level ----
    // may use RH from GeoVal
    std::vector<double> rh(nlevs, 0.0);
    for (std::size_t k = 0; k < nlevs; ++k) {
      if (tmpk[k] > 0.0 && press[k] > 0.0) {
        double es   = 611.0 * std::exp(L_over_Rv * (1.0 / 273.0 - 1.0 / tmpk[k]));
        double qsat = 0.622 * es / (0.5 * (press[k] + press[k+1]) - es);
        rh[k] = (qsat > 0.0) ? qv[k] / qsat : 0.0;
        if (rh[k] < 0.0) rh[k] = 0.0001;
      }
    }

    // ---- Integration of CWP from mink to maxk ----
    double cwp = 0.0;
    double lwp = 0.0;
    double iwp = 0.0;

    for (std::size_t k = mink; k < maxk; ++k) {
      // Validate pressure pair
      if (press[k] <= validMinPrs || press[k + 1] <= validMinPrs) continue;
      if (press[k] <= press[k + 1]) continue;  // pressure must decrease with height
      if (press[k] >= validMaxPrs) continue;

      double dP = press[k] - press[k + 1];

      if (GOESOBSType_ == "CWP") {
        double sumQ = qc[k] + qi[k] + qr[k] + qs[k] + qg[k] + qh[k];
        cwp += sumQ * dP;

      } else if (GOESOBSType_ == "LWP") {
        // Liquid species
        double sumQ = qc[k] + qr[k];
        lwp += sumQ * dP;

        // Also compute CWP if freezing check is enabled (CWP -> LWP)
        if (computeLayerTk_) {
          double sumQ = qc[k] + qi[k] + qr[k] + qs[k] + qg[k] + qh[k];
          cwp += sumQ * dP;
        }

      } else if (GOESOBSType_ == "IWP") {
        // Ice species
        double sumQ = qi[k] + qs[k] + qg[k] + qh[k];
        iwp += sumQ * dP;
      }
    }

    // ---- Layer-mean RH check (near-saturation -> skip obs) ----
    {
      double rhSum = 0.0;
      std::size_t rhCount = 0;
      for (std::size_t k = mink; k < maxk; ++k) {
        rhSum += rh[k];
        ++rhCount;
      }
      double layerRH = (rhCount > 0) ? rhSum / static_cast<double>(rhCount) : 0.0;

      if (layerRH > rhThreshold) {
        if ( GOESOBSType_ == "LWP" ) {
          lwp = missing;
	}
        if ( GOESOBSType_ == "IWP" ) {
          iwp = missing;
	}
	// This is practically not meaningful since CWP computed from bottom to 200 mb
        if ( GOESOBSType_ == "CWP" ) {
          cwp = missing;
	}
      }
    }

    // ---- Layer-mean temperature check (freezing: LWP -> CWP) ----
    if (computeLayerTk_ && GOESOBSType_ == "LWP") {
      double tkSum = 0.0;
      std::size_t tkCount = 0;
      for (std::size_t k = mink; k < maxk; ++k) {
        tkSum += tmpk[k];
        ++tkCount;
      }
      double layerTk = (tkCount > 0) ? tkSum / static_cast<double>(tkCount) : 300.0;

      if (layerTk < freezingTk) {
        // Below freezing: replace LWP with CWP
        lwp = cwp;
      }
    }

    // ---- Assign output based on obs type ----
    double pathValue = 0.0;
    if (GOESOBSType_ == "CWP") {
      pathValue = cwp;
    } else if (GOESOBSType_ == "LWP") {
      pathValue = lwp;
    } else if (GOESOBSType_ == "IWP") {
      pathValue = iwp;
    }

    if (pathValue != missing) {
      // Convert to kg/m2
      pathValue /= Constants::grav;

      // Clip negatives
      if (pathValue < 0.0) pathValue = 0.0;

      // Cap at maximum value
      if (pathValue > static_cast<double>(maximumHofx_)) {
        pathValue = static_cast<double>(maximumHofx_);
      }

      hofx[loc] = pathValue;
    } else {
      hofx[loc] = missing;
    }
  }

  oops::Log::trace() << "GOESCloudWaterPath::simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------
void GOESCloudWaterPath::print(std::ostream & os) const {
  os << "GOESCloudWaterPath: type=" << GOESOBSType_
     << ", addHail=" << addHailMxr_
     << ", pressureTop=" << pressureTop_
     << ", layerTk=" << computeLayerTk_
     << ", maxHofx=" << maximumHofx_
     << ", reverseCol=" << reverseColumnOrder_;
}

// -----------------------------------------------------------------------------
}  // namespace ufo
