/*
 * (C) Copyright
 *
 * Licensed under the terms of the Apache Licence Version 2.0
 * http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_GOESCLOUDWATERPATH_GOESCLOUDWATERPATH_H_
#define UFO_OPERATORS_GOESCLOUDWATERPATH_GOESCLOUDWATERPATH_H_

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variable.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/Parameter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {
  class ObsVector;
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
// Parameters for GOESCloudWaterPath
class GOESCloudWaterPathParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(GOESCloudWaterPathParameters, ObsOperatorParametersBase)

/*!
 * \brief Observation operator for GOES Cloud Water Path (CWP/LWP/IWP)
 *
 * Computes the vertically integrated cloud water path from model hydrometeor
 * mixing ratios, potentially using satellite-observed cloud base pressure (CBP) 
 * and cloud top pressure (CTP) as bounds.
 *
 * The operator supports three observation types controlled by `goes obs type`:
 * [ qh ] is optional
 *  - **CWP**: Total cloud water path (all hydrometeors: qc + qi + qr + qs + qg [+ qh])
 *  - **LWP**: Liquid water path (qc + qr). 
 *  - **IWP**: Ice water path (qi + qs + qg [+ qh])
 *
 * Integration uses the trapezoidal rule over pressure layers:
 *   path = SUM{ q[k] * (P[k] - P[k+1]) } / gravity
 *
 *  - If layer-mean RH > 0.975 (near-saturation), hofx is set to missing.
 *  - If layer-mean T < 273 K and obs type is LWP, LWP is replaced by CWP.
 *
 * ### YAML example
 *
 * \code{.yaml}
 * obs operator:
 *   name: GOESCloudWaterPath
 *   goes obs type: CWP
 *   add hail mixing ratio: false
 *   pressure top: 20000.0
 *   compute layer tk: true
 *   set maximum hofx value: 6.0
 *   reverse column order: true
 * \endcode
 *
 * ### Required MetaData (from obs netCDF):
 *  - `satcbp`: Satellite cloud base pressure (Pa)?
 *  - `satctp`: Satellite cloud top pressure (Pa)?
 *
 * Output:
 *  - The operator writes the computed water path (kg/m2) into hofx.
 */

 public:
  // GOES observation type: "CWP", "LWP", or "IWP"
  oops::Parameter<std::string> GOESOBSType{"goes obs type", "CWP", this};

  // Include hail mixing ratio in integration or not
  oops::Parameter<bool> addHailMxr{"add hail mixing ratio", false, this};

  // Top of integration column in Pa (levels above this will be zeroed out)
  oops::Parameter<float> pressureTop{"pressure top", 20000.0, this};

  // Maximum allowed hofx value (kg/m2); values above this are capped
  oops::Parameter<float> maximumHofx{"set maximum hofx value", 6.0, this};

  // Compute layer-mean temperature to check for freezing (T < 273K -> LWP becomes CWP)
  oops::Parameter<bool> computeLayerTk{"compute layer tk", true, this};

  // Reverse GeoVaL column order from top-to-bottom to bottom-to-top
  // Set to true if GeoVaLs are ordered top-to-bottom (index 0 = TOA)
  oops::Parameter<bool> reverseColumnOrder{"reverse column order", true, this};
};

// -----------------------------------------------------------------------------
// GOES Cloud Water Path observation operator
class GOESCloudWaterPath : public ObsOperatorBase {
 public:
  typedef GOESCloudWaterPathParameters Parameters_;

  static const std::string classname() {return "ufo::GOESCloudWaterPath";}

  GOESCloudWaterPath(const ioda::ObsSpace &odb, const Parameters_ &params);
  ~GOESCloudWaterPath() override;

  void simulateObs(const GeoVaLs &,
                   ioda::ObsVector &,
                   ObsDiagnostics &, const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override {return requiredVars_;}

 private:
  void print(std::ostream &) const override;

  const ioda::ObsSpace & odb_;
  std::string GOESOBSType_;
  bool addHailMxr_;
  float pressureTop_;
  bool computeLayerTk_;
  float maximumHofx_;
  bool reverseColumnOrder_;

  oops::Variables requiredVars_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GOESCLOUDWATERPATH_GOESCLOUDWATERPATH_H_
