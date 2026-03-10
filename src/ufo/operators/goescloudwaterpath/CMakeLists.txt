# (C) Copyright UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( goescloudwaterpath_files
    GOESCloudWaterPath.h
    GOESCloudWaterPath.cc
)
PREPEND( _p_goescloudwaterpath_files     "operators/goescloudwaterpath"     ${goescloudwaterpath_files} )

set ( goescloudwaterpath_src_files
      ${_p_goescloudwaterpath_files}
      PARENT_SCOPE
)
