SET(MODULE_NAME TOFReconstruction)

O2_SETUP(NAME ${MODULE_NAME})

set(SRCS
  src/DataReader.cxx
)

set(HEADERS
  include/TOFReconstruction/DataReader.h
)

SET(LINKDEF src/TOFReconstructionLinkDef.h)
SET(LIBRARY_NAME ${MODULE_NAME})
SET(BUCKET_NAME tof_reconstruction_bucket)

O2_GENERATE_LIBRARY()
