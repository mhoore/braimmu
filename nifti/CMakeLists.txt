set(NIFTI_SRCS
nifti1_io.cpp
znzlib.cpp
)
add_library(nifti ${NIFTI_SRCS})
