set(BRAIMMU_SRCS
  comm.cpp
  init.cpp
  input.cpp
  main.cpp
  output.cpp
  region.cpp
  scenario_connectome.cpp
  ScenarioConnectomeStrategyCPU.cpp
  scenario_geometry.cpp
)

if(CUDA_FOUND)
  set(BRAIMMU_SRCS ${BRAIMMU_SRCS}
    ScenarioConnectomeStrategyCUDA.cu
  )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_CUDA")
endif(CUDA_FOUND)

add_executable(braimmu ${BRAIMMU_SRCS})
target_link_libraries(braimmu nifti ${MPI_CXX_LIBRARIES})
