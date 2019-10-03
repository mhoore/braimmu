cmake_minimum_required(VERSION 3.8)
project(braimmu LANGUAGES CXX CUDA C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_C_FLAGS "-O3 -mtune=native -march=native")
# set(CMAKE_CXX_FLAGS "-O3 -mtune=native -march=native")
set(CMAKE_BUILD_TYPE Debug)

# find_package(Threads REQUIRED)
find_package(CUDA)

find_package(OpenMP)
if(OpenMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OpenMP_FOUND)

find_package(MPI)
if(MPI_CXX_COMPILE_FLAGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
endif(MPI_CXX_COMPILE_FLAGS)
if(MPI_CXX_LINK_FLAGS)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
endif(MPI_CXX_LINK_FLAGS)

include_directories(nifti/
  ${MPI_INCLUDE_PATH}
)

add_subdirectory(nifti)
add_subdirectory(src_cuda)
