cmake_minimum_required(VERSION 3.10)

project(molsim Fortran)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
endif()

set(FLIBMACRO "Src/molsim.lib")
add_definitions(-DFLIBMACRO="${FLIBMACRO}")

set(CMAKE_SKIP_BUILD_RPATH            FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH    TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_Fortran_FLAGS         "-DLOCAL_GFORTRAN -DF03_CBIND -D_NORMAL_ -I/usr/include -ffree-form -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -fimplicit-none -fbacktrace -g -ffpe-trap=zero,overflow -finit-real=nan")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2")


add_definitions(-DFFTW)
if(FFTW)
    set(fftw3 ${FFTW})
    message("FFTW given as option: ${FFTW}")
else()
    message("FFTW3 not specified. Searching...")
    find_library(fftw3 NAMES fftw fftw3 fftw3)
    if(fftw3)
       message("Found FFTW3 at ${fftw3}")
    else()
       message("FFTW not found!")
    endif()
endif()

include_directories(Src)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Src)
set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR}/Bin)

# add_executable(molsim_ser.exe
# add_library(serial
set(SERIAL_FILES
    Src/bd.F90
    Src/celllist.F90
    Src/coordinate.F90
    Src/denergy.F90
    Src/dump.F90
    Src/dynamic.F90
    Src/energy.F90
    Src/gc.f90
    Src/group.F90
    Src/image.F90
    Src/mc.F90
    Src/md.F90
    Src/mesh.F90
    Src/mixed.F90
    Src/mol.F90
    Src/molaux.F90
    Src/mollib.F90
    Src/molsim.F90
    Src/molsim.lib
    Src/molsim.ver
    Src/molsim_ser
    Src/moluser.F90
    Src/nlist.F90
    Src/particle.F90
    Src/potential.F90
    Src/sso.F90
    Src/static.F90
    Src/statistics.F90
    Src/testmesh.F90
    )

add_library(serial ${SERIAL_FILES})

add_executable(molsim_ser.exe Src/molsim.F90)
target_link_libraries(molsim_ser.exe serial fftw3)
