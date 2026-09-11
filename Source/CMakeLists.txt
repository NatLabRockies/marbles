target_sources(${marbles_lib_name}
  PRIVATE
    BC.H
    Bubble.H
    Bubble.cpp
    Constants.H
    EB.H
    EB.cpp
    Factory.H
    FillPatchOps.H
    Geometry.H
    Geometry.cpp
    IC.H
    IC.cpp
    LBM.H
    LBM.cpp
    Stencil.H
    Stencil.cpp
    Utilities.H
    Utilities.cpp
    VelocityBC.H
    VelocityBC.cpp
  )

target_sources(${marbles_exe_name} PRIVATE main.cpp)

target_include_directories(${marbles_lib_name} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

include(AMReXBuildInfo)
generate_buildinfo(${marbles_lib_name} ${CMAKE_SOURCE_DIR})

configure_file("${CMAKE_SOURCE_DIR}/CMake/MarblesVersion.H.in"
  "${CMAKE_CURRENT_BINARY_DIR}/MarblesVersion.H" @ONLY)

target_link_libraries_system(${marbles_lib_name} PUBLIC AMReX::amrex)
target_link_libraries(${marbles_exe_name} PRIVATE ${marbles_lib_name})

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/MarblesVersion.H
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
