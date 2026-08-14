#[[
Copyright (c) 2020-2024 Key4hep-Project.

This file is part of Key4hep.
See https://key4hep.github.io/key4hep-doc/ for further info.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
]]

set(PackageName TrackingPerformance)

project(${PackageName})

# ---------------------------------------------------------------------------
# External test data
# ---------------------------------------------------------------------------
#
#   https://key4hep.web.cern.ch/testFiles/k4DetectorPerformance/<hash>
#
# The test accesses it through:
#   DATA{.../MuGuns_CLD_o3_v01_2026_07_01.root}
# ---------------------------------------------------------------------------

include(ExternalData)

list(APPEND ExternalData_URL_TEMPLATES
  "https://key4hep.web.cern.ch:443/testFiles/k4DetectorPerformance/%(hash)"
  "https://key4hep.web.cern.ch:443/testFiles/k4DetectorPerformance/MuGuns_CLD_o3_v01_2026_07_01.root"
)

# ---------------------------------------------------------------------------
# Sources and headers
# ---------------------------------------------------------------------------

file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cc
  ${PROJECT_SOURCE_DIR}/src/*.cpp
  ${PROJECT_SOURCE_DIR}/src/components/*.cpp
)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree MathCore MathMore Graf Graf3d Hist)

# ---------------------------------------------------------------------------
# TrackingPerformance module
# ---------------------------------------------------------------------------

gaudi_add_module(${PackageName}
  SOURCES ${sources}
  LINK
    Gaudi::GaudiKernel
    EDM4HEP::edm4hep
    k4FWCore::k4FWCore
    ROOT::MathCore
    ROOT::MathMore
    ROOT::Physics
    ROOT::Graf
    ROOT::Graf3d
    ROOT::Hist
)

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(${PackageName} PROPERTIES
  PUBLIC_HEADER "${headers}"
)

# ---------------------------------------------------------------------------
# Test scripts
# ---------------------------------------------------------------------------
#
# testTrackingValidation.sh:
#   lightweight CI test
#   runs only TrackingValidation on a pre-produced CLD reco file
#
# testTrackingValidationFullReco.sh:
#   manual/example script
#   runs the full IDEA pipeline
#   intentionally not registered as a CI test
# ---------------------------------------------------------------------------

file(GLOB test_python_scripts
  ${PROJECT_SOURCE_DIR}/test/*.py
)

set(test_shell_scripts
  ${PROJECT_SOURCE_DIR}/test/testTrackingValidation.sh
  ${PROJECT_SOURCE_DIR}/test/testTrackingValidationFullReco.sh
)

set(test_scripts
  ${test_python_scripts}
  ${test_shell_scripts}
)

file(COPY ${test_scripts} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test)

# ---------------------------------------------------------------------------
# Installation
# ---------------------------------------------------------------------------

install(TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}" COMPONENT dev
)

install(FILES ${test_scripts} DESTINATION test)

set(test_name "testTrackingValidation")

ExternalData_Add_Test(
  k4DetectorPerformanceData
  NAME ${test_name}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test/testTrackingValidation.sh
          DATA{${CMAKE_CURRENT_SOURCE_DIR}/test/MuGuns_CLD_o3_v01_2026_07_01.root}
)

set_test_env(${test_name})

set_tests_properties(
  ${test_name}
  PROPERTIES
    TIMEOUT 300
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

ExternalData_Add_Target(k4DetectorPerformanceData)