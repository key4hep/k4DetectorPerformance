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
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)


project(k4Performance LANGUAGES CXX)

# please keep this layout for version setting, used by the automatic tagging script
set(PACKAGE_VERSION_MAJOR 1)
set(PACKAGE_VERSION_MINOR 2)
set(PACKAGE_VERSION_PATCH 3)

set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")

# dependencies are declared like this
find_package(ROOT COMPONENTS RIO Tree)
find_package(EDM4HEP)
find_package(k4FWCore)
find_package(Gaudi)
find_package(DD4hep)
#---------------------------------------------------------------

include(cmake/Key4hepConfig.cmake)

include(GNUInstallDirs)

include(CTest)
function(set_test_env _testname)
	set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT "ROOT_INCLUDE_PATH=$<$<TARGET_EXISTS:podio::podio>:$<TARGET_FILE_DIR:podio::podio>/../include>:$<$<TARGET_EXISTS:EDM4HEP::edm4hep>:$<TARGET_FILE_DIR:EDM4HEP::edm4hep>/../include>:$ENV{ROOT_INCLUDE_PATH}")
	set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}:${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}:${PROJECT_BINARY_DIR}/genConfDir/${PackageName}:$<$<TARGET_EXISTS:ROOT::Core>:$<TARGET_FILE_DIR:ROOT::Core>>:$<$<TARGET_EXISTS:EDM4HEP::edm4hep>:$<TARGET_FILE_DIR:EDM4HEP::edm4hep>>:$<$<TARGET_EXISTS:podio::podio>:$<TARGET_FILE_DIR:podio::podio>>:$ENV{LD_LIBRARY_PATH}")
	set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}:${PROJECT_BINARY_DIR}/genConfDir:$ENV{PYTHONPATH}")
endfunction()

add_subdirectory(RecoMCTruthLinkers)
add_subdirectory(Di-ParticlePerformance/Di-ParticleAnalysis)
add_subdirectory(TrackingPerformance)

include(cmake/CreateProjectConfig.cmake)
