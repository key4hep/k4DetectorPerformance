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

set(PackageName dQdx)

project(${PackageName})

file(GLOB sources ${PROJECT_SOURCE_DIR}/components/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/components/*.h)


gaudi_add_module(${PackageName}
  SOURCES ${sources}
  LINK
  Gaudi::GaudiKernel
  k4FWCore::k4FWCore
  EDM4HEP::edm4hep
  DD4hep::DDRec
  ROOT::Graf3d
)

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/components>
  $<INSTALL_INTERFACE:include>
  ${MarlinUtil_INCLUDE_DIRS}
)

install(TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}"
  COMPONENT dev)


set_target_properties(${PackageName} PROPERTIES
  PUBLIC_HEADER "${headers}"
)

SET(test_name "test_rundQdxAnalyser")
ADD_TEST(NAME ${test_name} COMMAND k4run test/rundQdxAnalyser.py)
set_test_env(${test_name})
