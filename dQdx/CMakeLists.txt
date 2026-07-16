# cmake .. -DCMAKE_INSTALL_PREFIX=../install

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
