cmake_minimum_required(VERSION 3.20)

project(
  ram
  VERSION 3.0.0
  LANGUAGES CXX
  DESCRIPTION
    "Ram is a c++ implementation of [minimap](https://github.com/lh3/minimap) with few modifications."
)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 23)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(RAM_MAIN_PROJECT ON)
endif()

option(ram_install "Generate library install target" OFF)
option(ram_build_exe "Build executable" ${RAM_MAIN_PROJECT})
option(ram_addres_sanitize
       "Build ram with ASAN support in Debug and RelWithDebInfo builds."
       ${RAM_MAIN_PROJECT})
option(ram_compile_diagnostics "Build ram with compile time diagnostics" OFF)
option(ram_profile_options "Build ram with profiling flags" OFF)

find_package(biosoup 0.10.0 QUIET)
find_package(glog 0.7.0 QUIET)
find_package(thread_pool 5.0.0 QUIET)
find_package(unordered_dense 4.4.0 QUIET)

if(NOT biosoup_FOUND)
  if(ram_install)
    message(FATAL_ERROR "Missing package biosoup!")
  endif()

  FetchContent_Declare(
    biosoup
    GIT_REPOSITORY https://github.com/rvaser/biosoup
    GIT_TAG 0.10.0)

  FetchContent_GetProperties(biosoup)
  if(NOT biosoup_POPULATED)
    FetchContent_Populate(biosoup)
    add_subdirectory(${biosoup_SOURCE_DIR} ${biosoup_BINARY_DIR}
                     EXCLUDE_FROM_ALL)
  endif()
endif()

if(NOT glog_FOUND)
  if(ram_install)
    message(FATAL_ERROR "Missing package glog!")
  endif()

  FetchContent_Declare(
    glog
    GIT_REPOSITORY https://github.com/google/glog/
    GIT_TAG v0.7.0)

  FetchContent_GetProperties(glog)
  if(NOT glog_POPULATED)
    FetchContent_Populate(glog)
    add_subdirectory(${glog_SOURCE_DIR} ${glog_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endif()

if(NOT thread_pool_FOUND)
  if(ram_install)
    message(FATAL_ERROR "Missing package thread_pool!")
  endif()

  FetchContent_Declare(
    thread_pool
    GIT_REPOSITORY https://github.com/tbrekalo/thread_pool
    GIT_TAG 5.0.0)

  FetchContent_GetProperties(thread_pool)
  if(NOT thread_pool_POPULATED)
    FetchContent_Populate(thread_pool)
    add_subdirectory(${thread_pool_SOURCE_DIR} ${thread_pool_BINARY_DIR}
                     EXCLUDE_FROM_ALL)
  endif()
endif()

if(NOT unordered_dense_FOUND)
  if(ram_install)
    message(FATAL_ERROR "Missing package unordered_dense!")
  endif()

  FetchContent_Declare(
    unordered_dense
    GIT_REPOSITORY https://github.com/martinus/unordered_dense
    GIT_TAG v4.4.0)

  FetchContent_GetProperties(unordered_dense)
  if(NOT unordered_dense_POPULATED)
    FetchContent_Populate(unordered_dense)
    add_subdirectory(${unordered_dense_SOURCE_DIR}
                     ${unordered_dense_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endif()

if(ram_build_exe)
  FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts
    GIT_TAG v3.1.1)

  FetchContent_GetProperties(cxxopts)
  if(NOT cxxopts_POPULATED)
    FetchContent_Populate(cxxopts)
    add_subdirectory(${cxxopts_SOURCE_DIR} ${cxxopts_BINARY_DIR}
                     EXCLUDE_FROM_ALL)
  endif()
endif()

if(ram_build_exe)
  find_package(bioparser 3.0.13 QUIET)
  if(NOT bioparser_FOUND)
    FetchContent_Declare(
      bioparser
      GIT_REPOSITORY https://github.com/rvaser/bioparser
      GIT_TAG 3.0.13)

    FetchContent_GetProperties(bioparser)
    if(NOT bioparser_POPULATED)
      FetchContent_Populate(bioparser)
      add_subdirectory(${bioparser_SOURCE_DIR} ${bioparser_BINARY_DIR}
                       EXCLUDE_FROM_ALL)
    endif()
  endif()
endif()

add_library(ram src/algorithm.cpp src/lcskpp.cpp src/minimizer_engine.cpp
                src/io.cpp)
add_library(ram::ram ALIAS ram)

target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include>)

target_link_libraries(ram bioparser::bioparser biosoup::biosoup glog::glog
                      thread_pool::thread_pool unordered_dense::unordered_dense)

target_compile_options(
  ram
  PRIVATE -Wall
          -Wcast-align
          -Wextra
          -Wold-style-cast
          -Wpedantic
          -Wshadow
          -Wunused
          -march=native)

if(ram_addres_sanitize)
  target_compile_options(
    ram PUBLIC $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>
               $<$<CONFIG:Debug,RelWithDebInfo>:-fno-omit-frame-pointer>)
  target_link_options(ram PUBLIC
                      $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>)
endif()

if(ram_compile_diagnostics)
  target_compile_options(
    ram
    PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wsuggest-attribute=const>
            $<$<CXX_COMPILER_ID:GNU>:-Wsuggest-attribute=noreturn>
            $<$<CXX_COMPILER_ID:GNU>:-Wsuggest-attribute=pure>
            $<$<CXX_COMPILER_ID:Clang>:-Rpass-analysis=loop-vectorize>)
endif()

if(ram_profile_options)
  target_compile_options(ram PUBLIC -O2 -g -fno-omit-frame-pointer)
endif()

if(ram_install)
  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)

  install(
    TARGETS ram
    EXPORT ${PROJECT_NAME}Targets
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY include/ram DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(
    EXPORT ${PROJECT_NAME}Targets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
endif()

if(ram_build_exe)
  add_executable(ram_exe src/main.cpp)

  target_link_libraries(ram_exe ram cxxopts::cxxopts)

  target_compile_definitions(ram_exe PRIVATE VERSION="${PROJECT_VERSION}")
  set_property(TARGET ram_exe PROPERTY OUTPUT_NAME ram)

  install(TARGETS ram_exe DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
