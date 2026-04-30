cmake_minimum_required(VERSION 3.15) # Recomendado 3.15+ para CMake moderno

# Silence CMake dev warning about removed FindBoost module
if(POLICY CMP0167)
  set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
  cmake_policy(SET CMP0167 NEW)
endif()

# Provide safe defaults when CMake is configured outside presets.
if(WIN32)
  if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND EXISTS "C:/vcpkg/scripts/buildsystems/vcpkg.cmake")
    set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")
  endif()

  if(NOT DEFINED VCPKG_TARGET_TRIPLET)
    set(VCPKG_TARGET_TRIPLET "x64-windows" CACHE STRING "Vcpkg target triplet")
  endif()
endif()

project(MinimizacionCortes VERSION 1.0 LANGUAGES CXX)

include(CheckIPOSupported)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(MINIMIZACIONCORTES_ENABLE_NATIVE_OPTIMIZATIONS "Enable CPU-specific optimization flags where available" OFF)
option(MINIMIZACIONCORTES_ALLOW_FETCHCONTENT "Allow downloading dependencies via FetchContent when missing" ON)

check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
if(ipo_supported)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()

# 1. BUSCAR LA LIBRERÍA
# Prefer package manager (vcpkg), fallback to FetchContent if not available.
find_package(raylib CONFIG QUIET)

if(NOT raylib_FOUND)
  if(MINIMIZACIONCORTES_ALLOW_FETCHCONTENT)
    include(FetchContent)
    FetchContent_Declare(
      raylib
      URL https://github.com/raysan5/raylib/archive/refs/tags/5.0.tar.gz
    )
    FetchContent_MakeAvailable(raylib)
  else()
    message(FATAL_ERROR "raylib was not found. Install raylib or configure with -DMINIMIZACIONCORTES_ALLOW_FETCHCONTENT=ON.")
  endif()
endif()

# 2. DECLARAR FUENTES
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")

# Exclude old OGDF-based FruchtermanReingold files (using custom implementation instead)
list(FILTER SOURCES EXCLUDE REGEX "FruchtermanReingold\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX "ogdf-master/")

# 3. CREAR EL EJECUTABLE
add_executable(${PROJECT_NAME} ${SOURCES})

# 4. CONFIGURAR EL TARGET (Incluir directorios y vincular librerías)
target_include_directories(${PROJECT_NAME}
  PUBLIC
    src
    include
)

target_link_libraries(${PROJECT_NAME} PRIVATE raylib)

target_compile_definitions(${PROJECT_NAME}
  PRIVATE
    $<$<CONFIG:Release>:NDEBUG>
)

if(MSVC)
  target_compile_options(${PROJECT_NAME}
    PRIVATE
      $<$<CONFIG:Release>:/O2>
      $<$<CONFIG:Release>:/Ob2>
      $<$<CONFIG:Release>:/Ot>
  )

  if(MINIMIZACIONCORTES_ENABLE_NATIVE_OPTIMIZATIONS)
    target_compile_options(${PROJECT_NAME}
      PRIVATE
        $<$<CONFIG:Release>:/favor:INTEL64>
    )
  endif()
else()
  target_compile_options(${PROJECT_NAME}
    PRIVATE
      $<$<CONFIG:Release>:-O3>
  )

  if(MINIMIZACIONCORTES_ENABLE_NATIVE_OPTIMIZATIONS)
    target_compile_options(${PROJECT_NAME}
      PRIVATE
        $<$<CONFIG:Release>:-march=native>
    )
  endif()
endif()