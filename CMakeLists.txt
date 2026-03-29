cmake_minimum_required(VERSION 3.15) # Recomendado 3.15+ para CMake moderno

# Silence CMake dev warning about removed FindBoost module
if(POLICY CMP0167)
  set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
  cmake_policy(SET CMP0167 NEW)
endif()

# Provide safe defaults when CMake is configured outside presets.
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND EXISTS "C:/vcpkg/scripts/buildsystems/vcpkg.cmake")
  set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")
endif()

if(NOT DEFINED VCPKG_TARGET_TRIPLET)
  set(VCPKG_TARGET_TRIPLET "x64-windows" CACHE STRING "Vcpkg target triplet")
endif()

project(MinimizacionCortes VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. BUSCAR LA LIBRERÍA
# Prefer package manager (vcpkg), fallback to FetchContent if not available.
find_package(raylib CONFIG QUIET)

if(NOT raylib_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    raylib
    URL https://github.com/raysan5/raylib/archive/refs/tags/5.0.tar.gz
  )
  FetchContent_MakeAvailable(raylib)
endif()

# 2. DECLARAR FUENTES
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")

# 3. CREAR EL EJECUTABLE
add_executable(${PROJECT_NAME} ${SOURCES})

# 4. CONFIGURAR EL TARGET (Incluir directorios y vincular librerías)
target_include_directories(${PROJECT_NAME}
  PUBLIC
    src
    include
)

target_link_libraries(${PROJECT_NAME} PRIVATE raylib)