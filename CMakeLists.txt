cmake_minimum_required(VERSION 3.15) # Recomendado 3.15+ para CGAL moderno

# Silence CMake dev warning about removed FindBoost module
if(POLICY CMP0167)
  set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
  cmake_policy(SET CMP0167 NEW)
endif()

project(MinimizacionCortes VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Suppress CGAL warning about Debug CMAKE_BUILD_TYPE
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

# 1. BUSCAR LA LIBRERÍA PRIMERO
# CGAL necesita encontrar sus dependencias (GMP, MPFR)
find_package(CGAL REQUIRED)

# 2. DECLARAR FUENTES
file(GLOB SOURCES "src/*.cpp" "main.cpp")

# 3. CREAR EL EJECUTABLE
add_executable(${PROJECT_NAME} ${SOURCES})

# 4. CONFIGURAR EL TARGET (Incluir directorios y vincular librerías)
target_include_directories(${PROJECT_NAME} PUBLIC include)

# Vincular CGAL al ejecutable
# Esto añade automáticamente los headers y las flags necesarias
target_link_libraries(${PROJECT_NAME} PRIVATE CGAL::CGAL)