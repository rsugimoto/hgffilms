cmake_minimum_required(VERSION 3.1)
project(main)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-Wno-c++11-narrowing")

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
# option(LIBIGL_WITH_TETGEN            "Use Tetgen"         ON)
# option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use IMGUI"          ON)

find_package(LIBIGL REQUIRED QUIET)


find_package (LAPACK REQUIRED)
add_definitions (-DHAVE_LAPACK)
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${LAPACK_LIBRARIES})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/LosTopos3D)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/newsparse)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/tunicate)

include_directories (${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory (LosTopos)

# Add your project files
file(GLOB SRCFILES *.cpp src/*cpp )
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_executable(${PROJECT_NAME} ${SRCFILES})
# target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::tetgen igl::triangle igl::opengl_glfw_imgui)
target_link_libraries(${PROJECT_NAME} LosTopos igl::core igl::opengl_glfw igl::opengl_glfw_imgui)