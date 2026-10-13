cmake_minimum_required(VERSION 3.21)
project(decay_chain_simulator LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_COLOR_DIAGNOSTICS ON)

find_package(OpenGL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(GLFW3 REQUIRED IMPORTED_TARGET glfw3)

# вендоренные зависимости
add_library(glad STATIC external/glad/src/gl.c)
target_include_directories(glad PUBLIC external/glad/include)

add_library(imgui STATIC
    external/imgui/imgui.cpp
    external/imgui/imgui_draw.cpp
    external/imgui/imgui_tables.cpp
    external/imgui/imgui_widgets.cpp
    external/imgui/backends/imgui_impl_glfw.cpp
    external/imgui/backends/imgui_impl_opengl3.cpp)
target_include_directories(imgui PUBLIC external/imgui external/imgui/backends)
target_link_libraries(imgui PUBLIC PkgConfig::GLFW3)
set_property(TARGET imgui PROPERTY COMPILE_OPTIONS -w) # чужой код: без предупреждений

# ядро физики (не зависит от графики — используется и консольным тестом)
add_library(decay_core STATIC
    src/core/Chains.cpp
    src/core/BatemanAnalytic.cpp
    src/core/DecaySim.cpp
    src/core/Util.cpp)
target_include_directories(decay_core PUBLIC src)

add_executable(decay_sim
    src/app/main.cpp
    src/gfx/Renderer.cpp
    src/ui/Ui.cpp)
target_link_libraries(decay_sim PRIVATE decay_core imgui glad OpenGL::GL)
set_property(TARGET decay_sim PROPERTY COMPILE_OPTIONS -Wall -Wextra -Wpedantic)

add_executable(bateman_test tests/test_bateman.cpp)
target_link_libraries(bateman_test PRIVATE decay_core)
set_property(TARGET bateman_test PROPERTY COMPILE_OPTIONS -Wall -Wextra -Wpedantic)