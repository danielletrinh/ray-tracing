cmake_minimum_required(VERSION 3.5)
project(RAY_TRACER)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenMP configuration (macOS specific)
if(APPLE)
    set(LIBOMP_H /opt/homebrew/opt/libomp/include)
    set(LIBOMP_LINK /opt/homebrew/opt/libomp/lib/libomp.dylib)

    include_directories(${LIBOMP_H})
    link_libraries(${LIBOMP_LINK})

    set(CMAKE_EXE_LINKER_FLAGS "-rpath /opt/homebrew/opt/libomp/lib")
endif()

set(SOURCE_FILES
    src/raytracer/main.cpp
    src/raytracer/raytracer.cpp
    src/helpers/camera.cpp
    src/raytracer/primitive.cpp
    src/raytracer/light.cpp
    src/helpers/scene.cpp
    src/helpers/vector3.cpp
    src/helpers/color.cpp
    src/helpers/bmp.cpp
)

file(GLOB ASSET_FILES "*.txt" "*.bmp")
foreach(ASSET_FILE ${ASSET_FILES})
    configure_file(${ASSET_FILE} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
add_executable(ray_tracer ${SOURCE_FILES})

target_link_libraries(ray_tracer
    ${LIBOMP_LINK}
)

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
include(CPack)