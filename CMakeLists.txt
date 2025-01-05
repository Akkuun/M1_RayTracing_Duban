cmake_minimum_required(VERSION 3.29)
project(Projet_Ray_Tracing)

set(CMAKE_CXX_STANDARD 14)

include_directories(HAI719I_Raytracer/src)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

add_executable(Projet_Ray_Tracing
        HAI719I_Raytracer/img/normalMaps/n1.ppm
        HAI719I_Raytracer/img/normalMaps/n2.ppm
        HAI719I_Raytracer/img/normalMaps/n3.ppm
        HAI719I_Raytracer/img/normalMaps/n4.ppm
        HAI719I_Raytracer/img/sphereTextures/s1.ppm
        HAI719I_Raytracer/img/sphereTextures/s2.ppm
        HAI719I_Raytracer/img/sphereTextures/s3.ppm
        HAI719I_Raytracer/img/sphereTextures/s4.ppm
        HAI719I_Raytracer/img/sphereTextures/s5.ppm
        HAI719I_Raytracer/img/sphereTextures/s6.ppm
        HAI719I_Raytracer/img/sphereTextures/s7.ppm
        HAI719I_Raytracer/src/Camera.cpp
        HAI719I_Raytracer/src/Camera.h
        HAI719I_Raytracer/src/Line.h
        HAI719I_Raytracer/src/Material.h
        HAI719I_Raytracer/src/matrixUtilities.h
        HAI719I_Raytracer/src/Mesh.cpp
        HAI719I_Raytracer/src/Mesh.h
        HAI719I_Raytracer/src/Plane.h
        HAI719I_Raytracer/src/Ray.h
        HAI719I_Raytracer/src/Scene.h
        HAI719I_Raytracer/src/Sphere.h
        HAI719I_Raytracer/src/Square.h
        HAI719I_Raytracer/src/Trackball.cpp
        HAI719I_Raytracer/src/Trackball.h
        HAI719I_Raytracer/src/Triangle.h
        HAI719I_Raytracer/src/Vec3.h
        HAI719I_Raytracer/main.cpp
        HAI719I_Raytracer/Makefile
        HAI719I_Raytracer/src/image_utils.cpp
)

target_link_libraries(Projet_Ray_Tracing ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})