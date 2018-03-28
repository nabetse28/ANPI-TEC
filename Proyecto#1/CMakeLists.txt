cmake_minimum_required(VERSION 3.9)
project(Proyecto1_ANPI)

set(CMAKE_CXX_STANDARD 11)

add_executable(Proyecto1_ANPI
        DeflationPolynomial.cpp
        Laguerre.cpp
        main.cpp
        Muller.cpp
        RootsFinder.cpp
        UserInterface.cpp)
