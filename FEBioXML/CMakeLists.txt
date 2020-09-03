cmake_minimum_required(VERSION 3.10)
project(FEBio)

add_definitions(${DEF})

include_directories(${INC})

file(GLOB SOURCES "*.cpp")

add_library(febioxml ${SOURCES})
