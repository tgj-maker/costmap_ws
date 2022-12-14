# ##############################################################################
# Setup version
cmake_minimum_required(VERSION 2.8.3)
project(virtual_costmap_layer)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)


# file(STRINGS ${CMAKE_SOURCE_DIR}/virtual_costmap_layer/VERSION version)
# message(STATUS "Read version from file: ${version}")

# cmake_minimum_required(VERSION 3.5.1)
# project(virtual_costmap_layer LANGUAGES C CXX VERSION ${version})

# set(CMAKE_CXX_STANDARD 14)
# set(CXX_STANDARD_REQUIRED ON)




# ##############################################################################
# Load and setup packages

find_package(catkin
             REQUIRED
             roscpp
             costmap_2d
             dynamic_reconfigure
             message_generation
             #custom_msgs
             tf)


# ##############################################################################
# Generate messages and configurations
generate_dynamic_reconfigure_options(cfg/VirtualLayer.cfg)

add_message_files(FILES
               Form.msg
               Obstacles.msg
               Zone.msg
               )
generate_messages(DEPENDENCIES geometry_msgs)
# ##############################################################################
# Setup catkin packages
catkin_package(INCLUDE_DIRS
               include
               LIBRARIES
               #virtual_layer
               CATKIN_DEPENDS
               costmap_2d
               dynamic_reconfigure
               message_runtime
               #custom_msgs
               tf
               DEPENDS)




# ##############################################################################
# Generate library

include_directories(include ${catkin_INCLUDE_DIRS})
add_library(virtual_layer src/virtual_layer.cpp)
target_link_libraries(virtual_layer ${catkin_LIBRARIES})
add_dependencies(virtual_layer ${PROJECT_NAME}_gencfg)
