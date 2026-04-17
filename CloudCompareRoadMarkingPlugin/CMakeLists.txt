cmake_minimum_required(VERSION 3.16)
project(qRoadMarkingPan2019)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# 外部构建方式（推荐）
# cmake -S CloudCompareRoadMarkingPlugin -B build-roadmarking \
#   -DCloudCompare_DIR=/path/to/CloudCompare/lib/cmake/CloudCompare

find_package(CloudCompare CONFIG REQUIRED)

find_package(Qt6 QUIET COMPONENTS Core Widgets)
if(Qt6_FOUND)
    set(QT_LIBS Qt6::Core Qt6::Widgets)
else()
    find_package(Qt5 REQUIRED COMPONENTS Core Widgets)
    set(QT_LIBS Qt5::Core Qt5::Widgets)
endif()

set(PLUGIN_SOURCES src/qRoadMarking.cpp)
set(PLUGIN_HEADERS include/qRoadMarking.h)

add_library(${PROJECT_NAME} SHARED ${PLUGIN_SOURCES} ${PLUGIN_HEADERS})

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# CloudCompare 在不同版本里的 target 名可能不同，这里做兼容处理。
foreach(candidate CCPluginAPI CCCoreLib qCC_db CloudCompare::CCPluginAPI CloudCompare::CCCoreLib CloudCompare::qCC_db)
    if(TARGET ${candidate})
        target_link_libraries(${PROJECT_NAME} PRIVATE ${candidate})
    endif()
endforeach()

target_link_libraries(${PROJECT_NAME} PRIVATE ${QT_LIBS})

set_target_properties(${PROJECT_NAME} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plugins
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plugins
)

