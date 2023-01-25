# Thirdparty libraries
include(FetchContent)
include(ExternalProject)


# Log directory for external content builds
make_directory(${CMAKE_CURRENT_BINARY_DIR}/logs)
set(INSTALL_FILE  ${CMAKE_CURRENT_BINARY_DIR}/install.cmake)

# Set environment

list(APPEND CMAKE_LIBRARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/lib ")
list(APPEND CMAKE_INCLUDE_PATH "${CMAKE_CURRENT_BINARY_DIR}/include ")
list(APPEND CMAKE_PREFIX_PATH  "${CMAKE_CURRENT_BINARY_DIR}")

# External Project resources
# ------------------------------------------------------------------------------

# DIVSUFSORT
find_library(DIVSUFSORT_LIB divsufsort HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
if(NOT DIVSUFSORT_LIB)
    message(STATUS "divsufsort library not found. Building as an external content")
else()
    message(STATUS "divsufsort library found at ${DIVSUFSORT_LIB}.")
endif()

# DIVSUFSORT64
find_library(DIVSUFSORT64_LIB divsufsort64 HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
if(NOT DIVSUFSORT64_LIB)
    message(STATUS "divsufsort64 library not found. Building as an external content")
else()
    message(STATUS "divsufsort64 library found at ${DIVSUFSORT64_LIB}.")
endif()

# SDSL
find_library(SDSL_LIB sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
find_path(SDSL_SRC sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/include PATHS ${CONDA_PREFIX}/include)
if(NOT SDSL_LIB)
    message(STATUS "sdsl library not found. Building as an external content")
    FetchContent_Declare(sdsl
            GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
            # GIT_TAG "v2.1.1"
    )

    if(NOT sdsl_POPULATED)
        FetchContent_Populate(sdsl)
        # Install in the current binary dir
        execute_process(COMMAND ./install.sh ${CMAKE_CURRENT_BINARY_DIR}
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${sdsl_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log)
        file(APPEND ${INSTALL_FILE} "MESSAGE(STATUS \"Installing sdsl.\")\n")
        file(APPEND ${INSTALL_FILE} "execute_process(COMMAND ./install.sh \${CMAKE_INSTALL_PREFIX}
                WORKING_DIRECTORY ${sdsl_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log)\n")
        if(result)
            message(FATAL_ERROR "CMake step for sdsl-lite failed: ${result}")
        endif()
    endif()

    set(SDSL_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(DIVSUFSORT_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(SDSL_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libsdsl.a)
    set(DIVSUFSORT_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort.a)
    set(DIVSUFSORT64_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort64.a)

else()
    message(STATUS "sdsl library found at ${SDSL_LIB}.")
    message(STATUS "sdsl sources found at ${SDSL_SRC}.")
endif()

add_library(sdsl STATIC IMPORTED GLOBAL)
set_property(TARGET sdsl PROPERTY IMPORTED_LOCATION ${SDSL_LIB})
target_include_directories(sdsl INTERFACE ${SDSL_SRC})

add_library(divsufsort STATIC IMPORTED GLOBAL)
set_property(TARGET divsufsort PROPERTY IMPORTED_LOCATION ${DIVSUFSORT_LIB})
    
add_library(divsufsort64 STATIC IMPORTED GLOBAL)
set_property(TARGET divsufsort64 PROPERTY IMPORTED_LOCATION ${DIVSUFSORT64_LIB})

## HTSLIB

#add HTSlib
FetchContent_Declare(
        htslib
        GIT_REPOSITORY https://github.com/samtools/htslib
)
FetchContent_GetProperties(htslib)
if(NOT htslib_POPULATED)
    set(disable_flags --disable-gcs --disable-s3 --disable-plugins --disable-lzma 
        --disable-bz2 --disable-libcurl --without-libdeflate)
    

    FetchContent_Populate(htslib)
    message("HTSLIB_SRC: ${htslib_SOURCE_DIR}")
    message("HTSLIB_BIN: ${htslib_BINARY_DIR}")
    execute_process(COMMAND autoreconf -i
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${htslib_SOURCE_DIR})
    if (result)
        message(FATAL_ERROR "autoreconf step for htslib failed: ${result}")
    endif ()
    execute_process(COMMAND ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR} ${disable_flags}
            WORKING_DIRECTORY ${htslib_SOURCE_DIR})
    if (result)
        message(FATAL_ERROR "configure step for htslib failed: ${result}")
    endif ()
    execute_process(COMMAND make lib-static
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${htslib_SOURCE_DIR})
    if (result)
        message(FATAL_ERROR "make step for htslib failed: ${result}")
    endif ()
    execute_process(COMMAND make install prefix=${CMAKE_CURRENT_BINARY_DIR}
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${htslib_SOURCE_DIR})
    if (result)
        message(FATAL_ERROR "install step for htslib failed: ${result}")    
    endif()
    set(htslib_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(htslib_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libhts.a)
endif()
add_library(htslib STATIC IMPORTED GLOBAL)
set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${htslib_LIB})
target_include_directories(htslib INTERFACE ${htslib_SRC})
