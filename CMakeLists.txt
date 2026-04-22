# ==============================================================================
# C-GEM (Carbon-Generic Estuarine Model) - Cross-platform CMake Build
# ==============================================================================
# Build Options:
#   cmake -DCGEM_ENABLE_CALIBRATION=OFF ..  (default, no NLopt needed)
#   cmake -DCGEM_ENABLE_CALIBRATION=ON ..   (requires NLopt library)
#
# Compatible with: Windows (MinGW/MSVC), Linux (GCC), macOS (Clang)
#
# Usage for Linux CI (no calibration):
#   mkdir build && cd build
#   cmake -DCMAKE_BUILD_TYPE=Release -DCGEM_ENABLE_CALIBRATION=OFF ..
#   cmake --build .
#
# Scientific guardrails preserved:
# - Hybrid forcing (hourly tide + daily NET discharge)
# - Fischer (1979) mass-source tributaries
# - Configuration-driven design (INPUT/*.txt)
# ==============================================================================

cmake_minimum_required(VERSION 3.15)
project(C_GEM_Tropical 
    VERSION 1.0.0
    DESCRIPTION "1D reactive-transport model for data-sparse tropical estuaries"
    LANGUAGES C
)

# C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ==============================================================================
# Build Options
# ==============================================================================

option(CGEM_ENABLE_CALIBRATION "Enable calibration module with NLopt" OFF)
option(CGEM_ENABLE_WARNINGS "Enable compiler warnings (helps remove dead/unused code)" ON)
option(CGEM_BUILD_TESTS "Build unit tests when test sources are present" ON)

# ==============================================================================
# Output directories
# ==============================================================================

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/Release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_SOURCE_DIR}/bin/Debug)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/bin/Release)

# ==============================================================================
# Core source files (always compiled)
# ==============================================================================

set(CORE_SOURCES
    src/bcforcing.c
    src/biogeo.c
    src/biogeout.c
    src/calibration_helpers.c
    src/diagnostics.c
    src/hydrodynamics.c
    src/init.c
    src/io_input.c
    src/io_output.c
    src/main.c
    src/transport.c
    src/utilities.c
    src/variables.c
)

# ==============================================================================
# Conditional calibration sources
# ==============================================================================

if(CGEM_ENABLE_CALIBRATION)
    message(STATUS "Calibration module: ENABLED (requires NLopt)")
    list(APPEND CORE_SOURCES src/calibration.c)
    
    # NLopt dependency
    # - On MSVC/Windows: the repo currently ships a runtime DLL but not an MSVC .lib import library.
    #   Use a small runtime shim (src/nlopt_winshim.c) that LoadLibrary()s libnlopt.dll.
    # - On MinGW/Unix: link normally via the provided import/static library.

    find_path(NLOPT_INCLUDE_DIR nlopt.h
        PATHS ${CMAKE_SOURCE_DIR}/external/nlopt/src/api
        NO_DEFAULT_PATH)

    if(NOT NLOPT_INCLUDE_DIR)
        message(FATAL_ERROR "NLopt headers not found. Expected external/nlopt/src/api/nlopt.h")
    endif()

    if(MSVC)
        message(STATUS "  NLopt: MSVC build will use runtime DLL shim")
        list(APPEND CORE_SOURCES src/nlopt_winshim.c)
        set(NLOPT_LIBRARY "")
    else()
        find_library(NLOPT_LIBRARY
            NAMES nlopt libnlopt
            PATHS ${CMAKE_SOURCE_DIR}/external/nlopt/build
            NO_DEFAULT_PATH)

        if(NOT NLOPT_LIBRARY)
            message(FATAL_ERROR "NLopt library not found. Disable CGEM_ENABLE_CALIBRATION or provide NLopt build artifacts.")
        endif()
        message(STATUS "  NLopt library: ${NLOPT_LIBRARY}")
    endif()

    # Runtime DLL (Windows): ensure the executable can find libnlopt.dll at runtime.
    if(WIN32)
        find_file(NLOPT_DLL
            NAMES libnlopt.dll
            PATHS ${CMAKE_SOURCE_DIR}/external/nlopt/build
            NO_DEFAULT_PATH)
        if(NLOPT_DLL)
            message(STATUS "  NLopt DLL: ${NLOPT_DLL}")
        endif()
    endif()

    message(STATUS "  NLopt include: ${NLOPT_INCLUDE_DIR}")
else()
    message(STATUS "Calibration module: DISABLED (using stub)")
    list(APPEND CORE_SOURCES src/calibration_stub.c)
endif()

# ==============================================================================
# Executable target
# ==============================================================================

add_executable(C_GEM_Daily ${CORE_SOURCES})

# ==============================================================================
# Compile definitions
# ==============================================================================

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(C_GEM_Daily PRIVATE CGEM_ENABLE_CALIBRATION=1)
else()
    target_compile_definitions(C_GEM_Daily PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

# ==============================================================================
# Include directories
# ==============================================================================

target_include_directories(C_GEM_Daily PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    # Treat NLopt headers as system/external when supported to avoid warning noise.
    target_include_directories(C_GEM_Daily SYSTEM PRIVATE ${NLOPT_INCLUDE_DIR})
endif()

# ==============================================================================
# Compiler flags (optimization and warnings)
# ==============================================================================

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    # GCC or Clang (Linux/MinGW)
    target_compile_options(C_GEM_Daily PRIVATE
        -O3
        -march=native
        -funroll-loops
        -fomit-frame-pointer
    )

    if(CGEM_ENABLE_WARNINGS)
        target_compile_options(C_GEM_Daily PRIVATE
            -Wall
            -Wextra
            -Wpedantic
        )
    endif()
    
    # Enable LTO if supported
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.9)
        set_property(TARGET C_GEM_Daily PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
    
elseif(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    # Microsoft Visual C++
    target_compile_options(C_GEM_Daily PRIVATE
        /O2
    )

    if(CGEM_ENABLE_WARNINGS)
        # /W4 surfaces unused/static issues. Do NOT use /WX (warnings-as-errors) for legacy code.
        target_compile_options(C_GEM_Daily PRIVATE /W4)
        # If available, reduce warnings from external/system headers.
        target_compile_options(C_GEM_Daily PRIVATE /external:anglebrackets /external:W0)
    else()
        target_compile_options(C_GEM_Daily PRIVATE /W0)
    endif()

    # Suppress "unsafe" function warnings (sscanf, strtok, etc.)
    target_compile_definitions(C_GEM_Daily PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

# ==============================================================================
# Linker flags and libraries
# ==============================================================================

# Math library (required on Linux, automatic on Windows)
if(UNIX)
    target_link_libraries(C_GEM_Daily m)
endif()

# NLopt library (if calibration enabled)
if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(C_GEM_Daily ${NLOPT_LIBRARY})
    endif()

    # Copy NLopt DLL next to executable on Windows so local runs work without PATH tweaks.
    if(WIN32 AND NLOPT_DLL)
        add_custom_command(TARGET C_GEM_Daily POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${NLOPT_DLL}"
                "$<TARGET_FILE_DIR:C_GEM_Daily>/libnlopt.dll"
            COMMENT "Copying NLopt runtime DLL next to C_GEM_Daily"
        )
    endif()
endif()

# ==============================================================================
# Testing (smoke test with environment variable)
# ==============================================================================

enable_testing()

add_test(NAME smoke_test_1day
    COMMAND ${CMAKE_COMMAND} -E env CGEM_FORCE_NO_CALIBRATION=1 CGEM_SKIP_VALIDATION=1 CGEM_MASS_DIAG_DAYS=1 CGEM_WARMUP_DAYS=0 CGEM_MAXT_DAYS=1
        $<TARGET_FILE:C_GEM_Daily>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(smoke_test_1day PROPERTIES
    TIMEOUT 180
    LABELS "smoke"
)

if(CGEM_BUILD_TESTS AND EXISTS ${CMAKE_SOURCE_DIR}/tests)

# ------------------------------------------------------------------------------
# Unit tests (no INPUT/ dependencies)
# ------------------------------------------------------------------------------

set(CORE_SOURCES_NO_MAIN ${CORE_SOURCES})
list(REMOVE_ITEM CORE_SOURCES_NO_MAIN src/main.c)

add_executable(test_transport_mass_conservation
    tests/test_transport_mass_conservation.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_transport_mass_conservation PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_transport_mass_conservation PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_transport_mass_conservation PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_transport_mass_conservation PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_transport_mass_conservation m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_transport_mass_conservation ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_transport_mass_conservation
    COMMAND $<TARGET_FILE:test_transport_mass_conservation>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_transport_mass_conservation PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)

# ------------------------------------------------------------------------------
# Additional unit tests (transport boundaries + tributary mass source)
# ------------------------------------------------------------------------------

add_executable(test_transport_openbound
    tests/test_transport_openbound.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_transport_openbound PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_transport_openbound PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_transport_openbound PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_transport_openbound PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_transport_openbound m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_transport_openbound ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_transport_openbound
    COMMAND $<TARGET_FILE:test_transport_openbound>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_transport_openbound PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_transport_dispcoef_discharge_parity
    tests/test_transport_dispcoef_discharge_parity.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_transport_dispcoef_discharge_parity PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_transport_dispcoef_discharge_parity PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_transport_dispcoef_discharge_parity PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_transport_dispcoef_discharge_parity PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_transport_dispcoef_discharge_parity m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_transport_dispcoef_discharge_parity ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_transport_dispcoef_discharge_parity
    COMMAND $<TARGET_FILE:test_transport_dispcoef_discharge_parity>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_transport_dispcoef_discharge_parity PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_transport_tributary_fischer
    tests/test_transport_tributary_fischer.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_transport_tributary_fischer PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_transport_tributary_fischer PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_transport_tributary_fischer PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_transport_tributary_fischer PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_transport_tributary_fischer m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_transport_tributary_fischer ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_transport_tributary_fischer
    COMMAND $<TARGET_FILE:test_transport_tributary_fischer>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_transport_tributary_fischer PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_transport_boundary_cv_length
    tests/test_transport_boundary_cv_length.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_transport_boundary_cv_length PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_transport_boundary_cv_length PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_transport_boundary_cv_length PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_transport_boundary_cv_length PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_transport_boundary_cv_length m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_transport_boundary_cv_length ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_transport_boundary_cv_length
    COMMAND $<TARGET_FILE:test_transport_boundary_cv_length>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_transport_boundary_cv_length PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_forcing_continuity_at_warmup
    tests/test_forcing_continuity_at_warmup.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_forcing_continuity_at_warmup PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_forcing_continuity_at_warmup PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_forcing_continuity_at_warmup PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_forcing_continuity_at_warmup PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_forcing_continuity_at_warmup m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_forcing_continuity_at_warmup ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_forcing_continuity_at_warmup
    COMMAND $<TARGET_FILE:test_forcing_continuity_at_warmup>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_forcing_continuity_at_warmup PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_light_temperature_sources_consistent_with_interpolator
    tests/test_light_temperature_sources_consistent_with_interpolator.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_light_temperature_sources_consistent_with_interpolator PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_light_temperature_sources_consistent_with_interpolator PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_light_temperature_sources_consistent_with_interpolator PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_light_temperature_sources_consistent_with_interpolator PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_light_temperature_sources_consistent_with_interpolator m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_light_temperature_sources_consistent_with_interpolator ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_light_temperature_sources_consistent_with_interpolator
    COMMAND $<TARGET_FILE:test_light_temperature_sources_consistent_with_interpolator>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_light_temperature_sources_consistent_with_interpolator PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_params_override_biogeo_rates
    tests/test_params_override_biogeo_rates.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_params_override_biogeo_rates PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_params_override_biogeo_rates PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_params_override_biogeo_rates PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_params_override_biogeo_rates PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_params_override_biogeo_rates m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_params_override_biogeo_rates ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_params_override_biogeo_rates
    COMMAND $<TARGET_FILE:test_params_override_biogeo_rates>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_params_override_biogeo_rates PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_missing_param_behavior_is_documented
    tests/test_missing_param_behavior_is_documented.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_missing_param_behavior_is_documented PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_missing_param_behavior_is_documented PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_missing_param_behavior_is_documented PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_missing_param_behavior_is_documented PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_missing_param_behavior_is_documented m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_missing_param_behavior_is_documented ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_missing_param_behavior_is_documented
    COMMAND $<TARGET_FILE:test_missing_param_behavior_is_documented>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_missing_param_behavior_is_documented PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_calibration_target_summary
    tests/test_calibration_target_summary.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_calibration_target_summary PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_calibration_target_summary PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_calibration_target_summary PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_calibration_target_summary PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_calibration_target_summary m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_calibration_target_summary ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_calibration_target_summary
    COMMAND $<TARGET_FILE:test_calibration_target_summary>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_calibration_target_summary PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_calibration_timebase_alignment
    tests/test_calibration_timebase_alignment.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_calibration_timebase_alignment PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_calibration_timebase_alignment PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_calibration_timebase_alignment PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_calibration_timebase_alignment PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_calibration_timebase_alignment m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_calibration_timebase_alignment ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_calibration_timebase_alignment
    COMMAND $<TARGET_FILE:test_calibration_timebase_alignment>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_calibration_timebase_alignment PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_calibration_score_breakdown
    tests/test_calibration_score_breakdown.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_calibration_score_breakdown PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_calibration_score_breakdown PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_calibration_score_breakdown PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_calibration_score_breakdown PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_calibration_score_breakdown m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_calibration_score_breakdown ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_calibration_score_breakdown
    COMMAND $<TARGET_FILE:test_calibration_score_breakdown>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_calibration_score_breakdown PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_biogeo_co2_exchange_sign
    tests/test_biogeo_co2_exchange_sign.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_biogeo_co2_exchange_sign PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_biogeo_co2_exchange_sign PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_biogeo_co2_exchange_sign PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_biogeo_co2_exchange_sign PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_biogeo_co2_exchange_sign m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_biogeo_co2_exchange_sign ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_biogeo_co2_exchange_sign
    COMMAND $<TARGET_FILE:test_biogeo_co2_exchange_sign>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_biogeo_co2_exchange_sign PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_carbonate_diagnostics_transport_only_invariant
    tests/test_carbonate_diagnostics_transport_only_invariant.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_carbonate_diagnostics_transport_only_invariant PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_carbonate_diagnostics_transport_only_invariant PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_carbonate_diagnostics_transport_only_invariant PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_carbonate_diagnostics_transport_only_invariant PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_carbonate_diagnostics_transport_only_invariant m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_carbonate_diagnostics_transport_only_invariant ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_carbonate_diagnostics_transport_only_invariant
    COMMAND $<TARGET_FILE:test_carbonate_diagnostics_transport_only_invariant>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_carbonate_diagnostics_transport_only_invariant PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_hydrodynamics_tributary_source_scaling
    tests/test_hydrodynamics_tributary_source_scaling.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_hydrodynamics_tributary_source_scaling PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_hydrodynamics_tributary_source_scaling PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_hydrodynamics_tributary_source_scaling PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_hydrodynamics_tributary_source_scaling PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_hydrodynamics_tributary_source_scaling m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_hydrodynamics_tributary_source_scaling ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_hydrodynamics_tributary_source_scaling
    COMMAND $<TARGET_FILE:test_hydrodynamics_tributary_source_scaling>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_hydrodynamics_tributary_source_scaling PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_biogeo_limiting_avoids_clamp
    tests/test_biogeo_limiting_avoids_clamp.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_biogeo_limiting_avoids_clamp PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_biogeo_limiting_avoids_clamp PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_biogeo_limiting_avoids_clamp PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_biogeo_limiting_avoids_clamp PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_biogeo_limiting_avoids_clamp m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_biogeo_limiting_avoids_clamp ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_biogeo_limiting_avoids_clamp
    COMMAND $<TARGET_FILE:test_biogeo_limiting_avoids_clamp>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_biogeo_limiting_avoids_clamp PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


add_executable(test_biogeo_o2_toc_phy_budget_parity
    tests/test_biogeo_o2_toc_phy_budget_parity.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_biogeo_o2_toc_phy_budget_parity PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_biogeo_o2_toc_phy_budget_parity PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_biogeo_o2_toc_phy_budget_parity PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_biogeo_o2_toc_phy_budget_parity PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_biogeo_o2_toc_phy_budget_parity m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_biogeo_o2_toc_phy_budget_parity ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_biogeo_o2_toc_phy_budget_parity
    COMMAND $<TARGET_FILE:test_biogeo_o2_toc_phy_budget_parity>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_biogeo_o2_toc_phy_budget_parity PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)


# This test source is not always present in downstream workspaces (e.g. when
# scientific branches prune experimental tests). Guard it so Release builds do
# not fail when the file is missing.
if(EXISTS ${CMAKE_SOURCE_DIR}/tests/test_biogeo_co2_limitation_prevents_runaway.c)
    add_executable(test_biogeo_co2_limitation_prevents_runaway
        tests/test_biogeo_co2_limitation_prevents_runaway.c
        ${CORE_SOURCES_NO_MAIN}
    )

    target_include_directories(test_biogeo_co2_limitation_prevents_runaway PRIVATE ${CMAKE_SOURCE_DIR}/src)

    if(CGEM_ENABLE_CALIBRATION)
        target_compile_definitions(test_biogeo_co2_limitation_prevents_runaway PRIVATE CGEM_ENABLE_CALIBRATION=1)
        target_include_directories(test_biogeo_co2_limitation_prevents_runaway PRIVATE ${NLOPT_INCLUDE_DIR})
    else()
        target_compile_definitions(test_biogeo_co2_limitation_prevents_runaway PRIVATE CGEM_ENABLE_CALIBRATION=0)
    endif()

    if(UNIX)
        target_link_libraries(test_biogeo_co2_limitation_prevents_runaway m)
    endif()

    if(CGEM_ENABLE_CALIBRATION)
        if(NLOPT_LIBRARY)
            target_link_libraries(test_biogeo_co2_limitation_prevents_runaway ${NLOPT_LIBRARY})
        endif()
    endif()

    add_test(NAME test_biogeo_co2_limitation_prevents_runaway
        COMMAND $<TARGET_FILE:test_biogeo_co2_limitation_prevents_runaway>
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

    set_tests_properties(test_biogeo_co2_limitation_prevents_runaway PROPERTIES
        TIMEOUT 30
        LABELS "unit"
    )
endif()


# ------------------------------------------------------------------------------
# Death tests (expected abort on invariant violation)
# ------------------------------------------------------------------------------

add_executable(test_abort_on_negative_concentration
    tests/test_abort_on_negative_concentration.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_abort_on_negative_concentration PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_abort_on_negative_concentration PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_abort_on_negative_concentration PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_abort_on_negative_concentration PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_abort_on_negative_concentration m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_abort_on_negative_concentration ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_abort_on_negative_concentration
    COMMAND $<TARGET_FILE:test_abort_on_negative_concentration>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_abort_on_negative_concentration PROPERTIES
    TIMEOUT 30
    LABELS "unit"
    WILL_FAIL TRUE
    PASS_REGULAR_EXPRESSION "FATAL TRACER ERROR:.*reason=negative"
    FAIL_REGULAR_EXPRESSION "FAIL: Expected abort on negative concentration"
)


add_executable(test_abort_on_nan_inf
    tests/test_abort_on_nan_inf.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_abort_on_nan_inf PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_abort_on_nan_inf PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_abort_on_nan_inf PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_abort_on_nan_inf PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_abort_on_nan_inf m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_abort_on_nan_inf ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_abort_on_nan_inf
    COMMAND $<TARGET_FILE:test_abort_on_nan_inf>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_abort_on_nan_inf PROPERTIES
    TIMEOUT 30
    LABELS "unit"
    WILL_FAIL TRUE
    PASS_REGULAR_EXPRESSION "FATAL TRACER ERROR:.*reason=nonfinite"
    FAIL_REGULAR_EXPRESSION "FAIL: Expected abort on NaN concentration"
)


add_executable(test_carbonate_nutrient_alkalinity
    tests/test_carbonate_nutrient_alkalinity.c
    ${CORE_SOURCES_NO_MAIN}
)

target_include_directories(test_carbonate_nutrient_alkalinity PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(CGEM_ENABLE_CALIBRATION)
    target_compile_definitions(test_carbonate_nutrient_alkalinity PRIVATE CGEM_ENABLE_CALIBRATION=1)
    target_include_directories(test_carbonate_nutrient_alkalinity PRIVATE ${NLOPT_INCLUDE_DIR})
else()
    target_compile_definitions(test_carbonate_nutrient_alkalinity PRIVATE CGEM_ENABLE_CALIBRATION=0)
endif()

if(UNIX)
    target_link_libraries(test_carbonate_nutrient_alkalinity m)
endif()

if(CGEM_ENABLE_CALIBRATION)
    if(NLOPT_LIBRARY)
        target_link_libraries(test_carbonate_nutrient_alkalinity ${NLOPT_LIBRARY})
    endif()
endif()

add_test(NAME test_carbonate_nutrient_alkalinity
    COMMAND $<TARGET_FILE:test_carbonate_nutrient_alkalinity>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(test_carbonate_nutrient_alkalinity PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)

elseif(CGEM_BUILD_TESTS)
    message(WARNING "CGEM_BUILD_TESTS=ON but ${CMAKE_SOURCE_DIR}/tests was not found; only smoke_test_1day will be available.")
endif()

# ==============================================================================
# Installation (optional)
# ==============================================================================

install(TARGETS C_GEM_Daily
    RUNTIME DESTINATION bin
)

# ==============================================================================
# Build summary
# ==============================================================================

message(STATUS "")
message(STATUS "========================================")
message(STATUS "C-GEM Build Configuration")
message(STATUS "========================================")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
message(STATUS "  Calibration: ${CGEM_ENABLE_CALIBRATION}")
if(CGEM_BUILD_TESTS AND EXISTS ${CMAKE_SOURCE_DIR}/tests)
    message(STATUS "  Unit tests: ENABLED")
elseif(CGEM_BUILD_TESTS)
    message(STATUS "  Unit tests: SKIPPED (tests/ directory not present)")
else()
    message(STATUS "  Unit tests: DISABLED")
endif()
message(STATUS "  Output dir: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "========================================")
message(STATUS "")
