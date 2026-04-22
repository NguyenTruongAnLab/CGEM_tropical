/**
 * @file diagnostics.h
 * @brief Consolidated diagnostics, validation, and debugging system for C-GEM model
 * @author Nguyen Truong An (refactored by GitHub Copilot)
 * @date Last updated: 08/2025
 * 
 * Key Features:
 * - Centralized debug logging with configurable levels
 * - Transport physics validation and monitoring
 * - Hydrodynamic-transport coupling diagnostics
 * - Performance monitoring and statistics
 * - End-of-simulation summaries
 * - Milestone tracking and periodic reporting
 */

#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include <stdbool.h>
#include <time.h>

// =====================================================================
// DEBUG SYSTEM - From debug.h/debug.c
// =====================================================================

// Debug level enumeration - improved 4-level system for better control
typedef enum {
    DEBUG_LEVEL_OFF = 0,        // No debug output
    DEBUG_LEVEL_ESSENTIAL = 1,  // Key informational checkpoints and critical errors (< 10 messages per simulation)
    DEBUG_LEVEL_WARNING = 2,    // Essential + important warnings (< 50 messages per simulation)
    DEBUG_LEVEL_DETAIL = 3      // Full diagnostic monitoring (throttled to prevent spam)
} DebugLevel;

// Module identifiers for debug messages
typedef enum {
    DEBUG_MODULE_GENERAL = 0,
    DEBUG_MODULE_INIT = 1,
    DEBUG_MODULE_HYDRO = 2,
    DEBUG_MODULE_TRANSPORT = 3,
    DEBUG_MODULE_BC = 4,     // Boundary Conditions
    DEBUG_MODULE_BIOGEO = 5,
    DEBUG_MODULE_FILE = 6,
    DEBUG_MODULE_COUPLING = 7,  // Hydrodynamic-transport coupling
    DEBUG_MODULE_COUNT
} DebugModule;

// External global variable to control debug level
extern int debug_level;

// Log frequency control (seconds)
#define DEBUG_FREQ_ALWAYS 0       // Print every time
#define DEBUG_FREQ_HOURLY 3600    // Once per hour
#define DEBUG_FREQ_DAILY 86400    // Once per day
#define DEBUG_FREQ_WEEKLY 604800  // Once per week

// =====================================================================
// HYDRODYNAMIC-TRANSPORT COUPLING DIAGNOSTICS 
// =====================================================================

// Diagnostic data structure to track hydrodynamic field updates
typedef struct {
    int last_update_time;           // Last time hydrodynamic fields were updated
    double last_velocity_checksum;  // Checksum of velocity field for change detection
    double last_depth_checksum;     // Checksum of depth field for change detection
    double last_area_checksum;      // Checksum of area field for change detection
    bool velocity_updated;          // Flag indicating if velocity was updated this timestep
    bool depth_updated;             // Flag indicating if depth was updated this timestep
    bool area_updated;              // Flag indicating if area was updated this timestep
    int validation_failures;       // Count of validation failures
} HydroTransportCouplingData;

// =====================================================================
// MASS BUDGET DIAGNOSTICS
// =====================================================================

typedef enum {
    MASS_BUDGET_STAGE_BOUNDARY = 0,   ///< Dirichlet boundary enforcement (Openbound)
    MASS_BUDGET_STAGE_TRIBUTARY = 1,  ///< Lateral source mixing (Fischer 1979)
    MASS_BUDGET_STAGE_ADVECTION = 2,  ///< TVD advection step
    MASS_BUDGET_STAGE_DIFFUSION = 3,  ///< Dispersion / diffusion step
    MASS_BUDGET_STAGE_COUNT
} MassBudgetStage;

// =====================================================================
// BIOGEOCHEMICAL REACTION MASS BUDGET DIAGNOSTICS
// =====================================================================

typedef enum {
    REACTION_PROCESS_PRIMARY_NH4 = 0,   ///< Phytoplankton growth on NH4 uptake
    REACTION_PROCESS_PRIMARY_NO3 = 1,   ///< Phytoplankton growth on NO3 uptake
    REACTION_PROCESS_PHYTO_DEATH = 2,   ///< Phytoplankton mortality / lysis
    REACTION_PROCESS_ADEGRAD = 3,       ///< Oxic mineralisation of organic matter
    REACTION_PROCESS_DENIT = 4,         ///< Heterotrophic denitrification
    REACTION_PROCESS_NITRIF = 5,        ///< Nitrification (NH4 → NO3)
    REACTION_PROCESS_O2_AIR = 6,        ///< Air-water gas exchange for O2
    REACTION_PROCESS_CO2_AIR = 7,       ///< Air-water gas exchange for CO2/DIC
    REACTION_PROCESS_BOUNDING = 8,      ///< Numerical bounding / clipping adjustments
    REACTION_PROCESS_EROSION = 9,       ///< SPM erosion from bed (Winterwerp 2002)
    REACTION_PROCESS_DEPOSITION = 10,   ///< SPM deposition to bed (Winterwerp 2002)
    REACTION_PROCESS_COUNT
} ReactionProcess;

// =====================================================================
// VALIDATION AND MONITORING STRUCTURES
// =====================================================================

// Structure to track validation metrics for end-of-simulation summary
typedef struct {
    // Salinity metrics
    double max_salt_intrusion_km;
    double min_salt_intrusion_km;
    double avg_salt_intrusion_km;
    
    // Velocity metrics  
    double max_velocity_ms;
    double avg_velocity_ms;
    
    // Dispersion metrics
    double max_dispersion_m2s;
    double min_dispersion_m2s;
    double avg_dispersion_m2s;
    
    // Tidal metrics
    double max_tidal_range_m;
    double avg_tidal_range_m;
    
    // Mass conservation
    double max_mass_error;
    double avg_mass_error;
    
    // Sample counts
    int validation_samples;
    int milestone_days_logged;
} ValidationSummary;

/**
 * @brief Initialize the complete diagnostics system
 * @param level Debug level to set (0=off, 1=essential, 2=detailed)
 */
void diagnostics_init(int level);

// -------------------------
// Debug Logging Functions
// -------------------------

/**
 * @brief Log a debug message if current debug level permits
 * @param level Minimum debug level required to show this message
 * @param module Source module identifier
 * @param format Printf-style format string
 * @param ... Additional arguments for format string
 */
void debug_log(DebugLevel level, DebugModule module, const char *format, ...);

/**
 * @brief Log a debug message at specified intervals to reduce spam
 * @param level Minimum debug level required to show this message
 * @param module Source module identifier
 * @param frequency_seconds Minimum seconds between log messages
 * @param key Unique identifier for this message type
 * @param format Printf-style format string
 * @param ... Additional arguments for format string
 */
void debug_log_periodic(DebugLevel level, DebugModule module, 
                        int frequency_seconds, const char *key, 
                        const char *format, ...);

/**
 * @brief Check if debug output should be shown based on current debug level
 * @param level Minimum required debug level for output
 * @return true if debug output should be shown
 */
bool debug_is_enabled(DebugLevel level);

/**
 * @brief Get string name for debug module
 * @param module Module enum value
 * @return const char* Module name string
 */
const char* debug_get_module_name(DebugModule module);

/**
 * @brief Initialize hydro-transport coupling diagnostics
 */
void init_coupling_diagnostics(void);

/**
 * @brief Checks for early warning signs of numerical instability
 * @details Designed to detect issues before they lead to model crashes
 * @param t Current simulation time
 * @return true if simulation should continue, false if critical issue detected
 */
bool check_numerical_stability(int t);

/**
 * @brief Validate prognostic tracer state arrays and hard-abort on invalid values.
 *
 * Enforces the project invariant (Memory Bank): all concentrations must remain
 * >= 0.0 and finite. This function terminates the simulation immediately if
 * any environmentally simulated tracer violates the constraint.
 *
 * @param stage Short context string (e.g., "Transport:TVD", "Biogeo:exit")
 * @param t Current simulation time (seconds)
 */
void diagnostics_validate_tracer_state(const char* stage, int t);

/**
 * @brief Report a fatal tracer invariant violation.
 *
 * In normal runs, this terminates the process. In calibration mode, it latches
 * a fatal request so the current iteration can be penalized without killing the
 * optimizer.
 */
void diagnostics_fatal_tracer_error(const char* stage, int t, int s, int i, double value, const char* reason);

/**
 * @brief In calibration mode, diagnostics can request aborting the current simulation run
 *        (without terminating the whole process).
 *
 * Normal runs still hard-exit on tracer invariant violations.
 */
bool diagnostics_fatal_requested(void);

/**
 * @brief Clear any latched fatal request (typically at the start of each calibration iteration).
 */
void diagnostics_clear_fatal_request(void);

/**
 * @brief Check if this time is a milestone day requiring debug output
 * @param t Current time in seconds
 * @return true if this is a milestone day for reporting
 */
bool is_milestone_day(int t);

/**
 * @brief Find maximum value in an array
 * @param array Array of values
 * @return Maximum value
 */
double find_max_value(double *array);

/**
 * @brief Find minimum value in an array
 * @param array Array of values  
 * @return Minimum value
 */
double find_min_value(double *array);

// =====================================================================
// OUTPUT CONTROL (merged from output_control.h)
// =====================================================================

extern int global_output_enabled;

void init_output_control(int enabled);
int printf_controlled(const char *format, ...);
void set_output_enabled(int enabled);
int is_output_enabled(void);

// Smart printf macros - suppressed during calibration
#define PRINTF_INIT(...) printf_controlled(__VA_ARGS__)
#define PRINTF_DATA(...) printf_controlled(__VA_ARGS__)
#define PRINTF_CALIB(...) printf(__VA_ARGS__)  // Always visible
#define PRINTF_ERROR(...) printf(__VA_ARGS__)  // Always visible

#endif /* DIAGNOSTICS_H */