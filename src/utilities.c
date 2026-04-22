/**
 * @file utilities.c
 * @brief Common utility functions using standard C libraries
 * @author An Nguyen
 * @date Last updated: 08/2025
 *
 * Provides efficient implementations of common operations:
 * - Math utilities (05/2022)
 * - String handling (05/2025)
 * - File operations (05/2025)
 * - Memory allocation helpers (08/2025)
 */

#include "define.h"
#include "variables.h"
#include "file.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <signal.h>
#include <ctype.h>  // For isspace function
#include <locale.h>  // For setlocale function
#include "diagnostics.h" 

double Min(double a, double b) {
    return (a < b) ? a : b;
}

double Max(double a, double b) {
    return (a > b) ? a : b;
}

/**
 * @brief Find the maximum value in a double array
 * @param array Pointer to the array
 * @return Maximum value found
 */
double find_max_value(double* array) {
    if (!array) return 0.0;
    double max_val = array[1]; // Start from index 1 as per model convention
    for (int i = 2; i <= M; i++) {
        if (array[i] > max_val) {
            max_val = array[i];
        }
    }
    return max_val;
}

/**
 * @brief Find the minimum value in a double array
 * @param array Pointer to the array  
 * @return Minimum value found
 */
double find_min_value(double* array) {
    if (!array) return 0.0;
    double min_val = array[1]; // Start from index 1 as per model convention
    for (int i = 2; i <= M; i++) {
        if (array[i] < min_val) {
            min_val = array[i];
        }
    }
    return min_val;
}

/**
 * @brief Safely allocates memory and checks for allocation failure.
 * 
 * @param size Size of memory to allocate in bytes
 * @param description Description for error reporting
 * @return Pointer to allocated memory
 */
void* safeMalloc(size_t size, const char* description) {
    void* ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "❌ ERROR: Memory allocation failed for %s (%zu bytes)\n", 
                description, size);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

/**
 * @brief Safely opens a file and checks for errors.
 * 
 * @param filename File to open
 * @param mode File opening mode ("r", "w", etc.)
 * @param description Description for error reporting
 * @return FILE pointer to opened file
 */
FILE* safeFileOpen(const char* filename, const char* mode, const char* description) {
    FILE* file = fopen(filename, mode);
    if (!file) {
        fprintf(stderr, "❌ ERROR: Could not open %s file '%s'\n", 
                description, filename);
        perror("Reason");
    }
    return file;
}

/**
 * @brief Trims whitespace from start and end of a string.
 * 
 * @param str String to trim
 * @return Pointer to the trimmed string (same as input)
 */
char* trimString(char* str) {
    if (!str) return NULL;
    
    // Skip leading whitespace
    char* start = str;
    while (isspace(*start)) start++;
    
    // All spaces?
    if (*start == 0) {
        *str = 0;
        return str;
    }
    
    // Find end of string
    char* end = start + strlen(start) - 1;
    
    // Trim trailing whitespace
    while (end > start && isspace(*end)) end--;
    
    // Terminate string after last non-whitespace character
    *(end + 1) = '\0';
    
    // Move string to beginning if needed
    if (start != str) {
        memmove(str, start, (end - start) + 2);
    }
    
    return str;
}

double linearInterpolate(double x1, double y1, double x2, double y2, double x) {
    if (fabs(x2 - x1) < 1e-10) return y1;  // Avoid division by zero
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

/**
 * @brief Gets current timestamp as a formatted string.
 * 
 * @param buffer Buffer to store timestamp string
 * @param size Size of buffer
 * @return Pointer to buffer for convenience
 */
char* getTimestamp(char* buffer, size_t size) {
    time_t now = time(NULL);
    struct tm* tm_info = localtime(&now);
    strftime(buffer, size, "%Y-%m-%d %H:%M:%S", tm_info);
    return buffer;
}


// =====================================================================
// SIGNAL HANDLING AND CLEANUP
// =====================================================================

/**
 * @brief Perform final cleanup operations
 */
void perform_final_cleanup(void) {
    close_all_model_output_files();
}

/**
 * @brief Signal handler to ensure cleanup on forced stop
 * @param sig Signal number received
 */
void handle_termination_signal(int sig) {
    printf("\n[INFO] Caught termination signal (%d). Cleaning up and exiting...\n", sig);
    perform_final_cleanup();
    fflush(stdout);
    exit(0);
}

/**
 * @brief Register signal handlers for clean termination
 */
void register_signal_handlers(void) {
    signal(SIGINT, handle_termination_signal);   // Ctrl+C
    signal(SIGTERM, handle_termination_signal);  // Termination
    #ifdef _WIN32
    signal(SIGBREAK, handle_termination_signal); // Windows Ctrl+Break
    #endif
}

// =====================================================================
// CONFIGURATION AND SETUP
// =====================================================================

/**
 * @brief Setup locale and register signal handlers
 */
void setup_locale_and_signals(void) {
    // This prevents locale-dependent comma decimal separators from corrupting CSV format
    setlocale(LC_NUMERIC, "C");
    
    // Register signal handlers for clean termination
    register_signal_handlers();
    
    printf("✅ Starting C-GEM Model...\n");
}

/**
 * @brief Initialize simulation variables and counters
 * @param start_time Simulation start time
 * @param steps_done Pointer to steps done counter
 * @param last_printed_day Pointer to last printed day
 * @param progress_pct Pointer to progress percentage
 * @param last_progress_time Pointer to last progress time
 * @param t_start Pointer to simulation start time
 */
void initialize_simulation_variables(clock_t start_time, long* steps_done, double* last_printed_day,
                                    double* progress_pct, clock_t* last_progress_time, long* t_start) {
    // Initialize counters
    *steps_done = 0;
    *progress_pct = 0;
    *last_printed_day = -5.0;
    *last_progress_time = start_time;
    *t_start = 0;
}

// =====================================================================
// PARAMETER FILE PATH RESOLUTION
// =====================================================================

const char *cgem_get_params_path(void) {
    // Note: getenv() returns a pointer that remains valid until the environment
    // is modified. We intentionally do not copy to avoid buffer management.
    const char *env = getenv("CGEM_PARAMS_PATH");
    if (env && env[0] != '\0') {
        return env;
    }
    return "INPUT/params.txt";
}
