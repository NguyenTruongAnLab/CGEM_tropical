/**
 * @file utilities.h
 * @brief Consolidated utilities: core utilities and performance cache
 * @author Nguyen Truong An
 * @date February 2025
 * 
 * This header consolidates all utility functions from:
 * - Core utilities (math, string, file operations)
 */

#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdbool.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>  // For clock_t type

// =====================================================================
// CORE UTILITIES
// =====================================================================

double Min(double a, double b);
double Max(double a, double b);
double find_max_value(double* array);
double find_min_value(double* array);
void* safeMalloc(size_t size, const char* description);
FILE* safeFileOpen(const char* filename, const char* mode, const char* description);
char* trimString(char* str);
double linearInterpolate(double x1, double y1, double x2, double y2, double x);
char* getTimestamp(char* buffer, size_t size);

// =====================================================================
// SIGNAL HANDLING AND CLEANUP
// =====================================================================

void handle_termination_signal(int sig);
void register_signal_handlers(void);
void perform_final_cleanup(void);

// =====================================================================
// CONFIGURATION AND SETUP
// =====================================================================

void setup_locale_and_signals(void);
void initialize_simulation_variables(clock_t start_time, long* steps_done, double* last_printed_day,
                                    double* progress_pct, clock_t* last_progress_time, long* t_start);

// =====================================================================
// PARAMETER FILE PATH RESOLUTION
// =====================================================================

// Return the path to the active params file.
// - Default: "INPUT/params.txt"
// - Override: environment variable CGEM_PARAMS_PATH
// The returned pointer remains valid for the lifetime of the process.
const char *cgem_get_params_path(void);

#endif // UTILITIES_H