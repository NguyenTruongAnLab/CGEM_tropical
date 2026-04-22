#ifndef FILE_H
#define FILE_H

#include <stdio.h>

// -----------------------------------------------------------------------------
// Output interfaces
// -----------------------------------------------------------------------------

void write_field(int t, int cell_count, double *data, const char *var_name);
int Hydwrite(int t);
void SaveAllWaterQuality(int t);
void Rates(double *co, const char *s, int t);
void Fluxwrite(int species, int t);

void close_all_model_output_files(void);
void reset_output_file_registry(void);

void init_simulation_date(void);
const char *get_simulation_start_date(void);

void ensure_directory_for_path(const char *filepath);

void run_postprocessing_export(void);

// -----------------------------------------------------------------------------
// Input interfaces
// -----------------------------------------------------------------------------

int read_calibration_flags(const char *filename,
						   int *calibration_mode_out,
						   int *debug_level_out);

void read_parameters(const char *filename);
void readBoundaryData(const char *filename);
void readGlobalConfigSettings(const char *filename);
void readTributaryData(const char *filename);
void readRiverbedProfile(double *riverbed_depth, const char *filename);

#endif // FILE_H
