#include "file.h"
#include "diagnostics.h"
#include "define.h"
#include "variables.h"
#include "utilities.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <ctype.h>
#include <stdint.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#include <io.h>
#else
#include <unistd.h>
#include <sys/wait.h>
#endif

#define MAX_OUTPUT_STREAMS 256
#define CGEMBIN_HEADER_SIZE (8 + (3 * sizeof(uint32_t)) + (2 * sizeof(double)))

typedef struct {
    char filename[256];
    FILE *handle;
    bool is_binary;
    bool header_written;
    int write_counter;
    void *io_buffer; // Optional heap buffer used by setvbuf(); freed on close
} OutputStream;

static OutputStream output_streams[MAX_OUTPUT_STREAMS];
static int output_stream_count = 0;
static bool output_directories_initialized = false;
static bool geometry_snapshot_written = false;
static char sim_start_date[32] = "2000-01-01";  // Overridden by SIM_START_DATE in params.txt
static bool sim_date_initialized = false;

// Environment-tunable I/O settings (initialized once)
static int env_buffer_size = -1;         // -1 means not initialized
static int env_flush_cadence_binary = 10;
static int env_flush_cadence_csv = 20;

static void init_io_env_settings(void) {
    if (env_buffer_size >= 0) {
        return;  // Already initialized
    }
    
    // Default: use system buffering (no setvbuf call)
    env_buffer_size = 0;
    
    const char *buffer_env = getenv("CGEM_OUTPUT_BUFFER_KB");
    if (buffer_env && buffer_env[0] != '\0') {
        int kb = atoi(buffer_env);
        if (kb > 0 && kb <= 1024) {
            env_buffer_size = kb * 1024;
        }
    }
    
    const char *flush_binary_env = getenv("CGEM_OUTPUT_FLUSH_BINARY");
    if (flush_binary_env && flush_binary_env[0] != '\0') {
        int cadence = atoi(flush_binary_env);
        if (cadence > 0) {
            env_flush_cadence_binary = cadence;
        }
    }
    
    const char *flush_csv_env = getenv("CGEM_OUTPUT_FLUSH_CSV");
    if (flush_csv_env && flush_csv_env[0] != '\0') {
        int cadence = atoi(flush_csv_env);
        if (cadence > 0) {
            env_flush_cadence_csv = cadence;
        }
    }
}

static void ensure_output_directories(void);
static void ensure_directory_for_path_internal(const char *filepath);
static OutputStream *acquire_output_stream(const char *path, bool is_binary);
static void write_binary_header(FILE *fp, uint32_t cell_count, uint32_t values_per_record);
static void write_csv_header(FILE *fp, int cell_count);
static void write_field_internal(int t, int cell_count, const double *data, const char *var_name, const char *csv_format);
static void format_simulation_date(int t, char *date_str);
static void write_geometry_snapshot(void);
static void build_output_path(char *dest, size_t len, const char *var_name, bool is_binary);
static void flush_stream(OutputStream *stream, bool force);
static void normalise_var_name(const char *input, char *output, size_t len);

static int is_centers_only_output_enabled(void);

#ifdef _WIN32
static int create_directory_if_missing(const char *path) {
    if (_mkdir(path) == 0) {
        return 1;
    }
    if (errno == EEXIST) {
        struct _stat st;
        if (_stat(path, &st) == 0 && (st.st_mode & _S_IFDIR)) {
            return 1;
        }
    }
    return 0;
}
#else
static int create_directory_if_missing(const char *path) {
    if (mkdir(path, 0775) == 0) {
        return 1;
    }
    if (errno == EEXIST) {
        struct stat st;
        if (stat(path, &st) == 0 && S_ISDIR(st.st_mode)) {
            return 1;
        }
    }
    return 0;
}
#endif

static void ensure_output_directories(void) {
    if (output_directories_initialized) {
        return;
    }

    (void)create_directory_if_missing("OUT");

    output_directories_initialized = true;
}

void ensure_directory_for_path(const char *filepath) {
    if (!filepath || *filepath == '\0') {
        return;
    }

    char buffer[512];
    strncpy(buffer, filepath, sizeof(buffer) - 1);
    buffer[sizeof(buffer) - 1] = '\0';

    size_t len = strlen(buffer);
    for (size_t i = 0; i < len; ++i) {
        if (buffer[i] == '\\') {
            buffer[i] = '/';
        }
    }

    char *last_sep = strrchr(buffer, '/');
    if (!last_sep) {
        return;
    }
    *last_sep = '\0';

    char partial[512];
    partial[0] = '\0';

    char *token = strtok(buffer, "/");
    while (token) {
        if (partial[0] != '\0') {
            strncat(partial, "/", sizeof(partial) - strlen(partial) - 1);
        }
        strncat(partial, token, sizeof(partial) - strlen(partial) - 1);
        create_directory_if_missing(partial);
        token = strtok(NULL, "/");
    }
}

static void ensure_directory_for_path_internal(const char *filepath) {
    ensure_directory_for_path(filepath);
}

static void normalise_var_name(const char *input, char *output, size_t len) {
    if (!input || *input == '\0') {
        snprintf(output, len, "unknown");
        return;
    }

    const char *start = input;
    if (strncmp(start, "OUT/", 4) == 0) {
        start += 4;
    }
    if (strncmp(start, "OUT\\", 4) == 0) {
        start += 4;
    }

    snprintf(output, len, "%s", start);
    output[len - 1] = '\0';

    for (char *p = output; *p; ++p) {
        if (*p == '\\' || *p == '/') {
            *p = '_';
        } else if (!(isalnum((unsigned char)*p) || *p == '_' || *p == '-')) {
            *p = '_';
        }
    }

    while (output[0] == '_') {
        memmove(output, output + 1, strlen(output));
    }

    char *dot = strrchr(output, '.');
    if (dot) {
        *dot = '\0';
    }

    if (output[0] == '\0') {
        snprintf(output, len, "unknown");
    }
}

static void build_output_path(char *dest, size_t len, const char *var_name, bool is_binary) {
    char name_buffer[256];
    memset(name_buffer, 0, sizeof(name_buffer));
    normalise_var_name(var_name, name_buffer, sizeof(name_buffer));
    const char *ext = is_binary ? "bin" : "csv";
    snprintf(dest, len, "OUT/%s.%s", name_buffer, ext);
    dest[len - 1] = '\0';
}

static OutputStream *acquire_output_stream(const char *path, bool is_binary) {


    for (int i = 0; i < output_stream_count; ++i) {
        if (strcmp(output_streams[i].filename, path) == 0) {
            return &output_streams[i];
        }
    }

    if (output_stream_count >= MAX_OUTPUT_STREAMS) {
        fprintf(stderr, "⚠️ Warning: Maximum number of output streams reached (%d). Skipping %s\n", MAX_OUTPUT_STREAMS, path);
        return NULL;
    }

    ensure_output_directories();
    ensure_directory_for_path_internal(path);

    const char *mode = is_binary ? "wb" : "w";
    FILE *fp = fopen(path, mode);
    if (!fp) {
        fprintf(stderr, "❌ Error: Unable to open %s for writing (%s)\n", path, strerror(errno));
        return NULL;
    }
    
    OutputStream *stream = &output_streams[output_stream_count++];
    strncpy(stream->filename, path, sizeof(stream->filename) - 1);
    stream->filename[sizeof(stream->filename) - 1] = '\0';
    stream->handle = fp;
    stream->is_binary = is_binary;
    stream->header_written = false;
    stream->write_counter = 0;
    stream->io_buffer = NULL;

    // Apply environment-tunable buffering (heap buffer must persist for stream lifetime).
    init_io_env_settings();
    if (env_buffer_size > 0) {
        void *buffer = malloc((size_t)env_buffer_size);
        if (buffer) {
            if (setvbuf(fp, (char *)buffer, _IOFBF, (size_t)env_buffer_size) == 0) {
                stream->io_buffer = buffer;
            } else {
                free(buffer);
            }
        }
    }
    return stream;
}

static void write_binary_header(FILE *fp, uint32_t cell_count, uint32_t values_per_record) {
    const char magic[8] = {'C','G','E','M','B','I','N','\0'};
    const uint32_t version = 2U;
    const double dt_seconds = (double)TS * (double)DELTI;
    const double dx_meters = (double)DELXI;

    fwrite(magic, sizeof(char), sizeof(magic), fp);
    fwrite(&version, sizeof(uint32_t), 1, fp);
    fwrite(&cell_count, sizeof(uint32_t), 1, fp);
    fwrite(&values_per_record, sizeof(uint32_t), 1, fp);
    fwrite(&dt_seconds, sizeof(double), 1, fp);
    fwrite(&dx_meters, sizeof(double), 1, fp);
}

static void write_csv_header(FILE *fp, int cell_count) {
    fprintf(fp, "Date");
    for (int i = 1; i <= cell_count; ++i) {
        fprintf(fp, ",Cell_%d", i);
    }
    fputc('\n', fp);
}

static void flush_stream(OutputStream *stream, bool force) {
    if (!stream || !stream->handle) {
        return;
    }
    if (force) {
        fflush(stream->handle);
        return;
    }
    
    init_io_env_settings();
    
    if (stream->is_binary) {
        if (stream->write_counter % env_flush_cadence_binary == 0) {
            fflush(stream->handle);
        }
    } else {
        if (stream->write_counter % env_flush_cadence_csv == 0) {
            fflush(stream->handle);
        }
    }
}

static void init_simulation_date_from_params(void) {
    if (sim_date_initialized) {
        return;
    }

    // Optional override: allow shifting the calendar without editing INPUT files.
    // Expected format: YYYY-MM-DD
    const char *env_start_date = getenv("CGEM_SIM_START_DATE");
    if (env_start_date && env_start_date[0] != '\0') {
        int y = 0, m = 0, d = 0;
        if (sscanf(env_start_date, "%d-%d-%d", &y, &m, &d) == 3 && y >= 1900 && m >= 1 && m <= 12 && d >= 1 && d <= 31) {
            size_t len = strlen(env_start_date);
            if (len < sizeof(sim_start_date)) {
                strncpy(sim_start_date, env_start_date, sizeof(sim_start_date) - 1);
                sim_start_date[sizeof(sim_start_date) - 1] = '\0';
            }
            sim_date_initialized = true;
            return;
        } else {
            printf("⚠️ Ignoring CGEM_SIM_START_DATE='%s' (expected YYYY-MM-DD)\n", env_start_date);
        }
    }

    FILE *fp = fopen(cgem_get_params_path(), "r");
    if (!fp) {
        sim_date_initialized = true;
        return;
    }

    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "SIM_START_DATE") && strchr(line, '=')) {
            char *date_start = strchr(line, '"');
            if (date_start) {
                date_start++;
                char *date_end = strchr(date_start, '"');
                if (date_end) {
                    size_t len = (size_t)(date_end - date_start);
                    if (len < sizeof(sim_start_date)) {
                        strncpy(sim_start_date, date_start, len);
                        sim_start_date[len] = '\0';
                    }
                }
            }
            break;
        }
    }
    fclose(fp);
    sim_date_initialized = true;
}

void init_simulation_date(void) {
    init_simulation_date_from_params();
}

static void format_simulation_date(int t, char *date_str) {
    init_simulation_date_from_params();

    int seconds_since_start = t - WARMUP;
    if (seconds_since_start < 0) {
        seconds_since_start = 0;
    }

    int total_days = seconds_since_start / 86400;
    int remaining_seconds = seconds_since_start % 86400;
    int hours = remaining_seconds / 3600;
    int minutes = (remaining_seconds % 3600) / 60;

    int start_year = 2000, start_month = 1, start_day = 1;  // fallback; overridden by sim_start_date
    (void)sscanf(sim_start_date, "%d-%d-%d", &start_year, &start_month, &start_day);

    int year = start_year;
    int month = start_month;
    int day = start_day + total_days;

    int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    while (day > days_in_month[month - 1]) {
        day -= days_in_month[month - 1];
        month++;
        if (month > 12) {
            month = 1;
            year++;
            if (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)) {
                days_in_month[1] = 29;
            } else {
                days_in_month[1] = 28;
            }
        }
    }

    snprintf(date_str, 32, "%04d-%02d-%02d %02d:%02d", year, month, day, hours, minutes);
}

const char *get_simulation_start_date(void) {
    init_simulation_date_from_params();
    return sim_start_date;
}

void write_field(int t, int cell_count, double *data, const char *var_name) {
    write_field_internal(t, cell_count, (const double *)data, var_name, NULL);

    // Optional centers-only output to avoid staggered-grid sawtooth artefacts in
    // longitudinal profiles: write odd indices only as a separate file.
    // Enabled via environment variable: CGEM_OUTPUT_CENTERS_ONLY=1
    if (!is_centers_only_output_enabled()) {
        return;
    }
    if (!data || !var_name || cell_count < 1) {
        return;
    }
    if (cell_count > MAXM) {
        return;
    }

    const int center_count = (cell_count + 1) / 2;
    if (center_count < 1) {
        return;
    }

    double centers[MAXM + 1];
    int k = 1;
    for (int i = 1; i <= cell_count; i += 2) {
        if (k > center_count) {
            break;
        }
        centers[k++] = data[i];
    }

    char centers_name[300];
    snprintf(centers_name, sizeof(centers_name), "%s_centers", var_name);
    write_field_internal(t, center_count, (const double *)centers, centers_name, NULL);
}

static int is_centers_only_output_enabled(void) {
    static int cached = -1;
    if (cached >= 0) {
        return cached;
    }
    const char *env = getenv("CGEM_OUTPUT_CENTERS_ONLY");
    cached = (env && env[0] != '\0' && atoi(env) != 0) ? 1 : 0;
    return cached;
}

static void write_field_internal(int t, int cell_count, const double *data, const char *var_name, const char *csv_format) {
    if (!data || !var_name) {
        return;
    }

    const bool is_binary = (output_storage_format == OUTPUT_FORMAT_BINARY);
    char path[512];
    memset(path, 0, sizeof(path));
    build_output_path(path, sizeof(path), var_name, is_binary);

    OutputStream *stream = acquire_output_stream(path, is_binary);
    if (!stream) {
        return;
    }

    FILE *fp = stream->handle;
    if (!stream->header_written) {
        if (is_binary) {
            write_binary_header(fp, (uint32_t)cell_count, (uint32_t)(cell_count + 1));
        } else {
            write_csv_header(fp, cell_count);
        }
        stream->header_written = true;
    }

    stream->write_counter++;

    if (is_binary) {
        // Write time as first value
        // IMPORTANT: For binary outputs, the time coordinate must align with
        // SIM_START_DATE, which is defined as the start of *post-warmup* forcing.
        // CSV outputs already use (t - WARMUP) via format_simulation_date().
        // Keep binary consistent so NetCDF export/validators don't get a +WARMUP shift.
        const int seconds_since_start = (t > (int)WARMUP) ? (t - (int)WARMUP) : 0;
        double time_value = (double)seconds_since_start;
        size_t written = fwrite(&time_value, sizeof(double), 1, fp);
        if (written != 1) {
            fprintf(stderr, "❌ Error: failed to write time to %s\n", path);
            perror("fwrite");
            exit(EXIT_FAILURE);
        }
        
        // Write contiguous data block directly (data[1..cell_count])
        written = fwrite(&data[1], sizeof(double), (size_t)cell_count, fp);
        if (written != (size_t)cell_count) {
            fprintf(stderr, "❌ Error: failed to write complete CGEMBIN record to %s (wrote %zu of %d values).\n",
                    path, written, cell_count);
            perror("fwrite");
            exit(EXIT_FAILURE);
        }
    } else {
        const char *format = csv_format ? csv_format : "%.10f";
        char date_buffer[32];
        format_simulation_date(t, date_buffer);
        fprintf(fp, "%s", date_buffer);
        for (int i = 1; i <= cell_count; ++i) {
            fputc(',', fp);
            fprintf(fp, format, data[i]);
        }
        fputc('\n', fp);
    }

    flush_stream(stream, false);
}

static void write_geometry_snapshot(void) {
    if (geometry_snapshot_written) {
        return;
    }

    geometry_snapshot_written = true;

    char path[256];
    snprintf(path, sizeof(path), "OUT/geometry.csv");
    ensure_output_directories();
    ensure_directory_for_path_internal(path);

    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "⚠️ Warning: Could not create %s (%s)\n", path, strerror(errno));
        return;
    }

    fprintf(fp, "Cell,Distance,RiverbedDepth,Width,BaseArea,Chezy,FRIC\n");
    for (int i = 1; i <= M; ++i) {
        double distance = (double)(i - 1) * DELXI;
        fprintf(fp, "%d,%.2f,%.4f,%.2f,%.4f,%.2f,%.6f\n",
                i,
                distance,
                riverbed_depth[i],
                width[i],
                baseArea[i],
                Chezy[i],
                FRIC[i]);
    }
    fclose(fp);
}

int Hydwrite(int t) {
    if (!is_output_enabled()) {
        return 1;
    }

    write_geometry_snapshot();

    static double velocity_buffer[MAXM + 1];
    static double discharge_buffer[MAXM + 1];

    for (int i = 1; i <= M; ++i) {
        if (i % 2 == 0) {
            velocity_buffer[i] = velocity[i];
        } else {
            if (i > 1 && i < M) {
                velocity_buffer[i] = 0.5 * (velocity[i - 1] + velocity[i + 1]);
            } else if (i == 1 && M >= 2) {
                velocity_buffer[i] = velocity[2];
            } else if (i == M && M > 1) {
                velocity_buffer[i] = velocity[M - 1];
            } else {
                velocity_buffer[i] = velocity[i];
            }
        }
        discharge_buffer[i] = Discharge(t, i, NULL);
    }

    const struct {
        const char *name;
        const double *data;
        const char *format;
    } hyd_fields[] = {
        {"velocity", velocity_buffer, "%.6f"},
        {"freeArea", freeArea, "%.6f"},
        {"waterDepth", waterDepth, "%.6f"},
        {"tau_b", tau_b, "%.6f"},
        {"eross", erosion_s, "%.8f"},
        {"deps", deposition_s, "%.8f"},
        {"erosv", erosion_v, "%.8f"},
        {"depv", deposition_v, "%.8f"},
        {"disp", disp, "%.6f"},
        {"totalArea", totalArea, "%.6f"},
        {"Discharge", discharge_buffer, "%.6f"}
    };

    const int field_count = (int)(sizeof(hyd_fields) / sizeof(hyd_fields[0]));
    for (int f = 0; f < field_count; ++f) {
        write_field_internal(t, M, hyd_fields[f].data, hyd_fields[f].name, hyd_fields[f].format);
    }

    return 1;
}

void SaveAllWaterQuality(int t) {
    if (!is_output_enabled()) {
        return;
    }

    for (int s = 0; s < CHEM_COUNT; ++s) {
        if (v[s].env != 1) {
            continue;
        }
        write_field(t, M, v[s].c, variableNames[s]);
    }
}

static void convert_path_to_varname(const char *path, char *buffer, size_t len) {
    normalise_var_name(path, buffer, len);
}

void Rates(double *co, const char *s, int t) {
    if (!enable_reaction_output) {
        return;
    }
    if (!is_output_enabled() || !co || !s) {
        return;
    }

    char var_name[256];
    convert_path_to_varname(s, var_name, sizeof(var_name));
    write_field_internal(t, M, co, var_name, NULL);
}

void Fluxwrite(int species, int t) {
    if (!enable_flux_output || !is_output_enabled()) {
        return;
    }
    if (species < 0 || species >= MAXV) {
        return;
    }
    if (v[species].env != 1 || species == pCO2 || species == PH || species == CO2) {
        return;
    }

    char adv_name[256];
    char disp_name[256];
    snprintf(adv_name, sizeof(adv_name), "Flux_Advection_%s", variableNames[species]);
    snprintf(disp_name, sizeof(disp_name), "Flux_Dispersion_%s", variableNames[species]);

    write_field_internal(t, M, v[species].advflux, adv_name, "%.10f");
    write_field_internal(t, M, v[species].disflux, disp_name, "%.10f");
}

void close_all_model_output_files(void) {
    for (int i = 0; i < output_stream_count; ++i) {
        if (output_streams[i].handle) {
            flush_stream(&output_streams[i], true);
            fclose(output_streams[i].handle);
            output_streams[i].handle = NULL;
        }

        if (output_streams[i].io_buffer) {
            free(output_streams[i].io_buffer);
            output_streams[i].io_buffer = NULL;
        }
    }
}

void reset_output_file_registry(void) {
    close_all_model_output_files();
    output_stream_count = 0;
    output_directories_initialized = false;
    geometry_snapshot_written = false;
}

// -----------------------------------------------------------------------------
// Post-processing: NetCDF export and validation scripts
// -----------------------------------------------------------------------------

#ifndef CGEM_PATH_MAX
#ifdef PATH_MAX
#define CGEM_PATH_MAX PATH_MAX
#else
#define CGEM_PATH_MAX 4096
#endif
#endif

#ifdef _WIN32
#include <process.h>
#endif

void run_postprocessing_export(void) {
    if (calibration_mode != 0) {
        return;  // NetCDF export reserved for production simulations
    }

    // Development default: do NOT export NetCDF unless explicitly requested.
    // Rationale: keeping per-variable *.bin files is faster for iteration and easier to debug.
    // Enable with: CGEM_ENABLE_NETCDF_EXPORT=1
    const char *enable_export_env = getenv("CGEM_ENABLE_NETCDF_EXPORT");
    const int enable_export = (enable_export_env && enable_export_env[0] != '\0' && atoi(enable_export_env) != 0);
    if (!enable_export) {
        if (debug_level > 0) {
            printf("ℹ️ NetCDF export skipped (set CGEM_ENABLE_NETCDF_EXPORT=1 to enable)\n");
        }
        return;
    }

    if (output_storage_format != OUTPUT_FORMAT_BINARY) {
        return;  // Nothing to convert when CSV output is enabled
    }

    char cwd[CGEM_PATH_MAX];
#ifdef _WIN32
    if (_getcwd(cwd, CGEM_PATH_MAX) == NULL) {
        if (debug_level > 0) {
            printf("⚠️ NetCDF export skipped: unable to determine working directory.\n");
        }
        return;
    }

    const char path_sep = '\\';
    const char *script_rel = "scripts\\export_to_netcdf.py";
    const char *bin_rel = "OUT";
    const char *netcdf_rel = "OUT\\cgem_outputs.nc";
    const char *params_rel = cgem_get_params_path();
#else
    if (getcwd(cwd, CGEM_PATH_MAX) == NULL) {
        if (debug_level > 0) {
            printf("⚠️ NetCDF export skipped: unable to determine working directory.\n");
        }
        return;
    }

    const char path_sep = '/';
    const char *script_rel = "scripts/export_to_netcdf.py";
    const char *bin_rel = "OUT";
    const char *netcdf_rel = "OUT/cgem_outputs.nc";
    const char *params_rel = cgem_get_params_path();
#endif

    char script_path[CGEM_PATH_MAX];
    char bin_root_path[CGEM_PATH_MAX];
    char netcdf_path[CGEM_PATH_MAX];
    char params_path[CGEM_PATH_MAX];

    if (snprintf(script_path, sizeof(script_path), "%s%c%s", cwd, path_sep, script_rel) >= (int)sizeof(script_path) ||
        snprintf(bin_root_path, sizeof(bin_root_path), "%s%c%s", cwd, path_sep, bin_rel) >= (int)sizeof(bin_root_path) ||
        snprintf(netcdf_path, sizeof(netcdf_path), "%s%c%s", cwd, path_sep, netcdf_rel) >= (int)sizeof(netcdf_path) ||
        snprintf(params_path, sizeof(params_path), "%s%c%s", cwd, path_sep, params_rel) >= (int)sizeof(params_path)) {
        if (debug_level > 0) {
            printf("⚠️ NetCDF export skipped: constructed path exceeded buffer length.\n");
        }
        return;
    }

#ifdef _WIN32
    if (_access(script_path, 0) != 0) {
#else
    if (access(script_path, F_OK) != 0) {
#endif
        if (debug_level > 0) {
            printf("⚠️ NetCDF export skipped: script not found at %s\n", script_path);
        }
        return;
    }

    const char *python_cmd = getenv("CGEM_PYTHON");
    if (!python_cmd || python_cmd[0] == '\0') {
        python_cmd = getenv("CGEM_CONDA_PYTHON");
    }
    if (!python_cmd || python_cmd[0] == '\0') {
        python_cmd = "python";
    }

#ifdef _WIN32
    char python_exec_buf[CGEM_PATH_MAX];
    strncpy(python_exec_buf, python_cmd, sizeof(python_exec_buf) - 1);
    python_exec_buf[sizeof(python_exec_buf) - 1] = '\0';

    size_t exec_len = strlen(python_exec_buf);
    if (exec_len >= 2 && python_exec_buf[0] == '"' && python_exec_buf[exec_len - 1] == '"') {
        memmove(python_exec_buf, python_exec_buf + 1, exec_len - 2);
        python_exec_buf[exec_len - 2] = '\0';
    }

    char script_arg[CGEM_PATH_MAX * 2];
    char bin_root_arg[CGEM_PATH_MAX * 2];
    char netcdf_arg[CGEM_PATH_MAX * 2];
    char params_arg[CGEM_PATH_MAX * 2];

    if (snprintf(script_arg, sizeof(script_arg), "\"%s\"", script_path) >= (int)sizeof(script_arg) ||
        snprintf(bin_root_arg, sizeof(bin_root_arg), "\"%s\"", bin_root_path) >= (int)sizeof(bin_root_arg) ||
        snprintf(netcdf_arg, sizeof(netcdf_arg), "\"%s\"", netcdf_path) >= (int)sizeof(netcdf_arg) ||
        snprintf(params_arg, sizeof(params_arg), "\"%s\"", params_path) >= (int)sizeof(params_arg)) {
        if (debug_level > 0) {
            printf("⚠️ NetCDF export skipped: argument quoting exceeded buffer length.\n");
        }
        return;
    }

    const char *argv_spawn[] = {
        python_exec_buf,
        script_arg,
        "--bin-root",
        bin_root_arg,
        "--netcdf",
        netcdf_arg,
        "--params",
        params_arg,
        NULL
    };

    int status = _spawnvp(_P_WAIT, python_exec_buf, (const char * const *)argv_spawn);
    if (status != 0) {
        if (debug_level > 0) {
            printf("⚠️ NetCDF export failed: _spawnvp returned %d (errno=%d).\n", status, errno);
        }
    }
#else
    /* Use fork+execvp instead of system() to avoid shell injection risks. */
    {
        pid_t pid = fork();
        if (pid < 0) {
            if (debug_level > 0) {
                printf("⚠️ NetCDF export skipped: fork() failed (errno=%d).\n", errno);
            }
            return;
        }
        if (pid == 0) {
            /* Child process */
            execlp(python_cmd, python_cmd,
                   script_path,
                   "--bin-root", bin_root_path,
                   "--netcdf",  netcdf_path,
                   "--params",  params_path,
                   (char *)NULL);
            _exit(127); /* exec failed */
        }
        /* Parent: wait for child */
        int wstatus = 0;
        waitpid(pid, &wstatus, 0);
        if (!WIFEXITED(wstatus) || WEXITSTATUS(wstatus) != 0) {
            if (debug_level > 0) {
                printf("⚠️ NetCDF export script returned non-zero status.\n");
            }
        }
    }
#endif
}
