# Chapter 6 — Implementation

This chapter describes the actual source code architecture of C-GEM: how the
files are organized, how data flows through the simulation, and how to extend
or modify the model.

---

## 6.1 Source files

| File | Purpose |
|------|---------|
| `main.c` | Entry point — parses environment variables, calls `Init()`, runs the main time loop |
| `init.c` / `init.h` | Model initialization: reads `params.txt` and `config_input.txt`, allocates arrays, sets default rate constants, builds the grid |
| `variables.c` / `variables.h` | Global state: all arrays, structs, enums, and extern declarations |
| `define.h` | Preprocessor constants (`MAXV`, `MAXM`, `MAX_TIMESERIES`, `TIDAL_PERIOD`, `G`, etc.) |
| `hydrodynamics.c` | Saint-Venant solver: `Hyd()` → `Newbc()` → `Coeffa()` → `NewUH()` → `Update()` → `Tridag()` |
| `transport.c` / `transport.h` | TVD advection + Crank–Nicolson dispersion: `Transport()` → `Adv()` / `Disp()` / `TVD()` |
| `biogeo.c` / `biogeo.h` | NPZD reactions, carbonate chemistry, gas exchange: `Biogeo()` |
| `biogeout.c` | Biogeochemical rate accessors (`Pbmax()`, `kmort()`, `Fnit()`, `O2sat()`, etc.) |
| `bcforcing.c` | Boundary-condition interpolation: `bgboundary()`, `Tide()`, `Discharge_ups()` |
| `io_input.c` | Parameter file reader (`read_parameters()`) and numeric two-column forcing loader (`readFile()`) |
| `io_output.c` | Output writers: per-variable CSV/BINARY streams, geometry snapshot, and diagnostics under `OUT/` |
| `diagnostics.c` / `diagnostics.h` | Runtime sanity checks: NaN detection, bound violations, numerical stability |
| `utilities.c` / `utilities.h` | Math helpers: interpolation, date parsing, moving averages |
| `file.h` | Cross-platform file path macros |
| `calibration.c` / `calibration.h` | Multi-stage NLopt calibration driver (optional, compile-time) |
| `calibration_helpers.c` / `calibration_helpers.h` | Observation matching and KGE/RMSE scoring for calibration |
| `calibration_stub.c` | No-op fallback when calibration is disabled (`-DCGEM_ENABLE_CALIBRATION=OFF`) |
| `nlopt_winshim.c` | Windows DLL runtime loader for NLopt |

---

## 6.2 Main simulation loop

The core loop in `main.c` runs at a fixed time step (`DELTI`, configured in
`INPUT/params.txt`; the bundled public example uses `360 s`):

```c
// From main.c (simplified for clarity)
Init();

for (long t = 0; t <= MAXT; t += DELTI) {
    Hyd(t);             // Saint-Venant: water level + velocity
    bgboundary(t);      // Interpolate boundary concentrations
    Transport(t);       // TVD advection + dispersion for all species

    check_numerical_stability(t);  // NaN / bound checks

    if (enable_biogeochemical_reactions || enable_carbonate_diagnostics) {
        Biogeo(t);      // NPZD reactions + carbonate equilibrium
    }

    updateSuspendedSediment(t);    // Erosion / deposition
    SaveAllWaterQuality(t);        // Write state output at configured intervals
    // Fluxwrite(species, t) is then called for each enabled transported variable
}
```

**Execution order matters**: hydrodynamics must run first (provides velocity
and water depth), then boundary conditions are interpolated, then transport
advects/disperses all species, and finally biogeochemistry modifies
concentrations in place.

---

## 6.3 Key data structures

### State variables: `struct Verb`

Each transported variable (Phy1, NO3, O2, DIC, etc.) is stored in a `Verb` struct:

```c
// From variables.h
struct Verb {
    char name[20];              // Display name (e.g., "Phy1")
    int env;                    // Active flag
    double c[MAXM + 1];        // Concentration at each grid cell
    double clb;                 // Lower boundary value
    double cub;                 // Upper boundary value
    double avg[MAXM + 1];      // Running average
    double concflux[MAXM + 1]; // Concentration flux
    double advflux[MAXM + 1];  // Advective flux
    double disflux[MAXM + 1];  // Dispersive flux
};

extern struct Verb v[MAXV];     // MAXV = CHEM_COUNT
```

Concentrations are accessed as `v[species].c[cell]`, where `species` is an
enum value (e.g., `v[NO3].c[50]` gives NO₃ at cell 50).

### Species enum

```c
// From variables.h
typedef enum {
    Phy1, Phy2, Si, NO3, NH4, PO4,
    PIP,   // Particulate Inorganic Phosphorus
    O2, TOC, Sal, SPM, DIC, AT,
    pCO2,  // auto-calculated from DIC + AT
    PH,    // auto-calculated
    CO2,   // auto-calculated
    CHEM_COUNT  // = 17
} Chem;
```

Four species (`pCO2`, `PH`, `CO2`, `PIP`) are auto-calculated and not
independently transported. The remaining 13 are solved by the transport
equation.

### Hydrodynamic arrays

Hydrodynamic state uses simple global arrays (not structs):

| Array | Description | Unit |
|-------|-------------|------|
| `level[MAXM+1]` | Water level | m |
| `velocity[MAXM+1]` | Flow velocity | m/s |
| `waterDepth[MAXM+1]` | Water depth | m |
| `totalArea[MAXM+1]` | Cross-sectional area | m² |
| `width[MAXM+1]` | Channel width | m |
| `Chezy[MAXM+1]` | Chézy roughness coefficient | m^½/s |
| `C[MAXM+1][5]` | Tridiagonal matrix coefficients | — |

### Boundary conditions

Each species has upstream and downstream boundary arrays:

```c
// From variables.h
typedef struct {
    double *time;       // Time stamps (seconds since simulation start)
    double *data;       // Concentration values
    const char *filePath;
    int required;
    int dataSize;
} BCArray;

extern BCArray upstreamBC[CHEM_COUNT];
extern BCArray downstreamBC[CHEM_COUNT];
```

At each time step, `bgboundary()` interpolates these arrays to get the current
boundary concentration.

### Tributaries

```c
// From variables.h
typedef struct {
    char name[50];
    int cellIndex;                         // Grid cell where tributary enters
    char dischargeFile[200];
    double *dischargeTime, *discharge;
    int dischargeDataSize;
    struct {
        char filePath[200];
        double *timeArray, *dataArray;
        int dataSize;
    } chemicalData[CHEM_COUNT];           // Per-species forcing
    double concentration[CHEM_COUNT];      // Current interpolated values
    int is_urban;
} Tributary;

extern Tributary *tributaries;
extern int numTributaries;
```

Tributaries inject mass using the Fischer (1979) instantaneous-mixing formula,
implemented in `applyTributarySourceTerms()` in `transport.c`.

---

## 6.4 Configuration files

### `params.txt` — model parameters

Key-value format, one parameter per line. Lines starting with `#` or `!` are
comments. Example entries:

```ini
# Grid
EL = 202000           # Estuary length [m]
DELXI = 2000          # Grid spacing [m]
DELTI = 360           # Time step [s]

# Tidal forcing
AMPL = 1.5            # Tidal amplitude at mouth [m]

# Segment geometry (up to 4 segments)
num_segments = 3
index_1 = 0
B1 = 5850.0           # Width at segment 1 start [m]
LC1 = 43000           # Convergence length [m]
Chezy1 = 57.0         # Chézy coefficient [m^(1/2)/s]
Rs1 = 1.0             # Storage width ratio [-]

# Biogeochemistry
Pb_Phy1 = 2.0e-4      # Max growth rate phytoplankton [s⁻¹]
KN_Phy1 = 25.0        # Half-saturation N [mmol N/m³]
kbg = 0.2             # Background light attenuation [m⁻¹]
```

The full parameter list is read by `read_parameters()` in `io_input.c`.
Any parameter not found in the file keeps its default value from
`assignBiogeochemicalRateConstants()` in `init.c`.

For most applications, adapting the model means editing `INPUT/params.txt` and
`INPUT/config_input.txt`, not changing the source code. Geometry, forcing-file
paths, runtime switches, calibration settings, and output format are all driven
from these configuration files.

### `config_input.txt` — forcing file paths

Specifies which forcing files provide data for boundaries and tributaries.
Format example:

```ini
# Upstream boundary
discharge=INPUT/Boundary/UB/discharge.csv
Sal=INPUT/Boundary/UB/Sal.csv
NO3=INPUT/Boundary/UB/NO3.csv

# Downstream boundary
elevation=INPUT/Boundary/LB/elevation.csv
light=INPUT/Boundary/LB/Light.csv

# Tributaries
numTributaries=4

name=DongNai
cellIndex=31
discharge=INPUT/Tributaries/Dongnai/discharge.csv
NO3=INPUT/Tributaries/Dongnai/NO3.csv
```

### Forcing CSV format

All forcing files use the same two-column format:

```csv
1,15.6
2,14.8
3,16.2
```

Column 1 is a **numeric time index**, not a calendar-date string; column 2 is
the value. The temporal meaning of that index comes from the corresponding
block in `config_input.txt` (`resolution=daily` or `resolution=hourly`). Units
must match what the model expects (see `INPUT/README.md`).

---

## 6.5 Output

All output is written to the `OUT/` directory. The public runtime writes one
top-level file per variable/diagnostic, using either `.csv` or `.bin`
depending on `output_storage_format`:

| File pattern | Content |
|-------------|---------|
| `OUT/<Var>.csv` or `OUT/<Var>.bin` | State variables and hydrodynamics (for example `Sal`, `O2`, `DIC`, `velocity`, `waterDepth`, `Discharge`) |
| `OUT/Flux_Advection_<Var>.*` and `OUT/Flux_Dispersion_<Var>.*` | Advective and dispersive flux diagnostics when `enable_flux_output = 1` |
| `OUT/Reaction_<Name>.*` and `OUT/Diag_<Name>.*` | Reaction and limitation diagnostics when `enable_reaction_output = 1` |
| `OUT/geometry.csv` | Geometry snapshot used by plotting tools |
| `OUT/figures/` and `OUT/tables/` | Analysis products generated later by the Python scripts |

The bundled Python tools read both CSV and CGEM binary output through
`tools/cgem_read.py`, so the post-processing workflow stays the same when you
switch storage format.

---

## 6.6 Build system

The project uses CMake with one optional feature flag:

```bash
# Without calibration
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCGEM_ENABLE_CALIBRATION=OFF

# With calibration (requires NLopt)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCGEM_ENABLE_CALIBRATION=ON
```

When calibration is disabled, `calibration_stub.c` provides no-op function
stubs so the rest of the code compiles unchanged. When enabled,
`calibration.c` and `calibration_helpers.c` are compiled and linked against
NLopt.

### Cross-platform support

- **Linux/macOS**: GCC or Clang, standard CMake build
- **Windows**: MSVC via Visual Studio or CMake. NLopt is loaded at runtime via
  `nlopt_winshim.c` (no link-time dependency)
- All source is C99. No platform-specific APIs are used in the core solver.

---

## 6.7 Testing

The project includes 21 CTest tests in the `tests/` directory:

```bash
ctest --test-dir build -C Release --output-on-failure
```

Tests cover:
- **Transport**: mass conservation, boundary conditions, dispersion–discharge
  parity, tributary Fischer mixing, control-volume length
- **Hydrodynamics**: tributary source-term scaling
- **Biogeochemistry**: growth limiting, CO₂ exchange sign, O₂/TOC/Phy budget
  parity, carbonate–nutrient alkalinity coupling
- **Calibration**: score breakdown, time-base alignment, target summary
- **Forcing**: continuity across warmup boundary, light/temperature consistency
- **Defensive**: abort on NaN/Inf, abort on negative concentrations

Each test is a standalone C file that links against the model source and
exercises one aspect of the solver.

---

## 6.8 Calibration (optional)

When compiled with `-DCGEM_ENABLE_CALIBRATION=ON`, the model supports
automated multi-stage parameter optimization using NLopt:

| Stage | Target | Parameters tuned |
|-------|--------|-----------------|
| 1 | Tidal range + salinity | Chézy, Rs, C_VDB, D0, AMPL correction (9 params) |
| 2 | SPM profiles | Erosion/deposition rates (15 params) |
| 3 | Nutrients, O₂, Chl-a | Biogeochemical rates (20 params) |
| 4 | pCO₂, pH, DIC, AT | Gas exchange scaling (1 param) |

**Calibration workflow:**

1. Set `calibration_mode` in `params.txt` to an NLopt algorithm (e.g., `5` for BOBYQA).
2. Set `calibration_stage` to the stage number (1–4).
3. Define parameter bounds in `INPUT/Calibration/calibration_parameters.csv`.
4. Provide observation targets in `INPUT/Calibration/targets_tidal.csv` (Stage 1)
   or `INPUT/Calibration/targets_seasonal_unified.csv` (Stages 2–4).
5. Run the model. The calibration driver (`run_hydrodynamic_calibration_from_params()`
   in `calibration.c`) iterates the solver, computes KGE/RMSE scores against
   observations, and writes the best parameters to `params_calibrated.txt`.

Calibration can be disabled at runtime with `CGEM_FORCE_NO_CALIBRATION=1`.

---

## 6.9 Extending the model

### Adding a new state variable

1. Add the species to the `Chem` enum in `variables.h` (before `CHEM_COUNT`).
2. Initialize its default parameters in `assignBiogeochemicalRateConstants()` in `init.c`.
3. Add reaction terms in `Biogeo()` in `biogeo.c`.
4. Provide upstream/downstream boundary CSVs and register paths in `config_input.txt`.
5. Add a test in `tests/` to verify the new species budget.

### Adding a new tributary

1. Create a folder under `INPUT/Tributaries/` with forcing CSVs.
2. Add a block in `INPUT/config_input.txt` with the tributary name, cell index,
   and file paths.
3. Increment `numTributaries` at the top of `config_input.txt`.

### Modifying the geometry

Edit `INPUT/Geometry/river_bed_profile.txt` (distance vs. depth) and adjust
segment parameters (`B`, `LC`, `Chezy`, `Rs`) and grid settings (`EL`, `DELXI`)
in `INPUT/params.txt`.
