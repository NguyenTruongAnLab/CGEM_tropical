# Getting Started

This guide walks you through building C-GEM, running the bundled example
configuration, and adapting the model to your own estuary. The key idea is that
the model is configured through `INPUT/params.txt` and
`INPUT/config_input.txt`, so most case changes do **not** require source-code
edits.

---

## 1. Prerequisites

| Tool | Version | Notes |
|------|---------|-------|
| C compiler | C99 | GCC, Clang, or MSVC |
| CMake | ≥ 3.15 | [cmake.org](https://cmake.org/) |
| NLopt | *(optional)* | Only needed for automated calibration |

## 2. Build

```bash
# Configure
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build build --config Release
```

The executable is placed in `bin/Release/C_GEM_Daily` (or `.exe` on Windows).

To enable calibration (requires NLopt):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCGEM_ENABLE_CALIBRATION=ON
cmake --build build --config Release
```

## 3. Run the Example

```bash
# Linux / macOS
./bin/Release/C_GEM_Daily

# Windows
bin\Release\C_GEM_Daily.exe
```

The model reads `INPUT/params.txt` and `INPUT/config_input.txt` by default.
The repository ships a Saigon River input dataset in `INPUT/`, but the solver
itself is case-agnostic: swapping estuaries is mainly a matter of editing these
configuration files and replacing the referenced forcing/geometry files.
Results are written to `OUT/`.

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `CGEM_PARAMS_PATH` | `INPUT/params.txt` | Alternate parameter file |
| `CGEM_TRANSPORT_ONLY` | `0` | Set to `1` for transport-only mode (skip biogeochemistry) |
| `CGEM_MASS_DIAG_DAYS` | *(all)* | Limit simulation to N days (useful for quick checks) |
| `CGEM_WARMUP_DAYS` | *(from params.txt)* | Override warmup duration (days) |
| `CGEM_MAXT_DAYS` | *(from params.txt)* | Override total simulation time (days) |
| `CGEM_FORCE_NO_CALIBRATION` | `0` | Set to `1` to skip calibration even if enabled |
| `CGEM_DEBUG_LEVEL` | `1` | Verbosity (0=silent, 1=essential, 2=detail, 3=verbose) |

## 4. Understand the Output

The model writes per-variable output files to `OUT/`, for example:

- **State variables and hydrodynamics** — `OUT/Sal.csv`, `OUT/O2.csv`,
   `OUT/velocity.csv`, `OUT/waterDepth.csv`, etc.
- **Flux diagnostics** — `OUT/Flux_Advection_<Var>.csv` and
   `OUT/Flux_Dispersion_<Var>.csv` when `enable_flux_output = 1`
- **Reaction diagnostics** — `OUT/Reaction_*.csv` and `OUT/Diag_*.csv` when
   `enable_reaction_output = 1`
- **Geometry snapshot** — `OUT/geometry.csv`
- **Derived figures/tables** — products created later by the Python tools under
   `OUT/figures/` and `OUT/tables/`

By default the public configuration writes CSV files. Advanced users can switch
to binary output with `output_storage_format = BINARY` in `INPUT/params.txt`.
The public Python tools read either CSV or CGEM binary output transparently.

## 5. Adapt to Your Estuary

### Step 1: Geometry

Create a bed-profile file (depth vs. distance from mouth) and place it in
`INPUT/Geometry/`.  Update `params.txt`:

```ini
EL = 150000             # Estuary length [m]
DELXI = 1500            # Grid spacing [m]
AMPL = 2.0              # Tidal amplitude at mouth [m]
riverbed_profile_file = "INPUT/Geometry/my_river_bed.txt"
```

### Step 2: Boundary Conditions

Prepare two-column numeric forcing files (`time_index,value`) for upstream and
downstream boundaries. The bundled example uses sequential counters in the
first column:

- **Upstream:** discharge, salinity, and water-quality concentrations
- **Downstream:** tidal elevation (hourly), salinity, marine WQ

Update file paths in `INPUT/config_input.txt`.

### Step 3: Tributaries

For each tributary, create a block in `config_input.txt`:

```ini
name=MyTributary
type=Tributary
resolution=daily
cellIndex=25          # Grid cell where tributary enters
discharge=INPUT/Tributaries/MyTrib/discharge.csv
NO3=INPUT/Tributaries/MyTrib/NO3.csv
# ... other variables
```

Set `numTributaries` at the top of the file.

### Step 4: Segment Geometry

Define channel segments in `params.txt`.  Each segment requires width at the
segment start (`B`), convergence length (`LC`), Chézy roughness, and storage
ratio:

```ini
num_segments = 2
index_1 = 0
B1 = 3000.0
LC1 = 40000
Chezy1 = 50.0
Rs1 = 1.2
index_2 = 30
B2 = 300.0
LC2 = 80000
Chezy2 = 35.0
Rs2 = 1.0
```

### Step 5: Calibrate

C-GEM uses a four-stage calibration workflow:

| Stage | Target | What is tuned |
|-------|--------|---------------|
| 1 | Tidal range + salinity | Chézy, Rs, C_VDB, D0, AMPL correction |
| 2 | SPM profiles | Erosion/deposition parameters |
| 3 | Nutrients, O₂, Chl-a | Biogeochemical rates |
| 4 | pCO₂, pH, DIC, AT | Gas exchange scaling |

Set `calibration_mode` in `params.txt` to an NLopt algorithm (1–7), set the
`calibration_stage`, and run.  Parameter bounds are defined in
`INPUT/Calibration/calibration_parameters.csv`.

### Step 6: Validate

Compare output profiles and time series against independent field data.
Place observation files in `INPUT/Validation/`.

## 6. Run Tests

```bash
ctest --test-dir build -C Release
```

The public repository always registers a smoke test. Some development branches
may also register additional unit tests when a `tests/` directory is present.

---

## 7. For Developers

### Adding a new biogeochemical variable

1. Add an entry to `enum Chem` in [src/variables.h](../../src/variables.h)
   and increment `CHEM_COUNT`.
2. Set the variable name in `Init()` inside [src/main.c](../../src/main.c)
   (e.g. `strcpy(v[MyVar].name, "MyVar")`).
3. Add boundary-condition loading in [src/bcforcing.c](../../src/bcforcing.c).
4. Add source/sink terms in `Biogeo()` inside [src/biogeo.c](../../src/biogeo.c).
5. Add a test under `tests/` and register it in [CMakeLists.txt](../../CMakeLists.txt).

### Key CMake options

| Option | Default | Effect |
|--------|---------|--------|
| `CGEM_ENABLE_CALIBRATION` | `OFF` | Link NLopt and compile calibration module |
| `CMAKE_BUILD_TYPE` | `Release` | Build type (`Debug` for asserts/symbols) |

### Code style

- C99, no external dependencies beyond standard library (and optional NLopt).
- One function per concern: `Hyd()` for hydrodynamics, `Transport()` for
  advection-dispersion, `Biogeo()` for reactions.
- All global state declared `extern` in `variables.h`, defined in `variables.c`.

---

## File Reference

| File | Purpose |
|------|---------|
| `INPUT/params.txt` | All model parameters |
| `INPUT/config_input.txt` | Boundary/tributary file mapping |
| `INPUT/README.md` | Detailed description of the example dataset |
| `src/README.md` | Source-code architecture overview |
