# C-GEM Tropical

A 1D reactive-transport model for **data-sparse tropical estuaries**, coupling hydrodynamics (Saint-Venant), advection–dispersion (TVD + Crank–Nicolson), and biogeochemistry (NPZD, carbonate chemistry, gas exchange). Built on the C-GEM framework (Volta et al., 2014; 2016).

The source code is **generic** — all site-specific information (geometry, tributary layout, forcing data, parameters) is defined in the `INPUT/` directory. The included configuration targets the **Saigon–Dong Nai estuary** (Vietnam), but the model can be applied to any tide-dominated tropical system by supplying new input files.

---

## Quick start

### Option 1: Portable release (no compiler needed)

Pre-built binaries for **Windows** and **Linux** are available on the [Releases](../../releases) page.
Download the archive for your platform, extract it, and run:

```bash
# Linux
./C_GEM_Daily

# Windows
C_GEM_Daily.exe
```

The release archives bundle the executable together with `INPUT/`, `tools/`,
and `docs/` so reviewers can run the model and make figures without rebuilding.

To adapt the model to a different estuary, edit the files in the `INPUT/` folder — no source code changes needed.

### Option 2: Build from source

Requires **CMake ≥ 3.15** and a C99 compiler (GCC, Clang, or MSVC).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

To enable the optional calibration module (requires [NLopt](https://github.com/stevengj/nlopt)):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCGEM_ENABLE_CALIBRATION=ON
```

> On Windows/MSVC, the explicit `-C Release` flag is required for `ctest`.
> The public repository always registers a smoke test; additional unit tests are
> built automatically when a `tests/` directory is present.

### Run a quick solver check

```bash
# Transport-only mode, 30 simulated days
cmake -E env CGEM_TRANSPORT_ONLY=1 CGEM_MASS_DIAG_DAYS=30 bin/Release/C_GEM_Daily
```

## Reviewer workflow

After building the model, a new user can reproduce the full software workflow in four steps:

1. Run the included Saigon–Dong Nai example.
2. Check forcing-data quality with `python tools/check_inputs.py`.
3. Plot outputs with the scripts in `tools/`:
  - `plot_profiles.py`
  - `plot_seasonal.py`
  - `plot_hydrodynamics.py`
  - `plot_reactions.py`
4. Overlay the bundled validation observations by adding `--with-obs`.

Install the analysis dependencies with:

```bash
pip install -r tools/requirements.txt
```

The observation tables already bundled in `INPUT/Validation/` are in model
units and work directly with `--with-obs`.

## Prepare input data

To adapt the model to a new estuary, the main work happens in `INPUT/`:

- `INPUT/params.txt` — physical, transport, sediment, and biogeochemical parameters
- `INPUT/config_input.txt` — file mapping for boundaries, forcings, and tributaries
- `INPUT/Boundary/` — upstream/downstream forcing files
- `INPUT/Tributaries/` — discharge and chemistry for each tributary
- `INPUT/Geometry/river_bed_profile.txt` — longitudinal bed profile
- `INPUT/Validation/` — independent observations for figure overlays and validation

The shipped forcing files use **two-column numeric time series**. See
`INPUT/README.md` before replacing them.

## Understand the outputs

After each run, the model writes outputs to `OUT/`:

- one file per simulated variable (for example `Sal.csv`, `O2.csv`, `DIC.csv`)
- hydrodynamic diagnostics such as `velocity.csv`, `waterDepth.csv`, and `disp.csv`
- `geometry.csv` for downstream plotting
- figures and tables produced by the Python tools under `OUT/figures/` and `OUT/tables/`

For a first review pass, the most useful next commands are:

```bash
python tools/plot_profiles.py --output-dir OUT --with-obs
python tools/plot_hydrodynamics.py --output-dir OUT --with-obs
python tools/plot_seasonal.py --output-dir OUT --with-obs
```

---

## Project structure

```
C_GEM_Tropical/
│
├── CMakeLists.txt              Build system (CMake)
├── README.md                   This file
├── LICENSE                     MIT License
├── docs/                       Extended scientific and implementation docs
├── src/                        Model source code (C99)
│   ├── main.c                  Entry point
│   ├── hydrodynamics.c         Saint-Venant tidal solver
│   ├── transport.c / .h        TVD advection–dispersion solver
│   ├── biogeo.c / .h           NPZD biogeochemistry + carbonate chemistry
│   ├── biogeout.c              Biogeochemical rate diagnostics
│   ├── bcforcing.c             Boundary-condition interpolation (tidal harmonics)
│   ├── init.c / .h             Model initialization + default rate constants
│   ├── io_input.c              Parameter & forcing file readers
│   ├── io_output.c             Time-series & profile output writers
│   ├── diagnostics.c / .h      Runtime sanity checks (NaN, negative, bounds)
│   ├── variables.c / .h        Global state arrays and counters
│   ├── utilities.c / .h        Math helpers (interpolation, date parsing)
│   ├── define.h                Preprocessor constants (array sizes, physics)
│   ├── file.h                  Cross-platform file-path helpers
│   ├── calibration.c / .h      Multi-stage calibration driver (optional)
│   ├── calibration_helpers.c/.h  Observation matching & scoring (optional)
│   ├── calibration_stub.c      No-op stub when calibration is disabled
│   └── nlopt_winshim.c         Windows DLL loader for NLopt
│
├── INPUT/                      Model configuration & forcing data
│   ├── params.txt              Calibrated parameters (grid, physics, biogeochem)
│   ├── config_input.txt        File paths for all boundary & tributary forcing
│   ├── README.md               Detailed description of every input file
│   ├── Boundary/
│   │   ├── UB/                 Upstream boundary (discharge, salinity, WQ, SPM)
│   │   ├── UB_stable/          Stable upstream endmember concentrations
│   │   ├── LB/                 Downstream boundary (tidal, light, WQ, carbonate)
│   │   └── wind.csv            Wind speed
│   ├── Geometry/
│   │   └── river_bed_profile.txt   Longitudinal bed elevation
│   ├── Tributaries/
│   │   ├── Dongnai/            Major tributary (cellIndex 31)
│   │   ├── Canals/             Urban canal system (cellIndex 37)
│   │   ├── VamThuat/           Urban tributary (cellIndex 45)
│   │   └── ThiTinh/            Peri-urban tributary (cellIndex 61)
│   ├── carbonate/              Derived DIC/AT from field pCO₂ + alkalinity
│   ├── Calibration/            Parameter bounds & observation targets
│   └── Validation/             Independent observation datasets
│
├── tests/                      Unit & integration tests (CTest, 21 tests)
│   ├── test_transport_*.c      Transport solver (mass conservation, boundaries,
│   │                           dispersion, tributary mixing)
│   ├── test_biogeo_*.c         Biogeochemistry (growth limiting, CO₂ exchange,
│   │                           O₂/TOC/Phy budgets, carbonate–nutrient coupling)
│   ├── test_hydrodynamics_*.c  Hydrodynamic source-term scaling
│   ├── test_calibration_*.c    Calibration scoring & time-base alignment
│   ├── test_forcing_*.c        Forcing continuity across warmup boundary
│   └── test_abort_on_*.c       Defensive checks (NaN, negative concentrations)
│
├── external/
│   └── nlopt/                  NLopt headers (optional; for calibration only)
│
├── tools/                      Python helpers for QC, plotting, validation
│
└── OUT/                        Model outputs (generated at runtime, not tracked)
```

## Source code architecture

The model executes a **daily time loop** with sub-daily tidal cycling:

```
main.c
  ├── io_input.c     → read params.txt, config_input.txt, all forcing CSVs
  ├── init.c         → allocate arrays, set default rate constants, build grid
  │
  ├── [daily loop]
  │   ├── bcforcing.c      → interpolate tidal & WQ boundary conditions
  │   ├── hydrodynamics.c  → solve Saint-Venant for water level & velocity
  │   ├── transport.c      → advect–disperse all state variables (TVD scheme)
  │   ├── biogeo.c         → NPZD reactions, carbonate equilibrium, gas exchange
  │   ├── diagnostics.c    → runtime checks (NaN, negative, bounds)
  │   └── biogeout.c       → log biogeochemical rate diagnostics
  │
  ├── io_output.c    → write time-series, profiles, mass budgets
  └── calibration.c  → (optional) multi-stage NLopt parameter optimization
```

**Key design principles:**

- All state is in global arrays (`variables.c`) — flat C, no OOP.

- Site-specific configuration lives entirely in `INPUT/` — the source code has no hardcoded estuary geometry.

- Boundary conditions are CSV time-series, interpolated at runtime by `bcforcing.c`.

- The calibration module is compile-time optional (`-DCGEM_ENABLE_CALIBRATION=ON/OFF`).

---

## Adapting to a new estuary

1. **Geometry**: edit `INPUT/Geometry/river_bed_profile.txt` with your cross-sections.

2. **Grid**: set `ncells`, `dx`, `L` in `INPUT/params.txt`.

3. **Boundaries**: replace CSV files in `INPUT/Boundary/UB/` and `INPUT/Boundary/LB/` with your upstream/downstream forcing.

4. **Tributaries**: add/remove tributary folders under `INPUT/Tributaries/` and update `config_input.txt` with file paths and cell indices.

5. **Parameters**: adjust biogeochemical rates in `INPUT/params.txt` (or run the calibration module).

See `docs/Documentation/Getting-Started.md` for a detailed walkthrough.

## Environment variables

| Variable                | Description                                               |
| ----------------------- | --------------------------------------------------------- |
| `CGEM_TRANSPORT_ONLY=1` | Run hydrodynamics + transport only (skip biogeochemistry) |
| `CGEM_MASS_DIAG_DAYS=N` | Stop after N simulated days (useful for quick tests)      |

## License

This project is licensed under the MIT License. See `LICENSE`.

## Status and citation

This repository is the **public software archive** for C-GEM Tropical.

- it contains the model source code, example input data, analysis tools, and technical documentation
- it is intended to be versioned through **GitHub releases** and archived through **Zenodo**
- the accompanying **Water Research** manuscript is being prepared/submitted separately and is **not** intended to be published as part of this public GitHub repository

For formal software citation, please cite the **Zenodo DOI corresponding to the exact release version used**.

The root `CITATION.cff` file is retained for GitHub metadata and machine-readable
software information, but the preferred journal-facing citation target is the
versioned **Zenodo record**.

Once the first Zenodo DOI is minted, add that DOI prominently to this section
and to the GitHub release notes.