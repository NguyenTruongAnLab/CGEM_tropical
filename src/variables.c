/**
 * @file variables.c
 * @brief Global variable definitions for C-GEM.
 */

#include "variables.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ============================
// Global Variable Definitions
// ============================

/** @brief Model flags: calibration_mode (algorithm: 0=none, 1=BOBYQA, 2=Nelder-Mead,
 *                    3=SBPLX, 4=COBYLA, 5=DIRECT-L, 6=CRS2-LM, 7=ISRES),
 *                    debug_level (0-2). */
int calibration_mode = 0;
int debug_level = 1;      // Default to ESSENTIAL level (errors and warnings)
OutputFormat output_storage_format = OUTPUT_FORMAT_CSV; // Default to CSV outputs
int enable_flux_output = 1;       // Default: emit flux diagnostics
int enable_reaction_output = 1;   // Default: emit reaction diagnostics


/** @brief Automatically align segment break with dominant tributary unless
 * disabled. */
int align_segment_break_with_major_tributary = 1;

/** @brief Calibration controller defaults (configuration-driven). */
double calibration_warmup_days = 100.0;     // Minimum 100 days for proper equilibration
double calibration_diagnostic_days = 200.0; // 200 days analysis period for stable statistics
int calibration_max_iterations = 200;       // Up to 200 iterations for thorough optimization



/** @brief Number of river segments (typical: 1-3). */
int num_segments = 1;
int segment_transition_width = 12; // Default 12 cells (24 km) half-width
double dispersion_abs_min = 10.0;  // Absolute floor [m2/s]
double dispersion_abs_max = 5000.0; // Absolute ceiling [m2/s]

// === Literature-Based Dispersion Physics Parameters (Gisen et al. 2015) ===

/** @brief Dispersion calibration parameters (Van der Burgh 1972 formulation).
 */
double C_VDB = 4.38;        // Shape factor (literature default)
double D0_CORRECTION = 2.0; // Magnitude correction (+100% for improved mouth mixing, MAINT-001A fix)
double K_DISCHARGE_SENSITIVITY = 0.0; // Discharge-dependent K correction (0 = disabled)
// === Tidal Physics Parameters (Astronomical Constants) ===
/** @brief M2 tidal period (hours) - astronomical constant. */
double TIDAL_PERIOD_M2_HOURS = 12.42;

/** @brief S2 tidal period (hours) - astronomical constant. */
double TIDAL_PERIOD_S2_HOURS = 12.00;

/** @brief M2 phase (π/2 radians) - configurable by location. */
double TIDAL_PHASE_M2_RAD = 1.5708;

/** @brief S2 phase (π/2 radians) - configurable by location. */
double TIDAL_PHASE_S2_RAD = 1.5708;

// --- Model Time & Grid ---

/** @brief Time settings: MAXT (total sim. time), WARMUP (equilibration), DELTI
 * (time step), TS (save interval) [seconds]. */
long MAXT = 0;
long WARMUP = 0;
int DELTI = 0;
int TS = 0;
long simulation_window_seconds = 0;

/** @brief Standard biogeochemical computation - no performance shortcuts. */

// Mass balance validation mode
int enable_biogeochemical_reactions = 1;
int enable_carbonate_diagnostics = 1;
int enable_suspended_sediment_dynamics = 1;
double biogeochem_cache_stability_threshold = 0.05; // Cache when <5% change

int enable_upstream_discharge_smoothing = 1;
double upstream_discharge_smoothing_timescale_days = 2.0;
double upstream_discharge_scale = 1.0;
int enable_tributary_subdaily_pulses = 0;
double tributary_pulse_shape = 0.0;
double urban_tributary_pulse_boost = 1.0;
double spm_upstream_scale = 1.0;
TributaryMassInjectionMode tributary_mass_injection_mode =
  TRIB_INJECTION_FISCHER;

// Quadratic mortality closure for implicit zooplankton grazing
// Following Fennel et al. (2006, JGR-Oceans), Edwards & Yool (2000)
double kmort2_Phy1 = 0.0;  // [s^-1 / (mmol C m^-3)] Quadratic closure rate for Phy1
double k_benth_graze = 0.0;  // [m/s] Benthic grazing clearance rate (Cloern 1982)
double kmort2_Phy2 = 0.0;  // [s^-1 / (mmol C m^-3)] Quadratic closure rate for Phy2
double kmort_sal = 0.0;    // [s^-1] Maximum salinity mortality rate (default: disabled)
double sal_tol_Phy1 = 5.0; // [PSU] Salinity tolerance for freshwater diatoms
double sal_tol_Phy2 = 15.0; // [PSU] Salinity tolerance for non-diatoms (more euryhaline)

// Calibration target files (strict canonical CSVs; stored under INPUT/Calibration)
char calibration_tidal_targets_file[256] = "INPUT/Calibration/targets_tidal.csv";

// Calibration staging (1=hydrodynamics/transport, 2=sediment+P, 3=eutrophication, 4=carbonate)
int calibration_stage = 1;

// Calibration parameter registry (strict canonical CSV)
char calibration_parameters_file[256] = "INPUT/Calibration/calibration_parameters.csv";

// Unified seasonal targets (strict canonical CSV)
char calibration_wq_seasonal_targets_file[256] = "INPUT/Calibration/targets_seasonal_unified.csv";

/** @brief Spatial settings: DELXI (step size [m]), EL (estuarine length [m]), M
 * (grid points), AMPL (tidal amplitude [m]). */
int DELXI = 0;
int EL = 0;
int M = 0;
double AMPL = 0.0;

/** @brief Tidal component parameters (universal - no hardcoding). */
double M2_fraction = 0.85; // Default M2 component fraction
double S2_fraction = 0.15; // Default S2 component fraction
double M4_ratio = 0.10;    // Default M4/M2 amplitude ratio
double AMPL_CORRECTION =
    1.0; // Dimensionless scaling for tidal elevation (Savenije 2012)
double TIDAL_CFL_TARGET =
    0.95; // Default Courant limit for tidal advection stability (LeVeque 2002)

/** @brief Grid point indices: M1 (last odd), M2 (second to last), M3 (third to
 * last). */
int M1 = 0;
int M2 = 0;
int M3 = 0;

/** @brief Path to riverbed profile file. */
char riverbed_profile_file[256];

// --- Segment Break Indices ---

/** @brief Segment break indices: index_1, index_2, index_3. */
int index_1 = 0;
int index_2 = -1;
int index_3 = -1;
int index_4 = -1;

// --- Segment-Specific Variables ---

/** @brief Segment 1 parameters: B1 (width), LC1 (convergence length), Chezy1
 * (coefficient), Rs1 (storage ratio). */
double B1 = 0.0;
double LC1 = 0.0;
double Chezy1 = 0.0;
double Rs1 = 0.0;

/** @brief Segment 2 parameters: B2, LC2, Chezy2, Rs2. */
double B2 = 0.0;
double LC2 = 0.0;
double Chezy2 = 0.0;
double Rs2 = 0.0;

/** @brief Segment 3 parameters: B3, LC3, Chezy3, Rs3. */
double B3 = 0.0;
double LC3 = 0.0;
double Chezy3 = 0.0;
double Rs3 = 0.0;

// --- Sediment & Hydrodynamics ---

/** @brief Sediment parameters for Segment 1: Mero1 (erosion), tau_ero1 (shear
 * stress erosion), tau_dep1 (shear stress deposition). */
double Mero1 = 0.0, tau_ero1 = 0.0, tau_dep1 = 0.0;

/** @brief Sediment parameters for Segment 2: Mero2, tau_ero2, tau_dep2. */
double Mero2 = 0.0, tau_ero2 = 0.0, tau_dep2 = 0.0;

/** @brief Sediment parameters for Segment 3: Mero3, tau_ero3, tau_dep3. */
double Mero3 = 0.0, tau_ero3 = 0.0, tau_dep3 = 0.0;

/** @brief Sediment parameters for Segment 4 (far upstream): Mero4, tau_ero4, tau_dep4. */
double Mero4 = 0.0, tau_ero4 = 0.0, tau_dep4 = 0.0;

// Per-species tributary mass source rate (mmol/s) accumulated in last step
double trib_mass_rate_last_step[CHEM_COUNT] = {0};
double last_boundary_flux_mouth[CHEM_COUNT] = {0};
double last_boundary_flux_upstream[CHEM_COUNT] = {0};
double last_boundary_flux_mouth_adv_in[CHEM_COUNT] = {0};
double last_boundary_flux_mouth_adv_out[CHEM_COUNT] = {0};
double last_boundary_flux_mouth_disp_in[CHEM_COUNT] = {0};
double last_boundary_flux_mouth_disp_out[CHEM_COUNT] = {0};
double last_boundary_flux_upstream_adv_in[CHEM_COUNT] = {0};
double last_boundary_flux_upstream_adv_out[CHEM_COUNT] = {0};
double last_boundary_flux_upstream_disp_in[CHEM_COUNT] = {0};
double last_boundary_flux_upstream_disp_out[CHEM_COUNT] = {0};

// Biogeochemistry diagnostics clamp counters (per species)
long biogeo_clamp_events[CHEM_COUNT] = {0};

// Array mapping enum indices to variable names
const char *variableNames[CHEM_COUNT] = {
    "Phy1", // Siliceous Phytoplankton
    "Phy2", // Non-Siliceous Phytoplankton
    "Si",   // Silica
    "NO3",  // Nitrate
    "NH4",  // Ammonium
    "PO4",  // Phosphate
    "PIP",  // Particulate Inorganic Phosphorus (Langmuir adsorption)
    "O2",   // Oxygen
    "TOC",  // Total Organic Carbon
    "Sal",  // Salinity
    "SPM",  // Suspended Particulate Matter
    "DIC",  // Dissolved Inorganic Carbon
    "AT",   // Total Alkalinity
    "pCO2", // pCO2
    "PH",   // pH
    "CO2"   // Dissolved CO2
};

// Array mapping forcing enum indices to forcing names
const char *forcingNames[FORCING_COUNT] = {
    "Discharge",   // Upstream river discharge
    "Elevation",   // Tidal elevation
    "Light",       // Light intensity
    "Temperature", // Water temperature
    "Wind"         // Wind speed
};

// Helper functions for data requirements
int isDataInputRequired(Chem var) {
  // Return 0 for variables that are auto-calculated
  return !((1 << var) & AUTO_CALCULATED_VARS);
}

int isTributaryDataRequired(Chem var) {
  // Return 0 for variables that are not needed for tributaries
  return !((1 << var) & TRIBUTARY_EXCLUDED_VARS);
}

// Allocate memory for boundary condition arrays
BCArray upstreamBC[CHEM_COUNT];
BCArray downstreamBC[CHEM_COUNT];

// Unified forcing data structure
ForcingArray forcingData[FORCING_COUNT];

// Tributary structure
Tributary *tributaries = NULL;
int numTributaries = 0;
int tributaryEnabled =
    0; // Definition of tributaryEnabled, initialized to 0 (disabled)

// Cumulative discharge tracking variables
double *Q_cumulative = NULL; // Cumulative discharge at each cell [m³/s]
double *Q_upstream = NULL;   // Upstream discharge at each cell [m³/s]
double *Q_tributaries_total =
    NULL; // Total tributary discharge at each cell [m³/s]
double *mixing_efficiency =
    NULL;                       // Tributary mixing efficiency at each cell [-]
double *momentum_factor = NULL; // Tributary momentum enhancement factor [-]

// Biogeochemical Variables
double K;
double NPP_NO3[MAXM + 1][2];
double NPP_NO3_tot[MAXM + 1];
double NPP_NH4[MAXM + 1][2];
double NPP_NH4_tot[MAXM + 1];
double GPP[MAXM + 1][2];
double phydeath[MAXM + 1][2];
double phydeath_tot[MAXM + 1];
double NPP[MAXM + 1];
double Si_consumption[MAXM + 1];
double adegrad[MAXM + 1];
double denit[MAXM + 1];
double nitrif[MAXM + 1];
double o2air[MAXM + 1];
double erosion_s[MAXM + 1];
double deposition_s[MAXM + 1];
double erosion_v[MAXM + 1];
double deposition_v[MAXM + 1];
double tau_b[MAXM + 1];
// Removed unused globals: kwind[MAXM], kflow[MAXM], vp[MAXM] - shadowed by locals
double co2air[MAXM + 1];
double reactionDIC[MAXM + 1];
double reactionTA[MAXM + 1];
// Removed unused: reactionHS[MAXM] - never referenced
double hco3[MAXM + 1];
double co3[MAXM + 1];
double pCO2atmo;
double piston_velocity_scale = 1.0;

// Diagnostic arrays for nutrient/light limitation and KD (saved at output stride)
double diag_fN[MAXM + 1];
double diag_fP[MAXM + 1];
double diag_fSi[MAXM + 1];
double diag_fI[MAXM + 1];
double diag_KD_total[MAXM + 1];

// Phase 4: override flags for params-driven biogeochemical constants.
// These are reset at the start of read_parameters() so repeated calls are deterministic.
int biogeo_override_Pb[2] = {0, 0};
int biogeo_override_alpha[2] = {0, 0};
int biogeo_override_kmortality[2] = {0, 0};
int biogeo_override_kox = 0;
int biogeo_override_knit = 0;
int biogeo_override_pCO2atmo = 0;
int biogeo_override_KN[2] = {0, 0};
int biogeo_override_KPO4[2] = {0, 0};
int biogeo_override_KSi = 0;

// Phase 5: carbonate solver failure policy (default preserves current behavior).
CarbonateFailurePolicy carbonate_solver_failure_policy = CARBONATE_FAIL_HOLD_LAST;

// Model parameters
// Constants
double alpha[2];
double Pb[2];
double kmax[2];
double kexcr[2];
double kgrowth[2];
double kmaint[2];
double kmortality[2];
double KSi[2];
double KN[2];
double KPO4[2];
double KTOC;
double KO2_ox;
double KO2_nit;
double KNO3;
double KNH4;
double KinO2;
double redsi;
double ksi_diss = 0.5;           // Fraction of biogenic Si from dead diatoms that dissolves [0-1]
double redn;
double redp;
double kox;
double kdenit;
double knit;
double kcbod_fast = 0.0;
double KO2_cbod = 31.0;
double KTOC_fast = 80.0;
double k_ww_bod = 0.0;           // Direct wastewater BOD O2 demand [mmol O2/m³/s], default 0 (off)
double KO2_ww_bod = 15.0;        // Half-saturation O2 for wastewater BOD [mmol/m³]
double SOD_rate = 0.0;           // Sediment O2 demand [mmol O2/m²/s], default 0 (off)
double KO2_SOD = 31.0;           // Half-saturation O2 for SOD [mmol/m³]
double kbg;
double kbg_fresh = 0.0;      // Freshwater excess attenuation [m⁻¹] (0 = disabled)
double S_cdom_threshold = 2.0; // Salinity midpoint for CDOM transition [PSU]
double S_cdom_width = 1.0;   // Width of salinity-CDOM transition [PSU]
double kbg_upstream = 0.0;   // Excess upstream background attenuation [m⁻¹] (0 = disabled)
double kbg_transition_km = 0.0; // Distance from mouth for turbidity transition [km] (0 = disabled)
double kbg_blend_km = 5.0;  // Spatial half-width of transition zone [km]
double kspm;
double kCDOM = 0.0;             // CDOM-specific light attenuation [m⁻¹/(mgC/L)]
double k_turb_Q = 0.0;      // Discharge-driven turbidity coefficient [m⁻¹]
double Q_turb_ref = 0.0;     // Reference discharge for turbidity [m³/s] (0 = auto-compute from data)
double alpha_turb = 1.5;    // Power-law exponent for Q-KD relationship [-]
double k_resus = 0.0;       // Depth-dependent resuspension turbidity [m⁻¹·m] (0 = disabled)
double sal_half_phy = -1.0;  // Salinity half-sat for freshwater Phy growth inhibition [PSU] (<0 = disabled)
double CN_toc = 6.625;           // Effective C:N ratio for bulk TOC mineralization (Redfield default)
double ws = 5.0e-4;  // Default 0.5 mm/s settling velocity (Winterwerp 2002), read from params.txt
double ws_phy = 0.1 / 86400.0; // Phyto settling ~0.1 m/day (RIVE/Literature default)
double rho_w;
double g;

// P-adsorption parameters (Langmuir isotherm, Nguyen et al. 2019)
double P_ac = 1.0;   // Adsorption capacity [mmol P / kg SPM], default from literature
double K_ps = 2.0;   // Half-saturation constant [mmol P / m³], default from literature
double k_ads = 1.0 / 3600.0;  // Adsorption rate constant [s⁻¹], ~1 hour equilibration

// Transport variables
// Biogeochemical tuning parameters
// KD_Phy = 0.005 m^2/mmolC; with C:Chl=15 this is approx 0.006 m^2/mgChl
double KD_Phy = 0.005;       // Default from hardcoded value (Phase 4)
double K_NH4_switch = 5.0;   // Default from hardcoded value (Phase 4)
double C_Chl_ratio = 15.0;   // Carbon to Chlorophyll ratio (g C / g Chl); 15 for shade-adapted tropical estuary (Geider 1997, Cloern 1995)

// Mechanistic CO2 limitation for phytoplankton primary production
int enable_co2_limitation = 0; // Default OFF (backwards compatible)
double K_CO2_Phy = 10.0;       // [mmolC/m^3] half-saturation for CO2 limitation

double fl[MAXM + 1];
double disp[MAXM + 1];
double watflux[MAXM + 1];

// Hydrodynamic variables
// Note: Y and E arrays are defined in hydrodynamics.c to avoid conflicts
double totalArea[MAXM + 1];
double baseArea[MAXM + 1];
double freeArea[MAXM + 1];
double totalAreaOld[MAXM + 1];
double riverbed_depth[MAXM + 1];
double elevation[MAXM + 1];
double velocity[MAXM + 1];
double width[MAXM + 1];
double tempFreeArea[MAXM + 1];
double tempVelocity[MAXM + 1];
double waterDepth[MAXM + 1];
double level[MAXM + 1]; // Water level [m]

// Hydrodynamic coefficients and properties
double Chezy[MAXM + 1];   // Chezy coefficient [m^(1/2)/s]
double FRIC[MAXM + 1];    // Friction coefficient [1/Chezy²] [-]
double Mero[MAXM + 1];    // Erosion coefficient [mg m-2 s-1]
double tau_ero[MAXM + 1]; // Critical shear stress for erosion [N m-2]
double tau_dep[MAXM + 1]; // Critical shear stress for deposition [N m-2]
double rs[MAXM + 1];      // Storage width ratio [-]

double C[MAXM + 1][5];
double Z[MAXM + 1];

// Tributary boundary data arrays
// Legacy trib_data/trib_time arrays removed — tributary data is now
// dynamically allocated via the Tributary struct in readTributaryData().

// Define the global variable structure
struct Verb v[MAXV];

/**
 * @brief Allocate memory for cumulative discharge tracking arrays
 */
void allocateCumulativeDischargeMemory() {
  if (Q_cumulative == NULL) {
    Q_cumulative = (double *)calloc(MAXM + 1, sizeof(double));
    Q_upstream = (double *)calloc(MAXM + 1, sizeof(double));
    Q_tributaries_total = (double *)calloc(MAXM + 1, sizeof(double));
    mixing_efficiency = (double *)calloc(MAXM + 1, sizeof(double));
    momentum_factor = (double *)calloc(MAXM + 1, sizeof(double));

    if (!Q_cumulative || !Q_upstream || !Q_tributaries_total ||
        !mixing_efficiency || !momentum_factor) {
      printf("❌ Error: Failed to allocate memory for cumulative discharge "
             "arrays\n");
      exit(1);
    }

    // Initialize with default values
    for (int i = 0; i <= MAXM; i++) {
      Q_cumulative[i] = 0.0;
      Q_upstream[i] = 0.0;
      Q_tributaries_total[i] = 0.0;
      mixing_efficiency[i] = 1.0; // Default complete mixing
      momentum_factor[i] = 1.0;   // Default no momentum enhancement
    }

    // Memory allocated successfully - silent operation
  }
}

/**
 * @brief Free memory for cumulative discharge tracking arrays
 */
void freeCumulativeDischargeMemory() {
  if (Q_cumulative) {
    free(Q_cumulative);
    Q_cumulative = NULL;
  }
  if (Q_upstream) {
    free(Q_upstream);
    Q_upstream = NULL;
  }
  if (Q_tributaries_total) {
    free(Q_tributaries_total);
    Q_tributaries_total = NULL;
  }
  if (mixing_efficiency) {
    free(mixing_efficiency);
    mixing_efficiency = NULL;
  }
  if (momentum_factor) {
    free(momentum_factor);
    momentum_factor = NULL;
  }
  // Memory freed successfully - silent operation
}
