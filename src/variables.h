/*____________________________________*/
/* variables.h                       */
/* Declare global variables, functions */
/* Last modified: 2/2025 An Nguyen    */
/*____________________________________*/

#ifndef VARIABLES_H
#define VARIABLES_H
#include "define.h"
#include <stdbool.h> // Add this for bool type
#include <stddef.h>  // Add this for offsetof macro


// Output storage format for simulation files
typedef enum {
  OUTPUT_FORMAT_CSV = 0,
  OUTPUT_FORMAT_BINARY = 1
} OutputFormat;

// Enum for chemical variables
typedef enum {
  Phy1,
  Phy2,
  Si,
  NO3,
  NH4,
  PO4,
  PIP,  // Particulate Inorganic Phosphorus (Langmuir adsorption onto SPM)
  O2,
  TOC,
  Sal,
  SPM,
  DIC,
  AT,
  pCO2, // auto-calculated
  PH,   // auto-calculated
  CO2,  // auto-calculated
  CHEM_COUNT
} Chem;

// MAXV must equal CHEM_COUNT so v[] fits all species
#define MAXV CHEM_COUNT

// Enum for forcing data types
typedef enum {
  FORCING_DISCHARGE,   // Upstream discharge
  FORCING_ELEVATION,   // Tidal elevation at downstream boundary
  FORCING_LIGHT,       // Light intensity
  FORCING_TEMPERATURE, // Water temperature
  FORCING_WIND,        // Wind speed
  FORCING_COUNT
} Forcing;

// Constants for data requirements, which are not needed for input
#define AUTO_CALCULATED_VARS (1 << pCO2 | 1 << PH | 1 << CO2 | 1 << PIP)
// SPM excluded: TSS mainly comes from the river (upstream), not from other tributaries
#define TRIBUTARY_EXCLUDED_VARS (1 << Sal | 1 << SPM | 1 << pCO2 | 1 << PH | 1 << CO2)

// Use CHEM_COUNT instead of NUM_BC_VARS for consistency
#define NUM_BC_VARS CHEM_COUNT

// Function to check data input requirements
int isDataInputRequired(Chem var);
int isTributaryDataRequired(Chem var);

// Global simulation parameters
extern int calibration_mode; // 0: none, 1: BOBYQA, 2: Nelder-Mead,
                             // 3: SBPLX, 4: COBYLA, 5: DIRECT-L,
                             // 6: CRS2-LM, 7: ISRES
extern int debug_level;      // Debug level (0=off, 1=error, 2=warn, 3=info,
                             // 4=detail, 5=trace)
extern OutputFormat output_storage_format; // File storage mode for model output (CSV/BINARY)
extern int enable_flux_output;   // 1=write flux diagnostics, 0=skip
extern int enable_reaction_output; // 1=write biogeochemical rates, 0=skip


// Geometry alignment controls
extern int
    align_segment_break_with_major_tributary; // 1=auto-align segment transition
                                              // to dominant tributary cell

// Dispersion calibration parameters (params.txt configuration)
extern double C_VDB;         // Van der Burgh shape factor [3.5-7.0]
extern double D0_CORRECTION; // D0 magnitude correction factor [0.5-2.0]
extern double K_DISCHARGE_SENSITIVITY;  // Seasonal K correction via Canter-Cremers: K_eff = K × (1 + α×(N/N_ref−1)) [0-1]

// Calibration controller (configuration-driven via params.txt)
extern double calibration_warmup_days;
extern double calibration_diagnostic_days;
extern int calibration_max_iterations;


// Time and spatial settings
extern long MAXT;
extern long WARMUP;
extern int DELTI, TS, DELXI, EL, M, M1, M2;
extern long simulation_window_seconds; // New global variable for simulation
                                       // window in seconds
extern double AMPL;

// Tidal component parameters (universal)
extern double M2_fraction; // M2 component fraction (0.7-0.9)
extern double S2_fraction; // S2 component fraction (0.1-0.3)
extern double M4_ratio;    // M4/M2 amplitude ratio (0.05-0.20)
extern double
    AMPL_CORRECTION;            // Dimensionless scaling for tidal elevation
                                // (Savenije 2012 calibration)
extern double TIDAL_CFL_TARGET; // Courant stability target for tidal advection
                                // (dimensionless)

// Mass balance validation mode
// Biogeochemical reaction modes
#define BIOGEO_OFF 0
#define BIOGEO_FULL 1
#define BIOGEO_P_ADSORPTION_ONLY 2

extern int
    enable_biogeochemical_reactions;

extern int enable_carbonate_diagnostics;
extern int enable_suspended_sediment_dynamics;

extern int
    warmup_species_minimal; // 1=Sal+SPM only during warmup, 0=all species
extern int enable_upstream_discharge_smoothing; // 1=apply exponential smoothing
                                                // to upstream NET discharge
extern double
    upstream_discharge_smoothing_timescale_days; // e-folding response time for
                                                 // smoothing [days]
extern double upstream_discharge_scale; // multiplier applied to upstream NET discharge (dimensionless, >0)
extern int enable_tributary_subdaily_pulses; // 1=apply mean-preserving subdaily modulation to tributary discharge
extern double tributary_pulse_shape;         // dimensionless pulse intensity "a" for exp(a*sin(theta))/I0(a)
extern double urban_tributary_pulse_boost;   // multiplier for urban/lateral tributaries (as defined in config)
extern double spm_upstream_scale; // multiplier applied to upstream SPM (dimensionless, >0)
extern TributaryMassInjectionMode
  tributary_mass_injection_mode; // Controls tributary mass source formulation

// Quadratic mortality closure for implicit zooplankton grazing
// Fennel et al. (2006, JGR-Oceans), Edwards & Yool (2000)
extern double kmort2_Phy1;  // [s^-1 / (mmol C m^-3)] Quadratic closure Phy1
extern double kmort2_Phy2;  // [s^-1 / (mmol C m^-3)] Quadratic closure Phy2

// Benthic grazing — depth-dependent phytoplankton removal
// Cloern (1982, Ecol. Monogr.); Alpine & Cloern (1992, Mar. Ecol. Prog. Ser.)
// Loss rate = k_benth_graze / H [s^-1]: shallow sections lose more biomass.
extern double k_benth_graze;  // [m/s] Benthic filtration clearance rate

// Salinity stress mortality for freshwater phytoplankton (Garnier et al., RIVE model)
// Freshwater diatoms (Aulacoseira, Melosira) cannot survive above ~5-10 PSU.
// f(S) = S^2 / (S^2 + sal_tol^2): sigmoid 0→1 as salinity increases.
extern double kmort_sal;       // [s^-1] Maximum salinity mortality rate
extern double sal_tol_Phy1;    // [PSU] Salinity tolerance threshold for Phy1 (freshwater diatoms)
extern double sal_tol_Phy2;    // [PSU] Salinity tolerance threshold for Phy2 (more euryhaline)

// Calibration target files (strict canonical CSVs; stored under INPUT/Calibration)
extern char calibration_tidal_targets_file[256];

// Calibration staging (1=hydrodynamics/transport, 2=sediment+P, 3=eutrophication, 4=carbonate)
extern int calibration_stage;

// Calibration parameter registry (strict canonical CSV; stored under INPUT/Calibration)
extern char calibration_parameters_file[256];

// Unified seasonal WQ targets (strict canonical CSV; stored under INPUT/Calibration)
extern char calibration_wq_seasonal_targets_file[256];

// Transport and calibration parameters
extern int num_segments;
extern int segment_transition_width; // Transition half-width in cells [cells]
extern double dispersion_abs_min;    // Absolute floor for dispersion [m2/s]
extern double dispersion_abs_max;    // Absolute ceiling for dispersion [m2/s]

// Biogeochemical tuning parameters
extern double KD_Phy;        // Light attenuation by phytoplankton [m^2/mmol]
extern double K_NH4_switch;  // NH4 preference half-saturation [mmol/m^3]
extern double C_Chl_ratio;   // Carbon to Chlorophyll ratio (g C / g Chl)

// Mechanistic CO2 limitation for primary production
extern int enable_co2_limitation; // 1=apply CO2(aq) limitation to GPP/NPP, 0=off
extern double K_CO2_Phy;          // CO2 half-saturation [mmolC/m^3]

// === Tidal Physics Parameters (Astronomical Constants) ===
extern double TIDAL_PERIOD_M2_HOURS; // M2 tidal period (hours)
extern double TIDAL_PERIOD_S2_HOURS; // S2 tidal period (hours)
extern double TIDAL_PHASE_M2_RAD;    // M2 phase (radians)
extern double TIDAL_PHASE_S2_RAD;    // S2 phase (radians)

// Segment 1 parameters
extern int index_1;
extern double B1, LC1;
extern double Chezy1;
extern double Rs1;

// Segment 2 parameters
extern int index_2;
extern double B2, LC2;
extern double Chezy2;
extern double Rs2;

// Segment 3 parameters
extern int index_3;
extern double B3, LC3;
extern double Chezy3;
extern double Rs3;

// Segment 4 (far upstream SPM) break index
extern int index_4;

// Sediment transport parameters (erosion & deposition)
extern double Mero1, tau_ero1, tau_dep1; // Segment 1
extern double Mero2, tau_ero2, tau_dep2; // Segment 2
extern double Mero3, tau_ero3, tau_dep3; // Segment 3
extern double Mero4, tau_ero4, tau_dep4; // Segment 4 (far upstream)

// Riverbed profile file path
extern char riverbed_profile_file[256];

extern double trib_mass_rate_last_step[CHEM_COUNT];
extern double last_boundary_flux_mouth[CHEM_COUNT];
extern double last_boundary_flux_upstream[CHEM_COUNT];
extern double last_boundary_flux_mouth_adv_in[CHEM_COUNT];
extern double last_boundary_flux_mouth_adv_out[CHEM_COUNT];
extern double last_boundary_flux_mouth_disp_in[CHEM_COUNT];
extern double last_boundary_flux_mouth_disp_out[CHEM_COUNT];
extern double last_boundary_flux_upstream_adv_in[CHEM_COUNT];
extern double last_boundary_flux_upstream_adv_out[CHEM_COUNT];
extern double last_boundary_flux_upstream_disp_in[CHEM_COUNT];
extern double last_boundary_flux_upstream_disp_out[CHEM_COUNT];

// Biogeochemistry diagnostics: count of non-negativity clamps applied per species.
// These counters are intended for audits/tests to detect when truncation is used.
extern long biogeo_clamp_events[CHEM_COUNT];

/*==============================*/
/* Boundary Condition Variables */
/*==============================*/

// Define a unified boundary data structure indexed by Chem
typedef struct {
  double *time;
  double *data;
  const char *filePath;
  int required;
  int dataSize; // Added field to track actual data size for safe interpolation
} BCArray;

// For upstream/downstream
extern BCArray upstreamBC[CHEM_COUNT];
extern BCArray downstreamBC[CHEM_COUNT];

/*==============================*/
/* Forcing  Inputs              */
/*==============================*/

// Define a unified forcing data structure indexed by Forcing enum
typedef struct {
  double *time;
  double *data;
  const char *filePath;
  int required;
  int dataSize; // Added field to track actual data size
} ForcingArray;

// Central array for all forcing data
extern ForcingArray forcingData[FORCING_COUNT];

// Names of forcing variables for easier access
extern const char *forcingNames[FORCING_COUNT];

/*==============================*/
/* Tributary Structure          */
/*==============================*/

typedef struct {
  char name[50];
  int cellIndex;

  // File path for discharge (always required)
  char dischargeFile[200]; // Use 200 length for consistency
  double *dischargeTime, *discharge;
  int dischargeDataSize; // Added to store the size of discharge data arrays

  // Chemical data - dynamically addressed by enum Chem values
  struct {
    char filePath[200]; // Use 200
    double *timeArray;
    double *dataArray;
    int dataSize; // Added to store size of individual chemical data arrays
  } chemicalData[CHEM_COUNT];

  // Current concentrations (used during simulation)
  double concentration[CHEM_COUNT];
  
  // Flag indicating if the tributary is an urban/wastewater source
  int is_urban;
} Tributary;

extern Tributary *tributaries; // Changed to pointer for dynamic allocation
extern int numTributaries;
extern int tributaryEnabled;

/*==============================*/
/* Cumulative Discharge Tracking */
/*==============================*/
extern double *Q_cumulative; // Cumulative discharge at each cell [m³/s]
extern double *Q_upstream;   // Upstream discharge at each cell [m³/s]
extern double
    *Q_tributaries_total; // Total tributary discharge at each cell [m³/s]
extern double
    *mixing_efficiency;         // Tributary mixing efficiency at each cell [-]
extern double *momentum_factor; // Tributary momentum enhancement factor [-]

/*==============================*/
/* Boundary Condition Structure */
/*==============================*/
typedef struct {
  const char *name;      // Displayed name of the variable
  const char *configKey; // Key to look for in config file
  const char **filePath; // Pointer to the file path variable
  double **timeArray;    // Pointer to the time array pointer
  double **valueArray;   // Pointer to the value array pointer
  int isRequired;        // Whether this parameter is required
} BoundaryCondition;

// External references to BoundaryCondition arrays from file.c
extern BoundaryCondition upstreamBCs[];
extern BoundaryCondition downstreamBCs[];
extern int numUpstreamBCs;
extern int numDownstreamBCs;

/*==============================*/
/* Function Declarations        */
/*==============================*/

// Declaration only (note the extern keyword)
extern const char *variableNames[CHEM_COUNT]; // Use CHEM_COUNT for consistency

/**
 * @brief Returns the name of a chemical species based on its enum index
 * @param s Species index from enum Chem
 * @return The name of the species as a string
 */
const char *getSpeciesName(int s);

/* init.c */
extern void Init();
extern void readConfigFile(const char *filename);
extern void allocateBoundaryMemory();
extern void freeBoundaryMemory();
extern void readBoundaryData(const char *filename);
extern void readTributaryConfigFile(const char *filename);
extern void allocateTributaryMemory();
extern void readGlobalConfigSettings(const char *filename);
extern void readTributaryData(const char *filename);
extern void allocateCumulativeDischargeMemory();
extern void freeCumulativeDischargeMemory();
extern void updateCumulativeDischarge(int t);
extern double calculateTributaryMixingEfficiency(int cell, double Q_trib,
                                                 double Q_main);
extern double calculateMomentumFactor(int cell, double Q_trib, double Q_main,
                                      double velocity_main);

/* file.c */
extern int Hydwrite(int);
extern void Rates(double *co, const char *s, int t);
extern void Fluxwrite(int species, int t);
extern int readFile(int datamax, double *gg, double *ff, const char *s);
extern int fileExists(const char *path);
extern int isValidPath(const char *path);
extern void checkFilePath(const char *filePath, const char *varName);
extern void readRiverbedProfile(double *slope, const char *filename);
extern void read_parameters(const char *filename);

extern void Hyd(int);
extern void export_hydro_iterations_summary(const char *path);
extern void calculateTidalDynamics(int); // Coupled analytical tidal dynamics (Cai et al., 2016)
extern void finalize_remaining_tidal_range(void);

/* transport.c */
extern void Transport(int);

/* biogeo.c */
extern void Biogeo(int);
extern double interpolateAtTime(int t, double *timeSeries, double *valueSeries);
extern double calc_co2_solubility(double Tk, double salinity);
extern double calc_co2_dry_wet_conversion(double temp_air, double pCO2atm, double salinity);
extern double calc_water_density(double temp_water, double salinity);
extern double compute_k600(double U_water, double U_wind, double h, double TSS, double temperature, double Sc_gas);
extern double disK1(double Sal, double t);
extern double disK2(double Sal, double t);
extern double disK3(double Sal, double t);
extern double disK4(double Sal, double t);
extern double disK5(double Sal, double t);

/* bcforcing.c */
extern void bgboundary(int);
extern double Tide(int);
extern double Discharge_ups(int);
extern double Discharge(int, int, double *);
extern double calculatePhysicsBasedTributaryDischarge(int, int);
extern double windspeed(int, int);

/* biogeout.c */
extern double waterT(int);
extern double I0(int);
extern double Pbmax(int, int);
extern double kmort(int, int);
extern double kmaintenance(int, int);
extern double Fhetox(int);
extern double Fhetden(int);
extern double Fnit(int);
extern double O2sat(int, int);

/*==============================*/
/* Hydrodynamic Variables       */
/*==============================*/

extern double Y[MAXM + 1]; // Convergence test array for velocities
extern double E[MAXM + 1]; // Convergence test array for water levels

extern double Chezy[MAXM + 1];   // Chezy coefficient [m^(1/2)/s]
extern double FRIC[MAXM + 1];    // Friction coefficient [1/Chezy²] [-]
extern double Mero[MAXM + 1];    // Erosion coefficient [mg m-2 s-1]
extern double tau_ero[MAXM + 1]; // Critical shear stress for erosion [N m-2]
extern double tau_dep[MAXM + 1]; // Critical shear stress for deposition [N m-2]
extern double rs[MAXM + 1];      // Storage width ratio [-]

extern double totalArea[MAXM + 1];      // Total cross-sectional area [m²]
extern double baseArea[MAXM + 1];       // Base cross-sectional area [m²]
extern double freeArea[MAXM + 1];       // Free surface cross-sectional area [m²]
extern double totalAreaOld[MAXM + 1];   // Previous total cross-sectional area [m²]
extern double riverbed_depth[MAXM + 1]; // Riverbed depth below datum [m]
extern double elevation[MAXM + 1];      // Bed elevation [m]
extern double velocity[MAXM + 1];       // Flow velocity [m/s]
extern double width[MAXM + 1];          // Channel width [m]
extern double tempFreeArea[MAXM + 1];   // Temporary free surface area [m²]
extern double tempVelocity[MAXM + 1];   // Temporary velocity [m/s]
extern double waterDepth[MAXM + 1];     // Water depth (totalArea/width) [m]
extern double level[MAXM + 1]; // Water level [m] - for tidal volume calculations
extern double tidal_excursion[MAXM + 1];    // Tidal excursion E(x) [m]
extern double tidal_velocity_amp[MAXM + 1]; // Velocity amplitude v(x) [m/s]
extern double phase_lag[MAXM + 1];          // Phase lag ε(x) [rad]
extern double damping_factor[MAXM + 1];     // Damping factor δ(x) [-]

extern double C[MAXM + 1][5]; // Coefficients for tridiagonal matrix
extern double Z[MAXM + 1];    // Right-hand side for tridiagonal matrix

// Historical tracking arrays for validation
extern double velocity_min[MAXM + 1];         // Minimum velocity tracking
extern double velocity_max[MAXM + 1];         // Maximum velocity tracking
extern double dispersion_min[MAXM + 1];       // Minimum dispersion tracking
extern double dispersion_max[MAXM + 1];       // Maximum dispersion tracking
extern double daily_min_level[MAXM + 1];      // Daily minimum water level
extern double daily_max_level[MAXM + 1];      // Daily maximum water level
extern double tidal_range_mean[MAXM + 1];     // Mean tidal range for each cell [m]

/*==============================*/
/* Transport & Biogechem       */
/*==============================*/

#ifndef STRUCT_VERB
#define STRUCT_VERB
struct Verb {
  char name[20];
  int env;
  double c[MAXM + 1];
  double clb;
  double cub;
  double avg[MAXM + 1];
  double concflux[MAXM + 1];
  double advflux[MAXM + 1];
  double disflux[MAXM + 1];
};
extern struct Verb v[MAXV];
#endif

/*==============================*/
/* Biogeochemical Rate Constants */
/*==============================*/

extern double Pb[2];         // Maximum specific photosynthetic rate
extern double alpha[2];      // Photosynthetic efficiency
extern double kexcr[2];      // Excretion rate
extern double kgrowth[2];    // Growth rate
extern double kmaint[2];     // Maintenance rate
extern double kmortality[2]; // Mortality rate
extern double KSi[2];        // Silica half-saturation constant
extern double KN[2];         // Nitrate half-saturation constant
extern double KPO4[2];       // Phosphate half-saturation constant
extern double KTOC;          // Total organic carbon half-saturation constant
extern double KO2_ox;        // Oxygen oxidation half-saturation constant
extern double KO2_nit;       // Oxygen nitrification half-saturation constant
extern double KinO2;         // Oxygen inhibition constant for denitrification
extern double KNO3;          // Nitrate half-saturation constant
extern double KNH4;          // Ammonium half-saturation constant
extern double redsi;         // Redfield silica ratio
extern double ksi_diss;      // Fraction of biogenic Si from dead diatoms dissolved [0-1]
extern double redn;          // Redfield nitrogen ratio
extern double redp;          // Redfield phosphorus ratio
extern double kox;           // Aerobic degradation rate
extern double kdenit;        // Denitrification rate
extern double knit;          // Nitrification rate
extern double kcbod_fast;    // Fast labile-C oxygen demand rate [s^-1]
extern double KO2_cbod;      // O2 half-saturation for fast CBOD [mmol/m3]
extern double KTOC_fast;     // TOC half-saturation for fast CBOD [mmolC/m3]
extern double k_ww_bod;      // Direct wastewater BOD O2 demand rate [mmol O2/m³/s]
extern double KO2_ww_bod;    // O2 half-saturation for wastewater BOD [mmol/m³]
extern double SOD_rate;      // Sediment O2 demand flux [mmol O2/m²/s]
extern double KO2_SOD;       // Half-saturation O2 for SOD [mmol/m³]
extern double kbg;           // Background attenuation coefficient
extern double kbg_fresh;     // Freshwater excess background attenuation [m⁻¹] (CDOM/humic)
extern double S_cdom_threshold; // Salinity midpoint for CDOM transition [PSU]
extern double S_cdom_width;  // Width of salinity-CDOM transition [PSU]
extern double kbg_upstream;  // Excess upstream background attenuation [m⁻¹] (terrigenous CDOM)
extern double kbg_transition_km; // Distance from mouth for turbidity transition [km]
extern double kbg_blend_km;  // Spatial half-width of transition zone [km]
extern double kspm; // Suspended particulate matter attenuation coefficient
extern double kCDOM;
extern double k_turb_Q;
extern double Q_turb_ref;    // Reference discharge for turbidity scaling [m³/s]
extern double alpha_turb;    // Power-law exponent for discharge-turbidity relationship [-]
extern double k_resus;       // Depth-dependent resuspension turbidity [m⁻¹·m] (0 = disabled)
extern double sal_half_phy;   // Salinity half-saturation for freshwater Phy inhibition [PSU]
extern double CN_toc;        // Effective C:N ratio for bulk TOC mineralization
extern double ws;    // Settling velocity (SPM)
extern double ws_phy; // Settling velocity (Phytoplankton)
extern double rho_w; // Water density
extern double g;     // Gravity

// P-adsorption parameters (Langmuir isotherm, Nguyen et al. 2019)
extern double P_ac;   // Adsorption capacity [mmol P / kg SPM]
extern double K_ps;   // Half-saturation constant for P adsorption [mmol P / m³]
extern double k_ads;  // Adsorption rate constant [s⁻¹]

extern double K;
extern double NPP_NO3[MAXM + 1][2];
extern double NPP_NO3_tot[MAXM + 1];
extern double NPP_NH4[MAXM + 1][2];
extern double NPP_NH4_tot[MAXM + 1];
extern double GPP[MAXM + 1][2];
extern double phydeath[MAXM + 1][2];
extern double phydeath_tot[MAXM + 1];
extern double NPP[MAXM + 1];
extern double Si_consumption[MAXM + 1];
extern double adegrad[MAXM + 1];
extern double denit[MAXM + 1];
extern double nitrif[MAXM + 1];
extern double o2air[MAXM + 1];
extern double erosion_s[MAXM + 1];
extern double deposition_s[MAXM + 1];
extern double erosion_v[MAXM + 1];
extern double deposition_v[MAXM + 1];
extern double tau_b[MAXM + 1];
extern double co2air[MAXM + 1];
extern double reactionDIC[MAXM + 1];
extern double reactionTA[MAXM + 1];
extern double hco3[MAXM + 1];
extern double co3[MAXM + 1];
extern double pCO2atmo;

// Diagnostic arrays for nutrient/light limitation and KD
extern double diag_fN[MAXM + 1];
extern double diag_fP[MAXM + 1];
extern double diag_fSi[MAXM + 1];
extern double diag_fI[MAXM + 1];
extern double diag_KD_total[MAXM + 1];

// Air–water gas exchange tuning (dimensionless). Default 1.0.
// Used by biogeo.c piston_velocity(); also set by calibration objective 3/4.
extern double piston_velocity_scale;

// ---------------------------------------------------------------------------
// Phase 4: params.txt overrides for key biogeochemical/carbonate constants
// ---------------------------------------------------------------------------
extern int biogeo_override_Pb[2];
extern int biogeo_override_alpha[2];
extern int biogeo_override_kmortality[2];
extern int biogeo_override_kox;
extern int biogeo_override_knit;
extern int biogeo_override_pCO2atmo;
extern int biogeo_override_KN[2];
extern int biogeo_override_KPO4[2];
extern int biogeo_override_KSi;

// ---------------------------------------------------------------------------
// Phase 5: Carbonate solver failure policy
// ---------------------------------------------------------------------------
typedef enum {
  CARBONATE_FAIL_HOLD_LAST = 0,         // Keep last diagnostics and proceed
  CARBONATE_FAIL_SKIP_CO2_EXCHANGE = 1, // Proceed, but set CO2 exchange to 0 when solver fails
  CARBONATE_FAIL_FATAL = 2              // Abort immediately on solver failure
} CarbonateFailurePolicy;

extern CarbonateFailurePolicy carbonate_solver_failure_policy;

// Transport variables
extern double fl[MAXM + 1];
extern double disp[MAXM + 1];
extern double watflux[MAXM + 1];

double Min(double a, double b);
double Max(double a, double b);

#endif // VARIABLES_H
