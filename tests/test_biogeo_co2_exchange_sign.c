/**
 * @file test_biogeo_co2_exchange_sign.c
 * @brief Test that CO2 air-water exchange has correct sign convention
 * 
 * This test verifies Phase 5 requirement: CO2 exchange should follow
 * physically correct sign convention:
 * - Positive flux = degassing (water to air) when pCO2_water > pCO2_atm
 * - Negative flux = invasion (air to water) when pCO2_water < pCO2_atm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "variables.h"

extern void Biogeo(int t);
extern void assignBiogeochemicalRateConstants(void);

int main(void)
{
    debug_level = 0;

    // 1-cell domain
    M = 1;
    M1 = 0;
    M2 = 0;
    M3 = 0;

    DELXI = 1000;
    EL = 1000;
    WARMUP = 0;
    TS = 1;
    DELTI = 3600; // 1 hour timestep

    enable_biogeochemical_reactions = 1;
    enable_carbonate_diagnostics = 0;

    assignBiogeochemicalRateConstants();
    
    // Set atmospheric CO2 to 400 µatm (typical current level)
    pCO2atmo = 400.0;

    // Provide forcing data
    static double temp_time[2] = {0.0, 10.0};
    static double temp_val[2] = {28.0, 28.0};
    forcingData[FORCING_TEMPERATURE].time = temp_time;
    forcingData[FORCING_TEMPERATURE].data = temp_val;
    forcingData[FORCING_TEMPERATURE].dataSize = 2;

    static double light_time[2] = {0.0, 10.0};
    static double light_val[2] = {100.0, 100.0}; // Low light to minimize NPP interference
    forcingData[FORCING_LIGHT].time = light_time;
    forcingData[FORCING_LIGHT].data = light_val;
    forcingData[FORCING_LIGHT].dataSize = 2;

    static double wind_time[2] = {0.0, 10.0};
    static double wind_val[2] = {5.0, 5.0}; // Moderate wind for gas exchange
    forcingData[FORCING_WIND].time = wind_time;
    forcingData[FORCING_WIND].data = wind_val;
    forcingData[FORCING_WIND].dataSize = 2;

    waterDepth[1] = 10.0;
    width[1] = 100.0;
    Chezy[1] = 50.0;
    velocity[1] = 0.5;

    // Test scenario 1: Supersaturated water (pCO2_water > pCO2_atm)
    // Should result in degassing (positive CO2 loss from water)
    printf("Test 1: Supersaturated water (degassing expected)\n");
    
    v[Sal].c[1] = 5.0;
    v[Phy1].c[1] = 5.0;
    v[Phy2].c[1] = 5.0;
    v[Si].c[1] = 50.0;
    v[PO4].c[1] = 5.0;
    v[NO3].c[1] = 50.0;
    v[NH4].c[1] = 10.0;
    v[TOC].c[1] = 300.0;
    v[O2].c[1] = 160.0;
    v[DIC].c[1] = 3000.0; // High DIC
    v[AT].c[1] = 2500.0;  // Lower AT/DIC ratio = lower pH = higher CO2
    v[SPM].c[1] = 10.0;

    double DIC_initial_1 = v[DIC].c[1];
    
    Biogeo(0);

    double DIC_change_1 = v[DIC].c[1] - DIC_initial_1;
    
    // In supersaturated conditions, CO2 should escape (DIC should decrease)
    if (DIC_change_1 >= 0.0) {
        fprintf(stderr, "FAIL: Expected DIC decrease (degassing) but got increase or no change (ΔDIC=%.6f)\n", 
                DIC_change_1);
        return 1;
    }
    printf("  PASS: DIC decreased (ΔDIC=%.6f mmol/m³), indicating CO2 degassing\n", DIC_change_1);

    // Test scenario 2: Undersaturated water (pCO2_water < pCO2_atm)
    // Should result in invasion (negative CO2 gain to water)
    printf("Test 2: Undersaturated water (invasion expected)\n");
    
    v[DIC].c[1] = 1000.0; // Low DIC
    v[AT].c[1] = 1500.0;  // Higher AT/DIC ratio = higher pH = lower CO2
    
    // Reset other species to avoid interference
    v[Phy1].c[1] = 1.0;
    v[Phy2].c[1] = 1.0;
    v[TOC].c[1] = 100.0;
    
    double DIC_initial_2 = v[DIC].c[1];
    
    Biogeo(DELTI);

    double DIC_change_2 = v[DIC].c[1] - DIC_initial_2;
    
    // In undersaturated conditions, CO2 should invade (DIC should increase from gas exchange alone)
    // Note: NPP will decrease DIC, so the net DIC change might be negative
    // But the CO2 exchange component (co2air[1]) should be positive
    extern double co2air[MAXM + 1];
    
    if (co2air[1] <= 0.0) {
        fprintf(stderr, "FAIL: Expected positive CO2 invasion flux but got co2air=%.6f mmol/m³/s\n", 
                co2air[1]);
        return 1;
    }
    printf("  PASS: CO2 invasion flux positive (co2air=%.6f mmol/m³/s)\n", co2air[1]);

    printf("PASS: CO2 air-water exchange sign convention correct\n");
    return 0;
}
