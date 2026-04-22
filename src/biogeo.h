/**
 * @file biogeo.h
 * @brief Function declarations for biogeochemical processes
 * @author An Nguyen
 * @date Last updated: 05/2025
 */

#ifndef BIOGEO_H
#define BIOGEO_H

#include "define.h"
#include "variables.h"

// Main biogeochemical calculation function
// Note: when enable_biogeochemical_reactions=0 and enable_carbonate_diagnostics=1,
// this function will still compute carbonate speciation diagnostics (pH/CO2/pCO2)
// from transported DIC+AT without changing DIC/AT mass.
void Biogeo(int t);

// Core biogeochemical process functions
void computePrimaryProduction(int t, int i);
void computeBiogeochemicalReactions(int t, int i);
void computeCarbonateChemistry(int t, int i);
void updateBiogeochemicalState(int t, int i);
void updateSuspendedSediment(int t);


// Helper and utility functions
double waterT(int t);
double I0(int t);
double Pbmax(int t, int s);
double kmort(int t, int s);
double kmaintenance(int t, int s);
double Fhetox(int t);
double Fhetden(int t);
double Fnit(int t);
double O2sat(int t, int i);
double Sc(int t, int i);

// Carbonate system calculations
double disK1(double Sal, double t);
double disK2(double Sal, double t);
double disK3(double Sal, double t);
double disK4(double Sal, double t);
double disK5(double Sal, double t);

double KH(int t, int i);
double Diff(int t);

#endif // BIOGEO_H
