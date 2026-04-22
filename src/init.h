#ifndef INIT_H
#define INIT_H

#include "define.h"

// Core initialization functions
void Init();
void cleanupOldOutputFiles();

// Initialization components
void initializeHydrodynamics();
void initializeTransportVariables();
void initializeBiogeochemistry();
void assignBiogeochemicalRateConstants();

#endif // INIT_H
