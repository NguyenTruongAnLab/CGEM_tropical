#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "define.h"
#include "variables.h"
#include <stdbool.h>

// Function declarations
void Transport(int t);
void Dispcoef(int t);
void Disp(double* co, int s);
void Openbound(double* co, int s, int t);
void TVD(double* co, int s);
void applyTributarySourceTerms(double* co, int s, int t);
void Boundflux(int s);
void export_limiter_events_summary(const char *path);

// Helper functions
bool is_milestone_day_for_logging(int t);
void InitializeDispersionAndSalinity(int t, bool resetSalinity);
void CalculateDispersionProfile(double D1, double K, double Qf, double A1, int t);
void WriteDerivedVariables(int t);

#endif // TRANSPORT_H
