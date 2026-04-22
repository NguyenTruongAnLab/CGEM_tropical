/**
 * @file biogeout.c
 * @brief Forcing functions for biogeochemistry (light, temperature responses)
 */

#include "define.h"
#include "variables.h"
#include <math.h>
#include <stdbool.h>

extern double interpolateInputData(int t, double* timeArray, double* dataArray, int dataSize,
                                  bool isHourly, const char* dataName, double defaultValue);

/**
 * @brief Tropical temperature stability cap.
 * Prevents runaway kinetics above 32°C via Gaussian decline (sigma=5°C).
 * @param T Water temperature [°C]
 * @return Scale factor [0, 1]
 */
static double tropical_temperature_limit(double T) {
    if (T > 32.0) {
        return exp(-pow(T - 32.0, 2) / 50.0);
    }
    return 1.0;
}

//Temperature
double waterT(int t)
{
    double default_temp = 20.0;  // Default temperature in °C

    if (!forcingData[FORCING_TEMPERATURE].data || forcingData[FORCING_TEMPERATURE].dataSize <= 0) {
        return default_temp;
    }

    return interpolateInputData(t,
                               forcingData[FORCING_TEMPERATURE].time,
                               forcingData[FORCING_TEMPERATURE].data,
                               forcingData[FORCING_TEMPERATURE].dataSize,
                               false,
                               "Temperature",
                               default_temp);
}

//Light intensity [mu E/ (m2 s)]
double I0(int t)
{
    double default_light = 0.0;
    
    // Check if we have light data available
    if (forcingData[FORCING_LIGHT].data && forcingData[FORCING_LIGHT].time && 
        forcingData[FORCING_LIGHT].dataSize > 0) {

        return interpolateInputData(t,
                                   forcingData[FORCING_LIGHT].time,
                                   forcingData[FORCING_LIGHT].data,
                                   forcingData[FORCING_LIGHT].dataSize,
                                   true,
                                   "Light",
                                   default_light);
    } else {
        // Fallback to analytical formula if no light data is available
        double day, l, r, ncloud;

        // Use the same (t - WARMUP) clock as the central interpolator.
        // Avoid floor() here to keep the forcing continuous across t=WARMUP.
        day = (double)(t - WARMUP) / (24.0 * 60.0 * 60.0);
        ncloud = 0.5 + 0.3 * sin((day - 80) * 2 * M_PI / 365);
        r = 12 + 4 * sin((day - 80) * 2 * M_PI / 365);
        l = (r - 12 + 24 * sin(2 * M_PI * (t - 6 * 3600) / (3600 * 24)));
        return 4600*ncloud*l/24;
    }
}

// Maximum photosynthetic production [1/s] — temperature-dependent (Q10 = 1.067)
double Pbmax(int t, int s)
{
    double T = waterT(t);
    double a = T - 20.0;
    return Pb[s] * pow(1.067, a) * tropical_temperature_limit(T);
}

// Mortality rate [1/s] — temperature-dependent (Q10 = 1.067)
double kmort(int t, int s)
{
    double T = waterT(t);
    double a = T - 20.0;
    return kmortality[s] * pow(1.067, a) * tropical_temperature_limit(T);
}

// Maintenance rate [1/s] — temperature-dependent
double kmaintenance(int t, int s)
{
    double T = waterT(t);
    return kmaint[s] * exp(0.0322 * (T - 20.0)) * tropical_temperature_limit(T);
}

// Aerobic degradation rate [mmolC/m³·s] — temperature-dependent (Q10 = 2.5)
double Fhetox(int t)
{
    double T = waterT(t);
    double a = (T - 20.0) / 10.0;
    return kox * pow(2.5, a) * tropical_temperature_limit(T);
}

// Denitrification rate [mmolC/m³·s] — temperature-dependent (Q10 = 1.07)
double Fhetden(int t)
{
    double T = waterT(t);
    double a = T - 20.0;
    return kdenit * pow(1.07, a) * tropical_temperature_limit(T);
}

// Nitrification rate [mmolN/m³·s] — temperature-dependent (Q10 = 1.08)
double Fnit(int t)
{
    double T = waterT(t);
    double a = T - 20.0;
    return knit * pow(1.08, a) * tropical_temperature_limit(T);
}


//O2 saturation
double O2sat(int t, int i)
{
	double O2, T, f, lnO2sat;

    T=waterT(t)+273.15;
	lnO2sat=-1.3529996*100.+157228.8/T - 66371490./(T*T) + 12436780000./(T*T*T) - 8621061.*100000./(T*T*T*T);
	f=-0.020573 + 12.142/T -2363.1/(T*T);

	O2=exp(lnO2sat+f*v[Sal].c[i]);

	return O2;
}
//Molecular diffusion coefficient for O2
double Diff(int t)
{
    double Diff, T;

    //T=Tabs(t);
    T=waterT(t)+273.15;

    Diff=(6.35*T-1664.)*1.0e-11;

    return Diff;
}

//Schmidt number
double Sc(int t, int i)
{
    double Sc0, Sc;

    Sc0=1800.6-120.1*waterT(t)+3.7818*(waterT(t)*waterT(t))-0.047608*(waterT(t)*waterT(t)*waterT(t));

    Sc=Sc0*(1+(3.14e-3)*v[Sal].c[i]);


    return Sc;
}

//Henry's constant for CO2
double KH(int t, int i)
{
	double KH, T, f, lnK0;

    T=waterT(t)+273.15;
	lnK0= -574.70126 + 21541.52/T - 0.000147759*(T*T) + 89.892*(log(T)) ;
	f= 0.029941 - 0.00027455*T + 0.00000053407*(T*T);

	KH=exp(lnK0+f*v[Sal].c[i]);   //mol/l.atm

	KH=KH*1e6;  //umol/l.atm (== mmol/m3.atm)

	return KH;
}

//==================================================================
// CHEMICAL EQUILIBRIUM CONSTANTS (Literature-Based)
//==================================================================

//_____Dissociation constant of CO2 (mmol m-3) < RTM : Cai & Wang
double disK1(double Sal, double t)
{
    double T = waterT(t) + 273.15;
    double pK1 = -14.8425 + (3404.71/T) + (0.032786*T);
    double f1 = -0.0230848 - (14.3456/T);
    double f2 = 0.000691881 + (0.429955/T);
    return pow(10, -(pK1 + (f1*pow(Sal,0.5)) + (f2*Sal)));
}

//_____Dissociation constant of CO2 (mmol m-3) <RTM : Cai & Wang
double disK2(double Sal, double t)
{
    double T = waterT(t) + 273.15;
    double pK2 = -6.4980 + (2902.39/T) + (0.02379*T);
    double f3 = -0.458898 + (41.24048/T);
    double f4 = 0.0284743 - (2.55895/T);
    return pow(10, -(pK2 + (f3*pow(Sal,0.5)) + (f4*Sal)));
}

//_____Dissociation constant of water (mmol m-3) <RTM : Cai & Wang
double disK3(double Sal, double t)
{
    double T = waterT(t) + 273.15;
    double lnKw = (-13847.26/T) + (148.9652) - (23.6521*log(T));
    double f = (118.67/T) - (5.977) + (1.0495*log(T));
    double g = -0.01615;
    return exp(lnKw + (f*sqrt(Sal)) + (g*Sal));
}

//_____Dissociation constant of sulfide (mmol m-3) :  Dickson 1990b
double disK4(double Sal, double t)
{
    double lnK4, lnK4sansunite, T, I;
    T = waterT(t) + 273.15;
    I = (19.919*Sal) / (1000 - 1.00198*Sal);
    lnK4 = (-4276.1/T) + 141.328 - (23.093*log(T));
    lnK4sansunite = (((-13856/T) + 324.57 - 47.986*log(T))*pow(I,0.5)) + 
                    (((35474/T) - 771.54 + 114.723*log(T))*I) - 
                    ((2698/T)*pow(I,1.5)) + ((1776/T)*pow(I,2)) + lnK4;
    return exp(lnK4sansunite + lnK4);
}

//_____Dissociation constant of boron (mmol m-3) : Dickson 1990a (idem Sandra)
double disK5(double Sal, double t)
{
    double Ts, lnK5, KK5;
    Ts = waterT(t) + 273.15;
    lnK5 = ((-8966.90 - 2890.53*sqrt(Sal) - 77.942*Sal + 1.728*pow(Sal,1.5) - 0.0996*Sal*Sal)/Ts) +
           (148.0248 + 137.1942*sqrt(Sal) + 1.62142*Sal) + 
           ((-24.4344 - 25.085*sqrt(Sal) - 0.2474*Sal)*log(Ts)) + 
           (0.053105*sqrt(Sal)*Ts);
    KK5 = exp(lnK5);
    return KK5;
}
