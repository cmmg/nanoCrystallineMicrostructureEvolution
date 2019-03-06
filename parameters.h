//problem geometry, mesh control
#define DIMS 3
#define problemWidth 1.0
#define refinementFactor 2

//mechanics properties
#define elasticModulus 1//2.0e11
#define PoissonsRatio 0.3

//time step controls
#define TimeStep 0.1
#define TotalTime 100*TimeStep

//output controls
#define outputFileName "solution"
#define pi 3.1415926535897932

#define isotropic_hardening 0
#define kinematic_hardening 0
#define Yield_stress 0.1
