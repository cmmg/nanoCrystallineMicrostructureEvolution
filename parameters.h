//problem geometry, mesh control
#define DIMS 3
#define problemWidth 1.0
#define refinementFactor 2

//mechanics properties
#define elasticModulus 1//2.0e11
#define PoissonsRatio 0.3
#define TimeStep 0.1
//time step controls

#define TotalTime 1000*TimeStep

//output controls
#define outputFileName "solution"


#define isotropic_hardening 0
#define kinematic_hardening 0
#define Yield_stress 0.1
#define n_slip_systems 1
