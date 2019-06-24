//problem geometry, mesh control
#define DIMS 3
#define problemWidth 1.0
#define refinementFactor 4

//mechanics properties
#define elasticModulus 4.2//2.0e11
#define PoissonsRatio 0.3
#define TimeStep 0.1
//time step controls

#define TotalTime 1000*TimeStep

//output controls
#define outputFileName "solution"


//#define Yield_stress 0.1
#define n_slip_systems 4
#define PI 3.14159265
#define self_hardening  0.0180
#define Ss 0.0148
#define n_seed_points 20
#define n_diff_grains 10

