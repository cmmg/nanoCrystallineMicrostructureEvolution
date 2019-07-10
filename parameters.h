//problem geometry, mesh control
#define DIMS 3
#define problemWidth 1.0
#define refinementFactor 2


//mechanics properties
#define elasticModulus 4//2.0e11
#define PoissonsRatio 0.3
#define TimeStep 0.1
//time step controls

#define TotalTime 1000*TimeStep

//output controls
#define outputFileName "solution"


//#define Yield_stress 0.1
#define n_slip_systems 1
#define PI 3.14159265
#define self_hardening  0.0180
#define Ss 0.0148
#define n_seed_points 30
#define n_diff_grains 2
#define totalDOF n_diff_grains+DIMS
#define InterfaceEnergyParameter {1.0e-3, 1.0e-3, 1.0e-3} //{Kx, Ky, Kz}
//#define InterfaceEnergyParameter {2, 2, 2} //{Kx, Ky, Kz}
#define L1 1.0
#define alpha1 1
#define beta1 1
#define gamma1 1
#define Mobility 1.0
