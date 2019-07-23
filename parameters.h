//problem geometry, mesh control
#define DIMS 3
#define problemWidth 10
#define problemHeight 100
#define NumMeshPoints 4
#define refinementFactor 1

//mechanics properties
#define elasticModulus 4//2.0e11
#define PoissonsRatio 0.3
#define TimeStep 0.1
//time step controls
#define TotalTime 1000*TimeStep

//solver controls
#define isDirectSolver true

//output controls
#define outputFileName "solution"


//#define Yield_stress 0.1
#define n_slip_systems 1
#define PI 3.14159265
#define self_hardening  0.0180
#define Ss 0.0148
#define n_seed_points 70
#define n_diff_grains 5
#define totalDOF n_diff_grains+DIMS
#define InterfaceEnergyParameter {1.0e-3, 1.0e-3, 1.0e-3} //{Kx, Ky, Kz}
//#define InterfaceEnergyParameter {2, 2, 2} //{Kx, Ky, Kz}
#define L1 0.2
#define alpha1 1
#define beta1 1
#define gamma1 1

