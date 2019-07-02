//problem geometry, mesh control
#define DIMS 3
#define problemWidth 1.0
#define refinementFactor 0


//mechanics properties
#define elasticModulus 4.2//2.0e11
#define PoissonsRatio 0.3
#define TimeStep 0.0001
//time step controls

#define TotalTime 1000*TimeStep

//output controls
#define outputFileName "solution"


//#define Yield_stress 0.1
#define n_slip_systems 1
#define PI 3.14159265
#define self_hardening  0.0180
#define Ss 0.0148
#define n_seed_points 1
#define n_diff_grains 1
#define totalDOF DIMS+n_diff_grains
#define InterfaceEnergyParameter {1.0e-2, 1.0e-2, 1.0e-2} //{Kx, Ky, Kz}
#define L1 5
#define alpha1 1
#define beta1 1
#define gamma1 1
//#define dFdC  400*c[q]*(c[q]-1.0)*(c[q]-0.5) //derivative of the free energy
#define Mobility 1.0
