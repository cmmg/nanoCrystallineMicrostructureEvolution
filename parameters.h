//problem geometry, mesh control
#define DIMS 2
#define problemWidth 1.0
#define refinementFactor 7

//mechanics properties
#define elasticModulus 20//2.0e11
#define PoissonsRatio 0.3

//time step controls
#define TimeStep 1.0e-2
#define TotalTime 1000*TimeStep
#define N_seed_points 10
#define n_diff_grains 2
//#define yield_stress 1.0
#define kappa1 5.0e-3
#define Vm 1.0
#define InterfaceEnergyParameter 5.0e-3
#define Mobility 10.0//10.0//10.0// 10.0//50.0
#define M_alpha 0.0//0.1//0.1
#define lambda1 0.1//1.0e-2
#define n_solute 1
#define n_chemical_potential 1
#define TotalDOF DIMS+n_diff_grains//+n_solute+n_chemical_potential
//output controls
#define outputFileName "solution"
#define alpha1 1000
#define beta1 5000
#define PI 3.14159265359
