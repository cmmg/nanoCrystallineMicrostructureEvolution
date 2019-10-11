
//problem geometry, mesh control
#define DIMS 3
#define problemWidth 1.0
#define refinementFactor 8

//phase field properties
#define InterfaceEnergyParameter {5.0e-6, 5.0e-6, 5.0e-6} //{Kx, Ky, Kz}
#define TimeStep 1.0e-2
#define TotalTime 10000*TimeStep


#define L 300.0
#define alpha1 1.0

//#define n_seed_points 150
#define n_diff_grains 5
#define N_seed_points 5
#define TotalDOF n_diff_grains
//output controls
#define outputFileName "solution"

