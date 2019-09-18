
//problem geometry, mesh control
#define DIMS 2
#define problemWidth 1.0
#define refinementFactor 7

//phase field properties
#define InterfaceEnergyParameter {1.0e-5, 1.0e-5, 1.0e-5} //{Kx, Ky, Kz}
#define TimeStep 1.0e-1
#define TotalTime 10000*TimeStep


#define L 50.0
#define alpha1 1.0

#define n_seed_points 20
#define n_diff_grains 5

#define TotalDOF n_diff_grains
//output controls
#define outputFileName "solution"

