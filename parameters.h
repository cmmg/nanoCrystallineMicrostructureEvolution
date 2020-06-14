//problem geometry, mesh control
#define DIMS 2
#define problemWidth 1.0
#define refinementFactor 5
#define maxRefinementLevel 7
#define minRefinementLevel 5

//flags mechanics control
#define isFiniteStrain false
#define isMechanics false
#define isTraction false
//time step controls
#define TimeStep 1.0e-2
#define TotalTime 1000*TimeStep
//grain-structure parameters
#define N_seed_points 40
#define n_diff_grains 6
//Allen Cahn parametrs
#define InterfaceEnergyParameter 5.0e-4
#define Mobility 30.0
#define Mobility_c 30.0
#define Mobility_m 30.0

#define TotalDOF n_diff_grains

#if isMechanics
#undef TotalDOF
#define TotalDOF DIMS+n_diff_grains
#endif

#define outputFileName "solution"
#define alpha1 2000
#define beta1 1000
#define PI 3.14159265359
