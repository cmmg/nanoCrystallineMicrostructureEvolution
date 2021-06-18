//mesh parameters
#define DIMS 3
#define refinementFactor 6
#define maxRefinementLevel 7
#define minRefinementLevel 3
#define problemWidth 1.0

// mechanics and solte drag control
#define isMechanics true
#define isSoluteDrag false
#define isFiniteStrain false

// grain growth parameters
#define N_seed_points 40
#define n_diff_grains 4
#define InterfaceEnergyParameter 2.50e-3

// Solute parameters
#define WA 0.9
#define WB 0.1
#define kappa 5.0e-3
#define n_solute 1
#define n_chemical_potential 1

//kinetic parameters
#define TimeStep 5.0e-7
#define TotalTime 1500*TimeStep
#define Mobility 5.0e5
#define M_alpha 1.0

//other parameters and variables
#define PI 3.1415
#define outputFileName "solution"
#define Vm 1.0
#define wellHeight 1.0

// elastic modulii
#define alpha1 20000
#define beta1 10000

// degree of freedom per node
#define TotalDOF n_diff_grains

#if isMechanics
#undef TotalDOF
#define TotalDOF n_diff_grains + DIMS
#endif

#if isSoluteDrag
#undef TotalDOF
#define TotalDOF n_diff_grains + n_solute + n_chemical_potential
#endif

#if (isMechanics && isSoluteDrag)
#undef TotalDOF
#define TotalDOF DIMS + n_diff_grains + n_solute + n_chemical_potential
#endif
