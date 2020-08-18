#define DIMS 3
#define problemWidth 1.0
#define refinementFactor 8
#define subdividedRectangle true

#define problemHeight 100.0
#define NumMeshPoints 3

#if subdividedRectangle
#undef problemWidth
#undef refinementFactor
#define problemWidth 10
#define refinementFactor 2
#endif

#define maxRefinementLevel 7
#define minRefinementLevel 4

#define isFiniteStrain true
#define isMechanics false
#define isTraction false
#define TimeStep 1.0e-2
#define TotalTime 1500*TimeStep
//grain-structure parameters                                                                  
#define N_seed_points 50
#define n_diff_grains 4
//Allen Cahn parametrs                                                                        
#define InterfaceEnergyParameter 5.0//.0e-4
#define Mobility 30.0
#define Mobility_c 30.0
#define Mobility_m 30.0
#define M_alpha 10.0
#define Vm 1.0
#define TotalDOF n_diff_grains

#if isMechanics
#undef TotalDOF
#define TotalDOF DIMS+n_diff_grains
#endif

#define outputFileName "solution"
#define alpha1 2000
#define beta1 1000
#define PI 3.14159265359
