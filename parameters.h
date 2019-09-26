//problem geometry, mesh control
#define DIMS 2
#define problemWidth 1
#define refinementFactor 7

//phase field properties
#define InterfaceEnergyParameter 5.0e-2 //2*8.0e-1// {1.0e-2, 1.0e-2, 1.0e-2} //{Kx, Ky, Kz}
//#define dFdC  400*c[q]*(c[q]-1.0)*(c[q]-0.5) //derivative of the free energy
#define kappa 1.0e-2
#define Mobility 20.0 //3.72e-4//0.5//*3.0e-12
#define M_gb     1e-4//2.64e-11//0.04
#define M_alpha  1 //1.0e-2//1.0e-4//1.0e-4 //0.001

//time step controls
#define TimeStep 5.0e-7
#define TotalTime 10000*TimeStep


#define n_solute 1
#define Vm 1.0//7.0e-6 //0.001
#define WA 100//1088//1000
#define WB -230.7615//-1000//900
#define gas_constant 8.314
#define temperature 900
#define GA 0
#define GB 0
#define N_seed_points 500
#define n_diff_grains 9
#define ch_potential 1
#define W_phi 50
#define zeta1 2
#define zeta2 -2
#define TotalDOF n_diff_grains+n_solute+ch_potential
//Vm is molar volumex
//output controls
#define outputFileName "solution"
