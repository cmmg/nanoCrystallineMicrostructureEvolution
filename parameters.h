//problem geometry, mesh control
#define DIMS 2
#define problemWidth 1.0
#define refinementFactor 7

//phase field properties
//#define InterfaceEnergyParameter 1.0e-2//8.0e-1// {1.0e-2, 1.0e-2, 1.0e-2} //{Kx, Ky, Kz}
//#define dFdC  400*c[q]*(c[q]-1.0)*(c[q]-0.5) //derivative of the free energy
//#define Mobility 2.0//3.72e-4//0.5//*3.0e-12
#define M_gb     1e-4//2.64e-11//0.04
//#define M_alpha  2.0e-2//0.0e-2 //0.001

//time step controls
#define TimeStep 1.0e-2
#define TotalTime 10000*TimeStep

#define n_phase 1
#define n_solute 1
#define TotalDOF n_phase+n_solute
#define Vm 1.0//7.0e-6 //0.001
//#define WA 100//1088//1000
//#define WB -230.7615//-1000//900
#define gas_constant 8.314
#define temperature 900
#define GA 0
#define GB 0
//Vm is molar volume
//output controls
#define outputFileName "solution"

#define Mobility 25.0
#define M_alpha 0.1//1e-1
#define InterfaceEnergyParameter 1.0e-3
#define WA 0.1
#define WB -0.9
