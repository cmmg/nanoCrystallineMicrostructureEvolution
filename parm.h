//problem geometry, mesh control                                                
#define DIMS 2
#define problemWidth 1
#define refinementFactor 7

   //phase field properties                                                        
#define InterfaceEnergyParameter 1.0e-3//2*8.0e-1//   {Kx, Ky, Kz}                                                                    
#define kappa  4.0e-3// 5.0e-2                                                  
   //#define dFdC  400*c[q]*(c[q]-1.0)*(c[q]-0.5) //derivative of the free energy  
#define Mobility 1.0//10//3.72e-4//0.5//*3.0e-12                                      
#define M_gb     1e-4//2.64e-11//0.04                                           
#define M_alpha  1.0//1.0e-4//1.0e-2//1.0e-4//1.0e-4 //0.001  
   //time step controls                                                            
#define TimeStep 5.0e-7
#define TotalTime 50000*TimeStep


#define N_seed_points 10
#define n_diff_grains 3
#define n_solute 1
#define ch_potential 1
#define TotalDOF n_diff_grains+n_solute+ch_potential


#define Vm 1.0//7.0e-6 //0.001                                                                
#define gas_constant 8.314
#define temperature 900
#define GA 0
#define GB 0   //Vm is molar volume                                                            
   //output controls                                                               
#define outputFileName "solution"
