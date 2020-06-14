//new
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2012
//authors: rudraa (2012, 2018)
//

#ifndef CHEMO_H_
#define CHEMO_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"



template<int dim>
void evaluateFieldAtQuadraturePoint(unsigned int q, Table<1, double >&phi, Table<1, double >&phi_conv, Table<2, double >&phi_j, double&sol, double&sol_conv, Table<1, double >&sol_j,double&mu ,Table<1, double >&mu_j,FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,double >& ULocal, dealii::Table<1, double>& ULocalConv, Vector<double>& R,FullMatrix<double>&local_matrix,unsigned int currentIncrement,unsigned int currentIteration ,std::vector<historyVariables<dim>* >& history, double &freeEnergyChemBulk, double & freeEnergyChemGB) {
  int var=0; if(isMechanics)var=dim;
  unsigned int dofs_per_cell=fe_values.dofs_per_cell;
  for(unsigned int i=0;i<dofs_per_cell;i++){
    
    int ci=fe_values.get_fe().system_to_component_index(i).first - var;
    if(ci>=0 && ci<n_diff_grains){
      phi[ci]+=fe_values.shape_value(i,q)*ULocal[i];
      phi_conv[ci]+=fe_values.shape_value(i,q)*ULocalConv[i];
      for(unsigned int j=0;j<dim;j++){
	phi_j[ci][j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];
      }
    }
    if(ci==n_diff_grains){
      sol+=fe_values.shape_value(i,q)*ULocal[i];
      sol_conv+=fe_values.shape_value(i,q)*ULocalConv[i];
      for(unsigned int j=0;j<dim;j++){
	sol_j[j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];
      }
    }
    if(ci==n_diff_grains+1){
      mu+=fe_values.shape_value(i,q)*ULocal[i];
      for(unsigned int j=0;j<dim;j++){
	mu_j[j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];
      }
    }

    //dof loop ends
  }

  double phiSQ=0.;
  for(unsigned int N=0;N<n_diff_grains;N++){
    phiSQ+=phi[N]*phi[N];
    freeEnergyChemBulk+=(4./3.)*(1.0-4.0* phi[N]*phi[N]*phi[N])*fe_values.JxW(q);
  }
  freeEnergyChemBulk+=4.0* phiSQ*phiSQ*fe_values.JxW(q);

  for(unsigned int N=0;N<n_diff_grains;N++){
    for(unsigned int j=0;j<dim;j++){
      freeEnergyChemGB+=(InterfaceEnergyParameter/2.)*(phi_j[N][j]*phi_j[N][j])*fe_values.JxW(q);   
    }
  }


  // definition ends
}

template<int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,double >& ULocal, dealii::Table<1, double>& ULocalConv, Vector<double>& R,FullMatrix<double>&local_matrix,unsigned int currentIncrement,  unsigned int currentIteration ,std::vector<historyVariables<dim>* >& history, double & freeEnergyChemBulk, double & freeEnergyChemGB){
 
  unsigned int dofs_per_cell=fe_values.dofs_per_cell;
  unsigned int n_q_points=fe_values.n_quadrature_points;
 
  for(unsigned int q=0;q<n_q_points;q++){
    Table<1, double >phi(n_diff_grains), phi_conv(n_diff_grains), sol_j(dim), mu_j(dim);
    Table<2,double >phi_j(n_diff_grains, dim);
    
    for(unsigned int i=0;i<n_diff_grains;i++){phi[i]=0.; phi_conv[i]=0.;}
    for(unsigned int i=0;i<dim;i++){sol_j[i]=0.; mu_j[i]=0.;}
    for(unsigned int i=0;i<n_diff_grains;i++){
      for(unsigned int j=0;j<dim;j++){
	phi_j[i][j]=0.;
      }
    }
   
    double mu=0.0, sol=0.0, sol_conv=0.0;
    evaluateFieldAtQuadraturePoint<dim>( q, phi, phi_conv, phi_j, sol, sol_conv, sol_j, mu , mu_j, fe_values,  DOF,fe_face_values,cell, dt,ULocal,  ULocalConv,  R,local_matrix,  currentIncrement,currentIteration, history, freeEnergyChemBulk, freeEnergyChemGB) ;
   
    double epsilon=InterfaceEnergyParameter;
    double M=M_alpha, M_phi=Mobility_c;
    
    
    for(unsigned int i=0;i<dofs_per_cell;i++){
      int ci=fe_values.get_fe().system_to_component_index(i).first - DOF;
      
      if(ci>=0 && ci<n_diff_grains){
	double phi2_sum=0.;
	double W_sol=0.;
	W_sol=1.0 ;
	for(unsigned int I=0;I<n_diff_grains;I++){
	  phi2_sum+=pow(phi[I],2);
	}
        
	R[i]+=(1/dt)*fe_values.shape_value(i,q)*(phi[ci]-phi_conv[ci])*fe_values.JxW(q);
	R[i]+=fe_values.shape_value(i,q)*(4./3.)*(M_phi/Vm)*(-12.0*phi[ci]*phi[ci]+12.0*phi[ci]*phi2_sum)*W_sol*fe_values.JxW(q);
	for(unsigned int j=0;j<dim;j++){
	  R[i]+=(M_phi)*epsilon*fe_values.shape_grad(i,q)[j]*phi_j[ci][j]*fe_values.JxW(q);
	}
	
      }
     
    }

    double g_phi=0., phi2_sum=0., phi3_sum=0.;
    for(unsigned int N=0;N<n_diff_grains;N++){
      phi2_sum+=pow(phi[N],2);
      phi3_sum+=pow(phi[N],3);
    }

    for(unsigned int A=0;A<dofs_per_cell;A++){
      for(unsigned int B=0;B<dofs_per_cell;B++){
	int ca=fe_values.get_fe().system_to_component_index(A).first- DOF;
	int cb=fe_values.get_fe().system_to_component_index(B).first - DOF;
	if(ca>=0 && ca<n_diff_grains ){
	  if(cb>=0 && cb<n_diff_grains){
	    if(ca==cb){
	      
	      local_matrix(A,B)+=fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*(1.0/dt)*fe_values.JxW(q);
	      local_matrix(A,B)+=M_phi*fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*(4./3.)*(-24.0*phi[ca]+12.0*phi2_sum+24.0*phi[ca]*phi[ca])*fe_values.JxW(q);
	      for(unsigned int i=0;i<dim;i++){
		local_matrix(A,B)+=epsilon*M_phi*fe_values.shape_grad(A,q)[i]*fe_values.shape_grad(B,q)[i]*fe_values.JxW(q);
	      }
	      
	    }
	    else{
	      local_matrix(A,B)+=M_phi*fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*(96.0/3.0)*(phi[ca]*phi[cb])*fe_values.JxW(q);
	    }
	    
	  }//if cj ends
	  
	}//if ci ends
	
	
	
      }//dofs J ends
      
      
    }//dofs I ends

  }
    
}

#endif /* CHEMO_H_ */

