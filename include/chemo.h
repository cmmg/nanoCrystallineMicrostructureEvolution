
//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2012
//authors: rudraa (2012, 2018)
//

#ifndef CHEMO_H_
#define CHEMO_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"

template<int dim>
struct historyVariables{
  historyVariables<dim>(){}
  Sacado::Fad::DFad<double> ID;
};

template<int dim>
void evaluateFieldAtQuadraturePoint(unsigned int q, Table<1, Sacado::Fad::DFad<double> >&phi, Table<1, Sacado::Fad::DFad<double> >&phi_conv, Table<2, Sacado::Fad::DFad<double> >&phi_j, Sacado::Fad::DFad<double>&sol, Sacado::Fad::DFad<double>&sol_conv, Table<1, Sacado::Fad::DFad<double> >&sol_j,Sacado::Fad::DFad<double>&mu ,Table<1, Sacado::Fad::DFad<double> >&mu_j,FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1,Sacado::Fad::DFad<double> >& R,FullMatrix<double>&local_matrix,unsigned int currentIncrement,unsigned int currentIteration ,Sacado::Fad::DFad<double>& free_energy,std::vector<historyVariables<dim>* >& history) {

  unsigned int dofs_per_cell=fe_values.dofs_per_cell;
  for(unsigned int i=0;i<dofs_per_cell;i++){
    
    unsigned int ci=fe_values.get_fe().system_to_component_index(i).first;
    if(ci<n_diff_grains){
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
    
  }
  
  //dof loop ends
  
  
  if(currentIteration==0){
    for(unsigned int i=0;i<n_diff_grains;i++){
      if(phi_conv[i]>0.93){history[q]->ID=i+1;break;}
      else{history[q]->ID=0;}
    }
  }
  
  //function definition ends
}

template<int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1,Sacado::Fad::DFad<double> >& R, /*double currentTime, double totalTime,*/ dealii::Table<2,Sacado::Fad::DFad<double> >& phi_conv,FullMatrix<double>&local_matrix,unsigned int currentIncrement, Sacado::Fad::DFad<double>& free_energy, unsigned int currentIteration ,std::vector<historyVariables<dim>* >& history,  Sacado::Fad::DFad<double>& faceEnergy ){

  unsigned int dofs_per_cell=fe_values.dofs_per_cell;
  unsigned int n_q_points=fe_values.n_quadrature_points;

  /*if(currentIncrement<20){
    M_alpha=0.;
  }*/
  
  for(unsigned int q=0;q<n_q_points;q++){
    Table<1, Sacado::Fad::DFad<double> >phi(n_diff_grains), phi_conv(n_diff_grains), sol_j(dim), mu_j(dim);
    Table<2,Sacado::Fad::DFad<double> >phi_j(n_diff_grains, dim);
    for(unsigned int i=0;i<n_diff_grains;i++){phi[i]=0.; phi_conv[i]=0.;}
    for(unsigned int i=0;i<dim;i++){sol_j[i]=0.; mu_j[i]=0.;}

    for(unsigned int i=0;i<n_diff_grains;i++){
      for(unsigned int j=0;j<dim;j++){
	phi_j[i][j]+=0.;
      }
    }
    
    Sacado::Fad::DFad<double> mu=0.0, sol=0.0, sol_conv=0.0;
    evaluateFieldAtQuadraturePoint<dim>( q, phi, phi_conv, phi_j, sol, sol_conv, sol_j, mu , mu_j, fe_values,  DOF,fe_face_values,cell, dt,ULocal,  ULocalConv,  R,local_matrix,  currentIncrement,currentIteration, free_energy, history) ;
    Sacado::Fad::DFad<double> epsilon=InterfaceEnergyParameter;
    Sacado::Fad::DFad<double> M=M_alpha, M_phi=Mobility;
    if(currentIncrement<20){
      M_phi=Mobility;
      M=0.;
      //dt=0.01;
     }
    if(currentIncrement>=20 && currentIncrement<50 ){
      M_phi=0.;
      M=50.0;
      //dt=0.01;
    }
    for(unsigned int i=0;i<dofs_per_cell;i++){
      unsigned int ci=fe_values.get_fe().system_to_component_index(i).first;
      
      if(ci<n_diff_grains){
	Sacado::Fad::DFad<double> phi2_sum=0., phi4_sum=0.,sq_sum=0.;
	Sacado::Fad::DFad<double> W_sol=0.;
	W_sol=WA*(1.0-sol)+WB*sol;
	if(currentIncrement<50){W_sol=1.0;}
	for(unsigned int I=0;I<n_diff_grains;I++){
	  phi2_sum+=pow(phi[I],2);
	}
        
	R[i]+=(1/dt)*fe_values.shape_value(i,q)*(phi[ci]-phi_conv[ci])*fe_values.JxW(q);
	R[i]+=fe_values.shape_value(i,q)*(4./3.)*(M_phi/Vm)*(-12.0*phi[ci]*phi[ci]+12.0*phi[ci]*phi2_sum)*W_sol*fe_values.JxW(q);

	for(unsigned int j=0;j<dim;j++){
	  R[i]+=(M_phi*W_sol)*epsilon*fe_values.shape_grad(i,q)[j]*phi_j[ci][j]*fe_values.JxW(q);
	}
	
      }
      if(ci==n_diff_grains){
      
	
	//pow(M_gb/M_alpha,4*c[q]*(1-c[q]));
	Sacado::Fad::DFad<double>g_phi=0.,phi2_sum=0.,phi3_sum=0.,phi4_sum=0., sq2_sum=0., sq_sum=0.;
	Sacado::Fad::DFad<double> dG_dX_2=0.,W_sol=0.;
	
	Sacado::Fad::DFad<double> g_d_phi=0.;
	Sacado::Fad::DFad<double> dG_dX_dPhi=0.;
	for(unsigned int I=0;I<n_diff_grains;I++){
	  phi2_sum+=pow(phi[I],2);
	  phi3_sum+=pow(phi[I],3);
	  phi4_sum+=pow(phi[I],4);
	}
	g_phi=(4./3.)*(1.0-4*phi3_sum+3.0*phi2_sum*phi2_sum);
	W_sol=WA*(1.0-sol)+WB*sol;
	dG_dX_2=1/sol- 1.0/(sol-1.0);
	Sacado::Fad::DFad<double> temp=0.;Table<1, Sacado::Fad::DFad<double> >grad_temp(dim);

	for(unsigned int j=0;j<dim;j++)grad_temp[j]=0.;

	for(unsigned int I=0;I<n_diff_grains;I++){
	  temp=0.;
	  temp=(4./3.)*(-12.0* phi[I]*phi[I] + 12.0* phi[I]*phi2_sum)*(WB-WA);
	  for(unsigned int j=0;j<dim;j++){
	    grad_temp[j]+=temp*phi_j[I][j];
	  }
	}
	
	Table<1, Sacado::Fad::DFad<double> >J(dim);
	for(int j=0;j<dim;j++)J[j]=0.;

	for(unsigned int j=0;j<dim;j++){
	  J[j]+=(-1.0/Vm)*M*(sol)*(1.0-sol)*(dG_dX_2*sol_j[j]+grad_temp[j])*fe_values.JxW(q);
	}
	
	R[i]+=(1/dt)*fe_values.shape_value(i,q)*(sol-sol_conv)*fe_values.JxW(q);
	for(unsigned int j=0;j<dim;j++){
	  R[i]-=Vm*fe_values.shape_grad(i,q)[j]*J[j]*fe_values.JxW(q);
	}	
      }	
  
    
      if(ci>=0 && ci <n_diff_grains){
	for(unsigned int j=0;j<dim;j++){
	  faceEnergy+= epsilon*phi_j[ci][j]*phi_j[ci][j]*fe_values.JxW(q);
	}
	
      }
      
      
    }
  }



}
  





#endif /* CHEMO_H_ */
