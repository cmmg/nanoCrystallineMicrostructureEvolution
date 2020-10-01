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
void evaluateFieldAtQuadraturePoint (unsigned int q, Table<1, Sacado::Fad::DFad<double> >&phi, Table<1, Sacado::Fad::DFad<double> >&phi_conv, Table<2, Sacado::Fad::DFad<double> >&phi_j, Sacado::Fad::DFad<double>& sol, Sacado::Fad::DFad<double>& sol_conv, Table<1, Sacado::Fad::DFad<double> >&sol_j, Sacado::Fad::DFad<double>& mu ,Table<1, Sacado::Fad::DFad<double> >& mu_j, FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double >& ULocalConv, dealii::Table<1,Sacado::Fad::DFad<double> >& R,unsigned int currentIncrement){

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
    if(ci==n_diff_grains+1){
      mu+=fe_values.shape_value(i,q)*ULocal[i];
      for(unsigned int j=0;j<dim;j++){
	mu_j[j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];
      }
    }

  }

  

}





//Chemistry residual implementation
template <int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1,Sacado::Fad::DFad<double> >& R, /*double currentTime, double totalTime,*/ FullMatrix<double>&local_matrix,unsigned int currentIncrement, Sacado::Fad::DFad<double>& free_energy){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;

  for (unsigned int q=0; q<n_q_points; ++q) {
    
    Table<1, Sacado::Fad::DFad<double> > phi(n_diff_grains), phi_conv(n_diff_grains), sol_j(dim), mu_j(dim);
    Table<2, Sacado::Fad::DFad<double> > phi_j(n_diff_grains, dim), phi_conv_j(n_diff_grains, dim);
    Sacado::Fad::DFad<double> sol, mu, sol_conv;
    sol=0.0; mu=0.0; sol_conv=0.0;
    
    
    for(unsigned int i=0;i<n_diff_grains;i++){phi[i]=0.0; phi_conv[i]=0.;}
    for(unsigned int i=0;i<dim;i++){sol_j[i]=0.0; mu_j[i]=0.0;}
    for(unsigned int i=0;i<n_diff_grains;i++){
      for(unsigned int j=0;j<dim;j++){
	phi_j[i][j]=0.0;
	phi_conv_j[i][j]=0.0;
      }
    }
    
    //evaluateFieldAtQuadraturePoint (unsigned int q, Table<1, Sacado::Fad::DFad<double> >&phi, Table<1, Sacado::Fad::DFad<double> >&phi_conv, Table<2, Sacado::Fad::DFad<double> >&phi_j, Sacado::Fad::DFad<double>& sol, Sacado::Fad::DFad<double>& sol_conv, Table<1, Sacado::Fad::DFad<double> >&sol_j, Sacado::Fad::DFad<double>& mu ,Table<1, Sacado::Fad::DFad<double> >& mu_j, FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocalConv, dealii::Table<1,Sacado::Fad::DFad<double> >& R,unsigned int currentIncrement)
    
    evaluateFieldAtQuadraturePoint<dim>(q, phi, phi_conv, phi_j, sol, sol_conv, sol_j, mu , mu_j, fe_values, DOF,fe_face_values,cell, dt,ULocal,  ULocalConv,  R,  currentIncrement); 
    //residual calculation
    for(unsigned int A=0; A<dofs_per_cell;A++){
      
      unsigned int ca=fe_values.get_fe().system_to_component_index(A).first;
      // for order parameters i.e. grains fields
      Sacado::Fad::DFad<double> M_phi=Mobility;
      Sacado::Fad::DFad<double> M_sol=M_alpha;
      Sacado::Fad::DFad<double> epsilon=InterfaceEnergyParameter;
      if(ca<n_diff_grains){
	Sacado::Fad::DFad<double> phi2Sum=0.;
	Sacado::Fad::DFad<double> W_sol=-2.0*sol*sol-5.0*sol+7.0;
	if(currentIncrement<30){
	  W_sol=1.0;
	  //M_sol=0.0;
	  dt=0.01;
	}
	for(unsigned int i=0;i<n_diff_grains;i++){
	  phi2Sum+=phi[i]*phi[i];
	}
	// double well for phi is g(phi_1, phi_2..)= (4/3)(1- 4*\Sigma \phi^3 + 3 * (\Sigma \phi^2)^2)
	R[A]+=(1.0/dt)*fe_values.shape_value(A,q)*(phi[ca]-phi_conv[ca])*fe_values.JxW(q);
	R[A]+=M_phi*W_sol* fe_values.shape_value(A,q)*(4./3.)*(-12.0*phi[ca]*phi[ca]+12.0* phi2Sum*phi[ca])*fe_values.JxW(q);
	for(unsigned int j=0;j<dim;j++){
	  R[A]+=M_phi*W_sol*epsilon*fe_values.shape_grad(A,q)[j]*phi_j[ca][j]*fe_values.JxW(q);
	}

      } 
      
      // for solute field                                                                                      
      if(ca==n_diff_grains){
	if(currentIncrement<30)M_sol=0.0;
	R[A]+=(1.0/dt)* (sol-sol_conv)*fe_values.shape_value(A,q)*fe_values.JxW(q);
	for(unsigned int j=0;j<dim;j++){
	  R[A]+=M_sol*fe_values.shape_grad(A,q)[j]*mu_j[j]*fe_values.JxW(q);
	}
      
      } 
      
      //for chemical potential                                                                                 
      if(ca==n_diff_grains+1){
	Sacado::Fad::DFad<double> g_phi=0., phi2Sum=0., phi3Sum=0.;
	for(unsigned int i=0;i<n_diff_grains;i++){
	  phi2Sum+=phi[i]*phi[i];
	  phi3Sum+=phi[i]*phi[i]*phi[i];
	}
	g_phi=(4./3.)*(1.0- 4.0* phi3Sum + 3.0 * phi2Sum*phi2Sum);
	//solute energy double well G
	Sacado::Fad::DFad<double> dG_dsol= 800*sol*(sol-1.0)*(sol-0.5)+ g_phi*(-4.0*sol-5);
	R[A]+=fe_values.shape_value(A,q)* (mu-dG_dsol)*fe_values.JxW(q);
	for(unsigned int j=0;j<dim;j++){
	  R[A]-=kappa*fe_values.shape_grad(A,q)[j]*sol_j[j]*fe_values.JxW(q);
	}
      } 
      
    }

  }
}

#endif /* CHEMO_H_ */
