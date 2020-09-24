//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2012
//authors: rudraa (2012, 2018)
//

#ifndef CHEMO_H_
#define CHEMO_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"

//Chemistry residual implementation
template <int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1,Sacado::Fad::DFad<double> >& R, /*double currentTime, double totalTime,*/ dealii::Table<1,double>& phi_conv,FullMatrix<double>&local_matrix,unsigned int currentIncrement, Sacado::Fad::DFad<double>& free_energy){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  //evaluate gradients 
  dealii::Table<1,Sacado::Fad::DFad<double> > phi(n_q_points), x(n_q_points), mu(n_q_points),mu_conv(n_q_points);//x_conv(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > phi_j(n_q_points, dim), x_j(n_q_points, dim), mu_j(n_q_points,dim);
  dealii::Table<2,double> phi_conv_j(n_q_points, dim),x_conv_j(n_q_points,dim), mu_conv_j(n_q_points,dim);
  Table<1, double> x_conv(n_q_points);
  for (unsigned int q=0; q<n_q_points; ++q) {
    phi[q]=0.0; phi_conv[q]=0.;
    x[q]=0.; x_conv[q]=0.; mu[q]=0.0;

    for (unsigned int j=0; j<dim; j++) {phi_j[q][j]=0.0; phi_conv_j[q][j]=0.0;  x_j[q][j]=0.0;x_conv_j[q][j]=0.;mu_j[q][j]=0.;mu_conv_j[q][j]=0; }
    
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if(ck==0) { phi[q]+=fe_values.shape_value(i,q)*ULocal[i]; phi_conv[q]+=fe_values.shape_value(i,q)*ULocalConv[i];}
      if(ck==1) { x[q]+=fe_values.shape_value(i,q)*ULocal[i]; x_conv[q]+=fe_values.shape_value(i,q)*ULocalConv[i];}
      if(ck==2) {mu[q]+=fe_values.shape_value(i,q)*ULocal[i];mu_conv[q]+=fe_values.shape_value(i,q)*ULocalConv[i];}

      for (unsigned int j=0; j<dim; j++) {
	unsigned int cj=fe_values.get_fe().system_to_component_index(i).first;
	if(cj==0){ phi_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];}	
	if(cj==1){x_j[q][j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];}
	if(cj==2){mu_j[q][j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];}
	
      }
    }
  }
 
Sacado::Fad::DFad<double> M_phi=Mobility;
//if(currentIncrement<35)M_phi=0.4;
Sacado::Fad::DFad<double> UGC=gas_constant; Sacado::Fad::DFad<double> Temp=temperature; Sacado::Fad::DFad<double> epsilon=InterfaceEnergyParameter;Sacado::Fad::DFad<double> M=0.;
  //evaluate Residual for allen cahn
  for(unsigned int q=0;q<n_q_points;q++){
    
    Sacado::Fad::DFad<double> Wx=-2.0*x[q]*x[q] - 5.0*x[q]+7.0;//WB*x[q]+WA*(1.0-x[q]);//WB*x[q]*exp(1.0-x[q]) + WA*(1.0-x[q])*exp(-x[q])+18706.5*(1.0-2.0*x[q]);
    //Sacado::Fad::DFad<double> Wx=WB*x[q] + WA*(1.0-x[q]);
    Sacado::Fad::DFad<double> g_phi=phi[q]*phi[q]*(1.0-phi[q])*(1.0-phi[q]);
    Sacado::Fad::DFad<double> g_d_phi=2.0*(2*pow(phi[q],3)-3.0*pow(phi[q],2)+phi[q]);
    //if(currentIncrement>=210){dt=dt/2.0;}
    if(currentIncrement<30){dt=0.01;Wx=1.0;M_phi=50;}
    for(unsigned int A=0;A<dofs_per_cell;A++){
      unsigned int ca=fe_values.get_fe().system_to_component_index(A).first;
      if(ca==0){
	R[A]+=(1/dt)*fe_values.shape_value(A,q)*(phi[q]-phi_conv[q])*fe_values.JxW(q);
	R[A]+=fe_values.shape_value(A,q)*(M_phi/Vm)*(g_d_phi*Wx)*fe_values.JxW(q);
	for(unsigned int i=0;i<dim;i++){
	  R[A]+=M_phi*Wx*epsilon*fe_values.shape_grad(A,q)[i]*phi_j[q][i]*fe_values.JxW(q);
	}
      }

      if(ca==1){
	// residual R2
	//Gronaghen treatment
	/*M=0.;
	  }*/

	// Fadi's Treatment
       
	M=M_alpha;
	if(currentIncrement<30)M=0.0;
	R[A]+=fe_values.shape_value(A,q)*(x[q]-x_conv[q])*(1.0/dt)*fe_values.JxW(q);
	for(unsigned int i=0;i<dim;i++){
	  R[A]+=M*fe_values.shape_grad(A,q)[i]*mu_j[q][i]*fe_values.JxW(q);
	}
      }
      if(ca==2){
	Sacado::Fad::DFad<double> dG_dX=0.;
	dG_dX=1000*x[q]*(x[q]-1.0)*(x[q]-0.5)+g_phi*(-4.0*x[q]-5.0);//(WB-WA);//( WB*(1.0-x[q])*exp(1.0-x[q]) +WA*(x[q]-2.0)*exp(-x[q])  );
	R[A]+=fe_values.shape_value(A,q)*(mu[q]- dG_dX)*fe_values.JxW(q);
	for(unsigned int i=0;i<dim;i++){
	  R[A]-=kappa*fe_values.shape_grad(A,q)[i]*x_j[q][i]*fe_values.JxW(q);
	}
      }
    }
  }
  

}

#endif /* CHEMO_H_ */
