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
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //evaluate gradients 
  
  Table<2,Sacado::Fad::DFad<double> > c(n_q_points,n_diff_grains);         
  Table<3,Sacado::Fad::DFad<double> > c_j(n_q_points,n_diff_grains, dim);
  Table<2,Sacado::Fad::DFad<double> >c_conv(n_q_points,n_diff_grains);
  for (unsigned int q=0;q<n_q_points;q++) {
    for(unsigned int p=0;p<n_diff_grains;p++){
      c[q][p]=0.0;c_conv[q][p]=0.0;
      for(unsigned int j=0;j<dim;j++){
	c_j[q][p][j]=0.0;
      }
    }
    for(unsigned int i=0;i<dofs_per_cell;i++){ 
      unsigned int ck=fe_values.get_fe().system_to_component_index(i).first;
      c[q][ck]+=fe_values.shape_value(i,q)*ULocal[i];
      c_conv[q][ck]+=fe_values.shape_value(i,q)*ULocalConv[i];
      for(unsigned int j=0;j<dim;j++){
	c_j[q][ck][j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];
      }
    }
    
  }
  double Kappa[]=InterfaceEnergyParameter;
  
  for(unsigned int i=0; i<dofs_per_cell; i++){
    unsigned int ck=fe_values.get_fe().system_to_component_index(i).first;
    if(ck>dim){
      for(unsigned int q=0;q<n_q_points;q++){
	double sq_sum=0.0;
	for(long int p=0;p<n_diff_grains;p++){
	  if(p==(ck-dim))continue;
	  else sq_sum+=(c_conv[q][p])*(c_conv[q][p]);
	}
	
	R[i]+=(1/dt)*fe_values.shape_value(i,q)*(c[q][ck]-c_conv[q][ck])*fe_values.JxW(q);
	R[i]-=L1*alpha1*fe_values.shape_value(i,q)*c_conv[q][ck]*fe_values.JxW(q);
	R[i]+=L1*beta1*fe_values.shape_value(i,q)*(std::pow(c_conv[q][ck],3))*fe_values.JxW(q);
	R[i]+=L1*2*gamma1*fe_values.shape_value(i,q)*c_conv[q][ck]*sq_sum*fe_values.JxW(q);
	
	for(unsigned int j=0;j<dim;j++){
	  Sacado::Fad::DFad<double> Kjj= Kappa[j];
	  Sacado::Fad::DFad<double> kc_j= c_j[q][ck][j]*Kjj; // Kjj*C_j
	  R[i]+=L1*fe_values.shape_grad(i,q)[j]*kc_j*fe_values.JxW(q);
	}
	
      }
    }
  }
}


#endif /* CHEMO_H_ */
