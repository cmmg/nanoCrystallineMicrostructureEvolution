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
void residualForChemo(FEValues<dim>& fe_values, unsigned int dof, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, double >& ULocal, dealii::Table<1,double >& ULocalConv, Vector<double> &local_rhs, FullMatrix<double>&local_matrix){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //evaluate gradients 
 
  Table<2,double > c(n_q_points,n_diff_grains);         
  Table<3,double > c_j(n_q_points,n_diff_grains, dim);
  Table<2,double >c_conv(n_q_points,n_diff_grains);
 
  for (unsigned int q=0;q<n_q_points;q++) {
    for(unsigned int p=0;p<n_diff_grains;p++){
      c[q][p]=0.0; c_conv[q][p]=0.0;
      for(unsigned int j=0;j<dim;j++){
	c_j[q][p][j]=0.0;
      }
    }

    for(unsigned int i=0;i<dofs_per_cell;i++){
      int ck=fe_values.get_fe().system_to_component_index(i).first - dof ;
      if(ck==0){
	c[q][ck]+=fe_values.shape_value(i,q)*ULocal[i];
	c_conv[q][ck]+=fe_values.shape_value(i,q)*ULocalConv[i];
	for(unsigned int j=0;j<dim;j++){
	  c_j[q][ck][j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];
	}
      }
    }
    
  }
  double Kappa[]=InterfaceEnergyParameter;

  for(unsigned int i=0; i<dofs_per_cell; i++){
    int ck=fe_values.get_fe().system_to_component_index(i).first - dof;
    if(ck==0){
      for(unsigned int q=0;q<n_q_points;q++){
	double sq_sum=0.0;
	/*for(long int p=0;p<n_diff_grains;p++){
	  if(p==(ck-dim))continue;
	  else sq_sum+=(c_conv[q][p])*(c_conv[q][p]);
	  }*/
	
	local_rhs[i]+=fe_values.shape_value(i,q)*(c[q][ck]-c_conv[q][ck])*fe_values.JxW(q);
	//local_rhs[i]+=dt*fe_values.shape_value(i,q)*(-0.0001)*fe_values.JxW(q);
	local_rhs[i]-=L1*alpha1*fe_values.shape_value(i,q)*c_conv[q][ck]*fe_values.JxW(q);
	local_rhs[i]+=L1*beta1*fe_values.shape_value(i,q)*(std::pow(c_conv[q][ck],3))*fe_values.JxW(q);
	local_rhs[i]+=L1*2*gamma1*fe_values.shape_value(i,q)*c_conv[q][ck]*sq_sum*fe_values.JxW(q);
	
	for(unsigned int j=0;j<dim;j++){
	  double Kjj= Kappa[j];
	  double kc_j= c_j[q][ck][j]*Kjj; // Kjj*C_j
	  local_rhs[i]+=L1*fe_values.shape_grad(i,q)[j]*kc_j*fe_values.JxW(q);
	}
	
      }
    }
  }
  

  //find jacobian matrix
  for(unsigned int q=0;q<n_q_points;q++){
    for(unsigned int A=0;A<dofs_per_cell;A++){
      int CA=fe_values.get_fe().system_to_component_index(A).first-dof;
      if(CA==0){
	for(unsigned int B=0;B<dofs_per_cell;B++){
	  int CB=fe_values.get_fe().system_to_component_index(B).first-dof;
	  if(CB==0){
	    double sq_sum=0.;
	    //for(unsigned int p=0;p<n_diff_grains;p++){
	      //if(p==(CB- (unsigned int)dim))continue;
	    //sq_sum+=std::pow(c_conv[q][0],2);
	    //}
	    //local_matrix(A,B)+= fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*(1.0/dt - L1*(alpha1+3*beta1*c_conv[q][0]*c_conv[q][0]+2*gamma1*sq_sum*c_conv[q][0]))*fe_values.JxW(q);
	    local_matrix(A,B)+=(1./dt)*fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*fe_values.JxW(q);
	  }
	}
      }
    }

    // second term
    for(unsigned int A=0;A<dofs_per_cell;A++){
      unsigned int CA=fe_values.get_fe().system_to_component_index(A).first-dof;
      if(CA==0){
	for(unsigned int I=0;I<dim;I++){
	  for(unsigned int B=0;B<dofs_per_cell;B++){
	    unsigned int CB=fe_values.get_fe().system_to_component_index(B).first-dof;
	    if(CB==0){
	      for(unsigned int J=0;J<dim;J++){
		local_matrix(A,B)+=L1*Kappa[J]*fe_values.shape_grad(A,q)[I]*fe_values.shape_grad(B,q)[J]*fe_values.JxW(q);
	      }
	    }
	  }	  
	}
      }
    }
  }
  
  
  
}


#endif /* CHEMO_H_ */
