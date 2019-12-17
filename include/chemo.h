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
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1,Sacado::Fad::DFad<double> >& R, /*double currentTime, double totalTime,*/ dealii::Table<1,double>& c_conv,FullMatrix<double>&local_matrix,unsigned int currentIncrement, Sacado::Fad::DFad<double>& free_energy, Sacado::Fad::DFad<double>& bulkEnergy, Sacado::Fad::DFad<double>& volume, Sacado::Fad::DFad<double>& faceEnergy){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
    
  //evaluate gradients 
  dealii::Table<1,Sacado::Fad::DFad<double> > c(n_q_points), x(n_q_points);//x_conv(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > c_j(n_q_points, dim), x_j(n_q_points, dim);
  dealii::Table<2,double> c_conv_j(n_q_points, dim),x_conv_j(n_q_points,dim);
  Table<1, double> x_conv(n_q_points);
  for (unsigned int q=0; q<n_q_points; ++q) {
    c[q]=0.0; c_conv[q]=0.;
    x[q]=0.; x_conv[q]=0.;// mu[q]=0.0; not needed 
    for (unsigned int j=0; j<dim; j++) {c_j[q][j]=0.0; c_conv_j[q][j]=0.0;  x_j[q][j]=0.0;x_conv_j[q][j]=0.; }
    
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if(ck==0) { c[q]+=fe_values.shape_value(i,q)*ULocal[i]; c_conv[q]+=fe_values.shape_value(i,q)*ULocalConv[i];}
      if(ck==1) { x[q]+=fe_values.shape_value(i,q)*ULocal[i]; x_conv[q]+=fe_values.shape_value(i,q)*ULocalConv[i];}
      
      for (unsigned int j=0; j<dim; j++) {
	unsigned int cj=fe_values.get_fe().system_to_component_index(i).first;
	if (cj==0){ c_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];}
	
	if(cj==1){x_j[q][j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];}
	
      }
    }
  }

  Sacado::Fad::DFad<double> M_phi=Mobility;
  Sacado::Fad::DFad<double> M=M_alpha;
  Sacado::Fad::DFad<double> UGC=gas_constant; Sacado::Fad::DFad<double> Temp=temperature; Sacado::Fad::DFad<double> epsilon=InterfaceEnergyParameter;
  //if(currentIncrement<40){M=0.;}
  //if(currentIncrement>40 && currentIncrement<40){M_phi=0.;M=4*M;}
  //if(currentIncrement>=30){M=0.;}
  for(unsigned int q=0;q<n_q_points;q++){
    
    Sacado::Fad::DFad<double> Wx=WB*(x[q])+WA*(1.0-x[q]);
  
    Sacado::Fad::DFad<double> g_phi=pow(c[q],2)*pow((1.0-c[q]),2);
    Sacado::Fad::DFad<double> g_d_phi=4.0*pow(c[q],3)-6.0*pow(c[q],2)+2.0*c[q]; 
      
    Sacado::Fad::DFad<double> dG_dX2=(1.0/x[q])-1.0/(x[q]-1.0);
    
    Sacado::Fad::DFad<double> dG_dX_dPhi=g_d_phi*(WB-WA);
    if(currentIncrement<30){Wx=1.0; dG_dX_dPhi=g_d_phi;M=0.;}
    
    Sacado::Fad::DFad<double> dG_dX2_dPhi=0.;
    for(unsigned int A=0;A<dofs_per_cell;A++){
      unsigned int ca=fe_values.get_fe().system_to_component_index(A).first;
      if(ca==0){
	R[A]+=(1/dt)*fe_values.shape_value(A,q)*(c[q]-c_conv[q])*fe_values.JxW(q);
	R[A]+=fe_values.shape_value(A,q)*(M_phi/Vm)*(g_d_phi*Wx)*fe_values.JxW(q);
	for(unsigned int i=0;i<dim;i++){
	  R[A]+=(M_phi*Wx)*epsilon*fe_values.shape_grad(A,q)[i]*c_j[q][i]*fe_values.JxW(q);
	}
      }

      if(ca==1){

	Table<1, Sacado::Fad::DFad<double>> J_B(dim);for(unsigned int i=0;i<dim;i++)J_B[i]=0.;
	for(unsigned int i=0;i<dim;i++){
	  
	  J_B[i]+=-(1/Vm)*M*(1.0-x[q])*x[q]*(dG_dX2*x_j[q][i]+dG_dX_dPhi*0.0/*c_j[q][i]*/);
	}
	R[A]+=(1/dt)*fe_values.shape_value(A,q)*(x[q]-x_conv[q])*fe_values.JxW(q);
	for(unsigned int i=0;i<dim;i++){
	  R[A]-=fe_values.shape_grad(A,q)[i]*J_B[i]*fe_values.JxW(q);
	}

      }
    }
  }

   for (unsigned int q=0;q<n_q_points;q++){
    for(unsigned int i=0;i<dofs_per_cell;i++){
      bulkEnergy+=(c[q]*c[q]*(1.0-c[q])*(1.0-c[q]))*fe_values.JxW(q);
    }
  }
  
  for (unsigned int q=0;q<n_q_points;q++){
    for(unsigned int i=0;i<dofs_per_cell;i++){
      volume+=c[q]*fe_values.JxW(q);
    }
  }
  
  for(unsigned int q=0;q<n_q_points;q++){
    for(unsigned int i=0;i<n_q_points;i++){
      for(unsigned int j=0;j<dim;j++){
	faceEnergy+=epsilon*c_j[q][j]*c_j[q][j]*fe_values.JxW(q);
      }
    }
  }
 
  
  
}

#endif /* CHEMO_H_ */
