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
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1,Sacado::Fad::DFad<double> >& R, /*double currentTime, double totalTime,*/ dealii::Table<1,double>& c_conv,FullMatrix<double>&local_matrix,unsigned int currentIncrement, Sacado::Fad::DFad<double>& free_energy){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
    
  //evaluate gradients 
  dealii::Table<1,Sacado::Fad::DFad<double> > c(n_q_points), x(n_q_points), mu(n_q_points),mu_conv(n_q_points);//x_conv(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > c_j(n_q_points, dim), x_j(n_q_points, dim), mu_j(n_q_points,dim);
  dealii::Table<2,double> c_conv_j(n_q_points, dim),x_conv_j(n_q_points,dim), mu_conv_j(n_q_points,dim);
  Table<1, double> x_conv(n_q_points);
  for (unsigned int q=0; q<n_q_points; ++q) {
    c[q]=0.0; c_conv[q]=0.;
    x[q]=0.; x_conv[q]=0.; mu[q]=0.0;

    for (unsigned int j=0; j<dim; j++) {c_j[q][j]=0.0; c_conv_j[q][j]=0.0;  x_j[q][j]=0.0;x_conv_j[q][j]=0.;mu_j[q][j]=0.;mu_conv_j[q][j]=0; }
    
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
      if(ck==0) { c[q]+=fe_values.shape_value(i,q)*ULocal[i]; c_conv[q]+=fe_values.shape_value(i,q)*ULocalConv[i];}
      if(ck==1) { x[q]+=fe_values.shape_value(i,q)*ULocal[i]; x_conv[q]+=fe_values.shape_value(i,q)*ULocalConv[i];}
      if(ck==2) {mu[q]+=fe_values.shape_value(i,q)*ULocal[i];mu_conv[q]+=fe_values.shape_value(i,q)*ULocalConv[i];}

      for (unsigned int j=0; j<dim; j++) {
	unsigned int cj=fe_values.get_fe().system_to_component_index(i).first;
	if(cj==0){ c_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];}	
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
    
    Sacado::Fad::DFad<double> Wx=WB*x[q]+WA*(1.0-x[q]);//WB*x[q]*exp(1.0-x[q]) + WA*(1.0-x[q])*exp(-x[q])+18706.5*(1.0-2.0*x[q]);
    //Sacado::Fad::DFad<double> Wx=WB*x[q] + WA*(1.0-x[q]);
    Sacado::Fad::DFad<double> g_phi=c[q]*c[q]*(1.0-c[q])*(1.0-c[q]);
    Sacado::Fad::DFad<double> g_d_phi=2.0*(2*pow(c[q],3)-3.0*pow(c[q],2)+c[q]);
    
    for(unsigned int A=0;A<dofs_per_cell;A++){
      unsigned int ca=fe_values.get_fe().system_to_component_index(A).first;
      if(ca==0){
	R[A]+=(1/dt)*fe_values.shape_value(A,q)*(c[q]-c_conv[q])*fe_values.JxW(q);
	R[A]+=fe_values.shape_value(A,q)*(M_phi/Vm)*(g_d_phi*Wx)*fe_values.JxW(q);
	for(unsigned int i=0;i<dim;i++){
	  R[A]+=(M_phi)*epsilon*fe_values.shape_grad(A,q)[i]*c_j[q][i]*fe_values.JxW(q);
	}
      }

      if(ca==1){
	// residual R2
	//Gronaghen treatment
	/*M=0.;
	M=M_alpha;//pow(M_gb/M_alpha,4*c[q]*(1-c[q]));
	Table<1, Sacado::Fad::DFad<double>> J_B(dim);for(unsigned int i=0;i<dim;i++)J_B[i]=0.;
	for(unsigned int i=0;i<dim;i++){
	  //Sacado::Fad::DFad<double> constt=M*(1.0-x[q])*x[q];
	  J_B[i]+=-(1/Vm)*M*(1.0-x[q])*x[q]*(dG_dX2*x_j[q][i]+dG_dX_dPhi*c_j[q][i]);
	}
	R[A]+=(1/dt)*fe_values.shape_value(A,q)*(x[q]-x_conv[q])*fe_values.JxW(q);
	for(unsigned int i=0;i<dim;i++){
	  R[A]-=fe_values.shape_grad(A,q)[i]*J_B[i]*fe_values.JxW(q);
	  }*/

	// Fadi's Treatment
       
	M=M_alpha;
	R[A]+=fe_values.shape_value(A,q)*(x[q]-x_conv[q])*(1.0/dt)*fe_values.JxW(q);
	for(unsigned int i=0;i<dim;i++){
	  R[A]+=M*fe_values.shape_grad(A,q)[i]*mu_j[q][i]*fe_values.JxW(q);
	}
      }
      if(ca==2){
	Sacado::Fad::DFad<double> dG_dX=0.;
	dG_dX=800*x[q]*(x[q]-1.0)*(x[q]-0.5)+g_phi*(WB-WA);//( WB*(1.0-x[q])*exp(1.0-x[q]) +WA*(x[q]-2.0)*exp(-x[q])  );
	R[A]+=fe_values.shape_value(A,q)*(mu[q]- dG_dX)*fe_values.JxW(q);
	for(unsigned int i=0;i<dim;i++){
	  R[A]-=kappa*fe_values.shape_grad(A,q)[i]*x_j[q][i]*fe_values.JxW(q);
	}
      }
    }
  }
  

  //define K_local matrix

  /* for(unsigned int q=0;q<n_q_points;q++){
    double Wx=WB*x[q]+WA*(1.0-x[q]);
    double g_phi=c[q]*c[q]*(1-c[q])*(1-c[q]);
    double g_d_phi=2.0*(2*pow(c[q],3)-3.0*pow(c[q],2)+c[q]);
    double dG_dX2=UGC*Temp/(x[q]*(1.0-x[q]));
    double dG_dX_dPhi=(4.0*pow(c[q],3)-6.0*pow(c[q],2)+2.0*c[q])*(WB-WA);
    double dG_dX3=-UGC*Temp/((1.0-x[q])*(1.0-x[q]));
    double dG3_dPhi2_dX=(12*c[q]*c[q]-12.0*c[q]+2.0)*(WB-WA);
    double dG_dX2_dPhi=0.;
    double dG_dPhi2=2.0*(6.0*c[q]*c[q]-6*c[q]+1)*( WB*x[q]+WA*(1.0-x[q]));
    for(unsigned int A=0;A<dofs_per_cell;A++){
      const unsigned int ca=fe_values.get_fe().system_to_component_index(A).first;
      if(ca==0){
	for(unsigned int B=0;B<dofs_per_cell;B++){
	  const unsigned int cb=fe_values.get_fe().system_to_component_index(B).first;
	  if(cb==0){
	    local_matrix(A,B)+=fe_values.shape_value(A,q)*(1/dt)*fe_values.JxW(q);
	    local_matrix(A,B)+=fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*(M_phi/Vm)*dG_dPhi2*fe_values.JxW(q);
	    for(unsigned int i=0;i<dim;i++){
	      local_matrix(A,B)-=M_phi*epsilon*epsilon*fe_values.shape_grad(A,q)[i]*fe_values.shape_grad(B,q)[i]*fe_values.JxW(q);
	    }
	  }
	}
      }
    }

  }

  for(unsigned int q=0;q<n_q_points;q++)  {
    M=0.;
    M=pow(M_gb/M_alpha,4*c[q]*(1-c[q]));
    double dG_dX_dPhi=2.0*(2*c[q]*c[q]*c[q]-3*c[q]*c[q]+c[q])*(WB-WA);
    for(unsigned int A=0;A<dofs_per_cell;A++){
      const unsigned int ca=fe_values.get_fe().system_to_component_index(A).first;
      if(ca==0){
	for(unsigned int B=0;B<dofs_per_cell;B++){
	  const unsigned int cb=fe_values.get_fe().system_to_component_index(B).first;
	  if(cb==1){
	    local_matrix(A,B)+=fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*(M_phi/Vm)*dG_dX_dPhi*fe_values.JxW(q);
	  }
	}
      }
    }
  }

  for(unsigned int q=0;q<dofs_per_cell;q++){
    M=0.;
    M=pow(M_gb/M_alpha,4*c[q]*(1-c[q]));
    double dG_dX_dPhi=2.0*(2*c[q]*c[q]*c[q]-3*c[q]*c[q]+c[q])*(WB-WA);
    for(unsigned int A=0;A<dofs_per_cell;A++){
      const unsigned int ca=fe_values.get_fe().system_to_component_index(A).first;
      if(ca==1){
	for(unsigned int B=0;B<dofs_per_cell;B++){
	  const unsigned int cb=fe_values.get_fe().system_to_component_index(B).first;
	  if(cb==0){
	    for(unsigned int i=0;i<dim;i++){
	      local_matrix(A,B)+=fe_values.shape_grad(A,q)[i]*fe_values.shape_value(B,q)*M*(1.0-x[q])*x[q]*dG_dX_dPhi*c_j[q][i]*fe_values.JxW(q);
	    }
	    for(unsigned int i=0;i<dim;i++){

	      local_matrix(A,B)+=M*x[q]*(1.0-x[q])*2.0*(2.0*c[q]*c[q]*c[q]-3.0*c[q]*c[q]+c[q])*(WB-WA)*fe_values.shape_grad(A,q)[i]*fe_values.shape_grad(B,q)[i]*fe_values.JxW(q);
	    }
	  }
	}
      }
    }
  }

  for(unsigned int q=0;q<n_q_points;q++){
    double Wx=WB*x[q]+WA*(1.0-x[q]);
    double g_phi=c[q]*c[q]*(1-c[q])*(1-c[q]);
    double g_d_phi=2.0*(2*pow(c[q],3)-3.0*pow(c[q],2)+c[q]);
    double dG_dX2=UGC*Temp/(x[q]*(1.0-x[q]));
    double dG_dX_dPhi=(4.0*pow(c[q],3)-6.0*pow(c[q],2)+2.0*c[q])*(WB-WA);
    double dG_dX3=-UGC*Temp/((1.0-x[q])*(1.0-x[q]));
    double dG3_dPhi2_dX=(12*c[q]*c[q]-12.0*c[q]+2.0)*(WB-WA);
    double dG_dX2_dPhi=0.;
    double dG_dPhi2=2.0*(6.0*c[q]*c[q]-6*c[q]+1)*( WB*x[q]+WA*(1.0-x[q]));

    M=0.;
    M=pow(M_gb/M_alpha,4*c[q]*(1-c[q]));
    for(unsigned int A=0;A<dofs_per_cell;A++){
      const unsigned int ca=fe_values.get_fe().system_to_component_index(A).first;
      if(ca==1){
	for(unsigned int B=0;B<dofs_per_cell;B++){
	  const unsigned int cb=fe_values.get_fe().system_to_component_index(B).first;
	  if(cb==1){
	    local_matrix(A,B)+=(1.0/dt)*fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*fe_values.JxW(q);
	    for(unsigned int i=0;i<dim;i++){
	      local_matrix(A,B)-=fe_values.shape_grad(A,q)[i]*fe_values.shape_value(B,q)*M*x[q]*((UGC*Temp/(x[q]*(1.0-x[q])))*x_j[q][i]+(2.0*(2.0*c[q]*c[q]*c[q]-3.0*c[q]*c[q]+c[q])*(WB-WA))*c_j[q][i])*fe_values.JxW(q);
	    }
	    for(unsigned int i=0;i<dim;i++){
	      local_matrix(A,B)+=M*(1.0-x[q])*fe_values.shape_grad(A,q)[i]*fe_values.shape_value(B,q)*((UGC*Temp/(x[q]*(1.0-x[q])))*x_j[q][i]+(2.0*(2.0*c[q]*c[q]*c[q]-3.0*c[q]*c[q]+c[q])*(WB-WA))*c_j[q][i])*fe_values.JxW(q);
	    }
	    for(unsigned int i=0;i<dim;i++){
	      local_matrix(A,B)+=M*(1.0-x[q])*x[q]*fe_values.shape_grad(A,q)[i]*fe_values.shape_value(B,q)*(dG_dX3*x_j[q][i])*fe_values.JxW(q);
	      }
	    for(unsigned int i=0;i<dim;i++){
	      local_matrix(A,B)+=M*(1.0-x[q])*x[q]*fe_values.shape_grad(A,q)[i]*fe_values.shape_grad(B,q)[i]*(UGC*Temp/(x[q]*(1.0-x[q])))*fe_values.JxW(q);
	      }
	  }
	}
      }
    }
  }*/
  
  
}

#endif /* CHEMO_H_ */
