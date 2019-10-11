//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2012
//authors: rudraa (2012, 2018)
//

#ifndef CHEMO_H_
#define CHEMO_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"


template<int  dim>
struct historyVariables{
public:
  historyVariables<dim>():grainID(0.0){}
  Sacado::Fad::DFad<double> grainID;
  
};


template<int dim>
void evaluatePhiValue(unsigned int q, Table<1, Sacado::Fad::DFad<double> >&c, Table<1, Sacado::Fad::DFad<double> >&c_conv,Table<2,Sacado::Fad::DFad<double> >&c_j ,FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocalConv,std::vector<historyVariables<dim>* >&history, unsigned int currentIteration){

  unsigned int dofs_per_cell=fe_values.dofs_per_cell;

  for(unsigned int i=0;i<dofs_per_cell;i++){
    unsigned int ck=fe_values.get_fe().system_to_component_index(i).first;
    if(ck<n_diff_grains){
      c[ck]+=fe_values.shape_value(i,q)*ULocal[i];
      c_conv[ck]+=fe_values.shape_value(i,q)*ULocalConv[i];

      for(unsigned int j=0;j<dim;j++){
	c_j[ck][j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];
      }
    }
  }
  
  if(currentIteration==0){
    unsigned int ctr=0, id=0;
    
    for(unsigned int i=0;i<n_diff_grains;i++){
      if(c_conv[i].val()>0.9){id=i+1;break;}
      else {ctr++;}
    }
    
    if(ctr==n_diff_grains){history[q]->grainID=n_diff_grains+10;}
    else {history[q]->grainID=id;}
    //if ends
  }
  
  // function definition ends 
}


template<int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, FullMatrix<double >&local_matrix,std::vector<historyVariables<dim>* >&history, unsigned int currentIteration){
  unsigned int dofs_per_cell=fe_values.dofs_per_cell;
  unsigned int n_q_points=fe_values.n_quadrature_points;
  double Kappa[]=InterfaceEnergyParameter;
  for(unsigned int q=0;q<n_q_points;q++){
    Table<1, Sacado::Fad::DFad<double> >c(n_diff_grains), c_conv(n_diff_grains);// phase variable at each quadrature point
    Table<2, Sacado::Fad::DFad<double> >c_j(n_diff_grains, dim);
    for(int i=0;i<n_diff_grains;i++){c[i]=0.; c_conv[i]=0.;}
    for(int i=0;i<n_diff_grains;i++)for(int j=0;j<dim;j++)c_j[i][j]=0.;
    evaluatePhiValue<dim>( q, c, c_conv, c_j ,fe_values, DOF, fe_face_values,cell,dt, ULocal, ULocalConv,history, currentIteration);//calculate phi values and gradient at this quadrature point

    for(unsigned int i=0;i<dofs_per_cell;i++){
      unsigned int ck=fe_values.get_fe().system_to_component_index(i).first;
      Sacado::Fad::DFad<double> sq_sum=0.;
      for(int I=0;I<n_diff_grains;I++){if(ck==I)continue; else sq_sum+=c[I]*c[I];}
      R[i]+=(1.0/dt)*fe_values.shape_value(i,q)*(c[ck]-c_conv[ck])*fe_values.JxW(q);
      R[i]-=L*fe_values.shape_value(i,q)*(c[ck]-pow(c[ck],3)-2.0*alpha1*c[ck]*sq_sum)*fe_values.JxW(q);

      for(unsigned int j=0;j<dim;j++){
	Sacado::Fad::DFad<double> Kjj=Kappa[j];
	Sacado::Fad::DFad<double> kc_j=c_j[ck][j]*Kjj;
	R[i]+=L*fe_values.shape_grad(i,q)[j]*kc_j*fe_values.JxW(q);
      }
    }
    
  }
  
}


//Chemistry residual implementation
/*template <int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, dealii::Table<2, Sacado::Fad::DFad<double> >& c_conv, FullMatrix<double >&local_matrix){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  Table<2, Sacado::Fad::DFad<double> >LocalMatrix(dofs_per_cell, dofs_per_cell);
  //evaluate gradients 
  
  dealii::Table<2, Sacado::Fad::DFad<double> > c(n_q_points,n_diff_grains);         
  dealii::Table<3, Sacado::Fad::DFad<double> > c_j(n_q_points,n_diff_grains, dim);
  for (unsigned int q=0;q<n_q_points;q++) {
    for(unsigned int p=0;p<n_diff_grains;p++){
      c[q][p]=0.0;c_conv[q][p]=0.0;
      for(unsigned int j=0;j<dim;j++){
	c_j[q][p][j]=0.0;
      }
    }
    for(unsigned int i=0;i<dofs_per_cell;i++){ 
      unsigned int ck=fe_values.get_fe().system_to_component_index(i).first- DOF;
      c[q][ck]+=fe_values.shape_value(i,q)*ULocal[i];
      c_conv[q][ck]+=fe_values.shape_value(i,q)*ULocalConv[i];
      for(unsigned int j=0;j<dim;j++){
	c_j[q][ck][j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];
      }
    }
    
  }
  double Kappa[]=InterfaceEnergyParameter;
  
  for(unsigned int i=0; i<dofs_per_cell; i++){
    unsigned int ck=fe_values.get_fe().system_to_component_index(i).first - DOF;
    for(unsigned int q=0;q<n_q_points;q++){
      Sacado::Fad::DFad<double> sq_sum=0.0;
      for(int p=0;p<n_diff_grains;p++){ if(p==ck)continue; else sq_sum+=(c_conv[q][p])*(c_conv[q][p]); }
      R[i]+=(1.0/dt)*fe_values.shape_value(i,q)*(c[q][ck]-c_conv[q][ck])*fe_values.JxW(q);
      R[i]-=L*fe_values.shape_value(i,q)*(c[q][ck]-pow(c[q][ck],3)-2.0*alpha1*c[q][ck]*sq_sum )*fe_values.JxW(q);
      
      for(unsigned int j=0;j<dim;j++){
	Sacado::Fad::DFad<double> Kjj= Kappa[j];
	Sacado::Fad::DFad<double> kc_j= c_j[q][ck][j]*Kjj; // Kjj*C_j
	R[i]+=L*fe_values.shape_grad(i,q)[j]*kc_j*fe_values.JxW(q);
      }
      
    }
  }
  



  
}*/


#endif /* CHEMO_H_ */
