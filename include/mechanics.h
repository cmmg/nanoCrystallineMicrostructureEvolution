//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2011
//authors: rudraa (2011, 2018)
//
#ifndef MECHANICS_H_
#define MECHANICS_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include "deformationMap.h"
#include "variations.cc"
#include "tangentmodulus.cc"

template <int dim>
struct historyVariables{
public:
  historyVariables<dim>():  Fp_previous(dim, dim), Fp_previousIteration(dim, dim), Gamma(n_slip_systems),slipRates(n_slip_systems), CRSS(n_slip_systems), CRSS_iteration(n_slip_systems), Gamma_iteration(n_slip_systems),shear_Stress(n_slip_systems), shear_Stress_iteration(n_slip_systems),Gamma_total(n_slip_systems){}
  Table<1, std::vector<double> >slipRates;
  double eqvstress, equvStrain;
  Table<1, double>CRSS, CRSS_iteration;
  Table<1, double> Gamma,Gamma_iteration,Gamma_total ;
  Table<1, double>shear_Stress, shear_Stress_iteration;
  dealii::FullMatrix<double> Fp_previous, Fp_previousIteration;
};

template<int dim>
void crystal_rotation(Table<2, double>&slip_normal, Table<2, double>&slip_direction,double &theta){
  unsigned int n=dim;
  theta=30.0;
  FullMatrix<double> Rx, Ry, Rz;
  Rx=IdentityMatrix(n);Ry=IdentityMatrix(n); Rz=IdentityMatrix(n);
  double param=theta*PI/180;
  Rx(1,1)=cos(param); Rx(2,2)=cos(param); Rx(1,2)=-sin(param); Rx(2,1)=sin(param);
  Ry(0,0)=cos(param); Rx(2,2)=cos(param); Rx(0,2)=sin(param); Rx(2,0)=-sin(param);
  Rz(0,0)=cos(param); Rx(1,1)=cos(param); Rx(0,1)=-sin(param); Rx(1,0)=sin(param);
  Table<2, double>normals(n_slip_systems,dim), directions(n_slip_systems, dim);
  for(unsigned int i=0;i<n_slip_systems;i++){
    for(unsigned int j=0;j<dim;j++){
      normals[i][j]=slip_normal[i][j];
      directions[i][j]=slip_direction[i][j];
      slip_normal[i][j]=0.;
      slip_direction[i][j]=0.;
    }
  }
  for(unsigned int i=0;i<n_slip_systems;i++){
    for(unsigned int j=0;j<dim;j++){
      for(unsigned int k=0;k<dim;k++){
	for(unsigned int l=0;l<dim;l++){
	  for(unsigned int m=0;m<dim;m++){
	    slip_normal[i][j]+=Rx(j,k)*Ry(k,l)*Rz(l,m)*normals[i][m];
	    slip_direction[i][j]+=Rx(j,k)*Ry(k,l)*Rz(l,m)*directions[i][m];
	  }
	}
      }
    }
    
  }
  

  
  
}

template<class T, int dim>
  void tensor_dyadic(Table<2, T>& A, Table<2, T>& B, Table<4, T>& C){

  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  C[i][j][k][l]=0.;
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  C[i][j][k][l]+=A[i][j]*B[k][l];
}
template<class T, int dim>
  void tensor_contraction(Table<2, T>& A, Table<2, T>& B, T & norm){
  norm=0.;
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      norm+=A[i][j]*B[i][j];
}

template<class T, int dim>
  void tensor_multiply(Table<4, T>& A,Table<2, T>& B, Table<2, T>& C){
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      C[i][j]=0.;
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  C[i][j]+=A[i][j][k][l]*B[k][l];
}


template<class T,int dim>
  void calculate_stress_strain(FullMatrix<T>& F, FullMatrix<T>& Fp,FullMatrix<T>& F_e,FullMatrix<T>& R_cauchy_green, Table<2,T>& E, Table<2,T>& SPK, Table<2,T>& Kirchhoff_stress,Table<2,T>& CauchyStress, Table<2,T>& Piola_stress, Table<4,T>& ElasticModulii){
  FullMatrix<double>Fp_inv(dim, dim),F_inv(dim, dim);
  Fp_inv.invert(Fp);F_inv.invert(F);
  double det_F;
  det_F=F.determinant();
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	F_e(i,j)+=F(i,k)*Fp_inv(k,j);
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	R_cauchy_green(i,j)+=F_e(k,i)*F_e(k,j);
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      E[i][j]+=0.5*(R_cauchy_green[i][j]-(double)(i==j));
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  SPK(i,j)+=ElasticModulii[i][j][k][l]*E[k][l];
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  Kirchhoff_stress[i][l]+=F_e(i,j)*SPK(j,k)*F_e(l,k);
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      CauchyStress[i][j]=(1/det_F)*Kirchhoff_stress[i][j];
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	Piola_stress[i][j]+=Kirchhoff_stress[i][k]*F_inv(j,k);
}

template<class T, int dim>
  void calculate_schmidt_tensors(FullMatrix<T>& Rotation, Table<2,T>& slip_normal, Table<2, T>& slip_direction,std::vector<Table<2, T> >& SchmidtTensor, std::vector<Table<2, T> >& SchmidtTensorUnsym,std::vector<Table<2, T> >&OriginalSchmidt ){
  Table<1, double>n(dim), m(dim);Table<2, double> TT(dim, dim); Table<2, double> TT1(dim, dim);
  std::vector<Table<1, double> > normal, direction;
  Rotation=IdentityMatrix(3);
  for(unsigned int i=0;i<n_slip_systems;i++){
    for(unsigned int A=0;A<dim;A++) for(unsigned int B=0;B<dim;B++){ TT[A][B]=0.;TT1[A][B]=0.;}
    for(unsigned int j=0;j<dim;j++){
      for(unsigned int k=0;k<dim;k++){
	TT[j][k]+=slip_direction[i][j]*slip_normal[i][k];
      }
    }
    OriginalSchmidt.push_back(TT);
  }  
  for(unsigned int i=0;i<n_slip_systems;i++){
    for(unsigned int j=0;j<dim;j++){
      n[j]=0.;m[j]=0.;
    }
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++){
	n[j]+=Rotation(j,k)*slip_normal[i][k];
	m[j]+=Rotation(j,k)*slip_direction[i][k];
      }
    normal.push_back(n); direction.push_back(m);
  }
  //calculate schmidt tensor for every slip system
  for(unsigned int i=0;i<n_slip_systems; i++){
    for(unsigned int A=0;A<dim;A++) for(unsigned int B=0;B<dim;B++){ TT[A][B]=0.;TT1[A][B]=0.;}
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim; k++){
	TT[j][k]=((normal[i][j]*direction[i][k])+(direction[i][j]*normal[i][k]))/2.0;
	TT1[j][k]=(direction[i][j]*normal[i][k]);
      }
    SchmidtTensor.push_back(TT);
    SchmidtTensorUnsym.push_back(TT1);
  }
  
}

template<int dim>
void update_CRSS(Table<1, double>&crss, Vector<double>& gamma_total, Table<1, double>&CRSSprevious,Table<2, double>&slip_normal){
  double h0=self_hardening;
  double q;double a=2.25;
  for(unsigned int i=0;i<n_slip_systems;i++){
    double sum=0.;
    crss[i]=0;
    for(unsigned int j=0;j<n_slip_systems;j++){
      // determine q=1 for coplanar slip systems and 1.4 for others
      double coplanarity=0.;
      coplanarity=std::pow((slip_normal[i][0]-slip_normal[j][0]),2)+std::pow((slip_normal[i][1]-slip_normal[j][1]),2)+std::pow((slip_normal[i][2]-slip_normal[j][2]),2);
      
      if(coplanarity<1e-6)q=1;
      else q=1.4;
      double h_alpha_beta= h0*q*std::pow((1- crss[j]/Ss),a);
      sum+=h_alpha_beta*gamma_total[j];
    }
  crss[i]=(CRSSprevious[i])+sum;
  }
  
}



template<class T, int dim>
  void calculate_shear(Table<2, T>& SPK,std::vector<Table<2, T> >& SchmidtTensorUnsym,Table<1, T>& shearStress,Table<1, T>& crss,std::vector<int>& _ActiveSystems){
  
  for(int I=0;I<n_slip_systems;I++){
    tensor_contraction<double, dim>(SchmidtTensorUnsym[I], SPK, shearStress[I]);
   
   
    if((std::abs(shearStress[I]))-crss[I]>std::abs(1e-3*crss[I])){
      _ActiveSystems.push_back(I); 
    }
  }
}
template<class T, int dim>
  void tensor_triple_contraction(Table<4, T>& EM, Table<2, T>& A, Table<2, T>& B, T & norm){
  norm=0.;
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  norm+=EM[i][j][k][l]*A[i][j]*B[k][l];
      
}
template<class T,int dim>
  void find_polar_decomposition(FullMatrix<T>& F_e, FullMatrix<T>& Rotation){
  FullMatrix<T> F1(dim, dim),invF1(dim, dim), trial(dim, dim);
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      F1(i,j)=F_e(i,j);
  double norm=F1.determinant(); 
  for(unsigned int i=0;i<40;i++){
    invF1.invert(F1);
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++){
	Rotation(j,k)=(F1(j,k)+invF1(k,j))/2.0;
      }
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	F1(j,k)=Rotation(j,k);
  }
}

//mechanics implementation
template <class T, int dim>
  void evaluateStress(unsigned int q, const unsigned int DOF, const Table<1, T>& ULocal, const deformationMap<T, dim>& defMap, unsigned int currentIteration,unsigned int currentIncrement,  std::vector<historyVariables<dim>*>& history, Table<2, double>& Kirchhoff_stress, Table<2, double>& Piola_stress, Table<4,double>& C,  Table<2,double>& invF, FEValues<dim>& fe_values){
 
  
  //update history variables at the start of increment
  if (currentIteration==0){
    history[q]->Fp_previous=history[q]->Fp_previousIteration;
    history[q]->CRSS=history[q]->CRSS_iteration;
    history[q]->Gamma=history[q]->Gamma_iteration;
    history[q]->shear_Stress=history[q]->shear_Stress_iteration;
    for(int i=0;i<n_slip_systems;i++){
      history[q]->Gamma_total[i]+=history[q]->Gamma[i];
    }
    
    if(q==0){
      std::cout<<"\nCRSS\n";
      for(int i=0;i<n_slip_systems;i++)std::cout<<history[q]->CRSS[i]<<" ";
      std::cout<<"\nShearStress\n";
      for(int i=0;i<n_slip_systems;i++)std::cout<<history[q]->shear_Stress[i]<<" ";
      std::cout<<"\nGamma\n";
      for(int i=0;i<n_slip_systems;i++)std::cout<<history[q]->Gamma[i]<<" ";
      std::cout<<"\nGammaTotal\n";
      for(int i=0;i<n_slip_systems;i++)std::cout<<history[q]->Gamma_total[i]<<" ";
      std::cout<<"\n\n";
    }
  }
  
  Table<2, double> slip_normal(n_slip_systems, dim);
  Table<2, double> slip_direction(n_slip_systems, dim);
  /*declare slip plane normals and slip directions*/
 
  slip_normal[0][0]=0.577;    slip_normal[0][1]=0.577;    slip_normal[0][2]=0.577; slip_direction[0][0]=0.707;  slip_direction[0][1]=-0.707;   slip_direction[0][2]=0.0;
  /*slip_normal[1][0]=0.577;    slip_normal[1][1]=0.577;    slip_normal[1][2]=0.577; slip_direction[1][0]=-0.707; slip_direction[1][1]=0.0;      slip_direction[1][2]=0.707;
  /*slip_normal[2][0]=0.577;    slip_normal[2][1]=0.577;    slip_normal[2][2]=0.577; slip_direction[2][0]=0.0;    slip_direction[2][1]=0.707;    slip_direction[2][2]=-0.707;
  /*slip_normal[3][0]=-0.577;   slip_normal[3][1]=0.577;    slip_normal[3][2]=0.577; slip_direction[3][0]=0.707;  slip_direction[3][1]=0.0;      slip_direction[3][2]=0.707;
  slip_normal[4][0]=-0.577;   slip_normal[4][1]=0.577;    slip_normal[4][2]=0.577; slip_direction[4][0]=-0.707; slip_direction[4][1]=-0.707;   slip_direction[4][2]=0.0;
  slip_normal[5][0]=-0.577;   slip_normal[5][1]=0.577;    slip_normal[5][2]=0.577; slip_direction[5][0]=0.0;    slip_direction[5][1]=0.707;    slip_direction[5][2]=-0.707;
    /*slip_normal[6][0]=0.577;    slip_normal[6][1]=-0.577;   slip_normal[6][2]=0.577; slip_direction[6][0]=-0.707; slip_direction[6][1]=0.0;      slip_direction[6][2]=0.707;
  slip_normal[7][0]=0.577;    slip_normal[7][1]=-0.577;   slip_normal[7][2]=0.577; slip_direction[7][0]=0.0;    slip_direction[7][1]=-0.707;   slip_direction[7][2]=-0.707;
  slip_normal[8][0]=0.577;    slip_normal[8][1]=-0.577;   slip_normal[8][2]=0.577; slip_direction[8][0]=0.707;  slip_direction[8][1]=0.707;    slip_direction[8][2]=0.0;
  slip_normal[9][0]=-0.577;   slip_normal[9][1]=-0.577;   slip_normal[9][2]=0.577; slip_direction[9][0]=-0.707; slip_direction[9][1]=0.707;    slip_direction[9][2]=0.0;
  slip_normal[10][0]=-0.577;  slip_normal[10][1]=-0.577;  slip_normal[10][2]=0.577;slip_direction[10][0]=0.707; slip_direction[10][1]=0.0;     slip_direction[10][2]=0.707;
  slip_normal[11][0]=-0.577;  slip_normal[11][1]=-0.577;  slip_normal[11][2]=0.577;slip_direction[11][0]=0.0;   slip_direction[11][1]=-0.707;  slip_direction[11][2]=-0.707;*/

 
 
  
  FullMatrix<double> F(dim, dim),F_inv(dim, dim),Fp(dim, dim),Fp_inv(dim, dim),F_e(dim, dim),Fe_inv(dim, dim), Rotation(dim, dim), R_cauchy_green(dim, dim);
  
  Table<2, double> Stress(dim, dim), E(dim, dim), SPK(dim, dim), CauchyStress(dim, dim), B_beta(dim, dim);
  Table<4, double>ElasticModulii(dim, dim, dim, dim), _ElasticModulii(dim, dim, dim, dim);
  std::vector<Table<2,double> > SchmidtTensor, SchmidtTensorUnsym, OriginalSchmidt;
  Table<1, double>   shearStress(n_slip_systems), crss(n_slip_systems);
  //std::vector<unsigned int>ActiveSystems;
  Vector<int>ActiveSystems,activesystemoriginal;
  std::vector<int>_ActiveSystems, _activesystemoriginal;
  unsigned int activesystem=0;
  
  double nu=PoissonsRatio;double Y=elasticModulus;
  double lambda=(nu*Y)/((1+nu)*(1-2*nu)), mu=Y/(2*(1+nu));
  double kappa= lambda+(2.0/3.0)*mu;
  double det_F=0.;
  double theta=0.0;
  Table<2, double> Identity(dim, dim);
  Table<2, double> seedPoint(n_seed_points, dim);
  std::srand(5);
  for(unsigned int i=0;i<n_seed_points;i++){
    for(unsigned int j=0;j<dim;j++){
      seedPoint[i][j]=((double)(std::rand() %100)/100)-0.5;

    }
  }
  Point<dim> quadPoint=fe_values.quadrature_point(q);
  //assign grain id to the quadrature point
  Table<1, double>distance(n_seed_points);
  for(unsigned int i=0;i<n_seed_points;i++){
    distance[i]=0.;
    for(unsigned int j=0;j<dim;j++){
      distance[i]+=std::pow((quadPoint[j]-seedPoint[i][j]),2);
    }
  }
  double mind=distance[0];unsigned int grain_id=0;
  for(unsigned int i=1;i<n_seed_points;i++){
    if(distance[i]<mind)grain_id=i;
  }
  switch(grain_id){
  case 0: theta=22;break;
  case 1: theta=44;break;
  case 2: theta=66;break;
  case 3: theta=88;break;
  }
  // crystal_rotation<dim>(slip_normal, slip_direction,theta);
  //define ElasticModulii
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++){
	for(unsigned int l=0;l<dim;l++){
	  ElasticModulii[i][j][k][l]=lambda*(i==j)*(k==l) + mu*((i==k)*(j==l) + (i==l)*(j==k));
	}

      }
  
  /*calculate deformation gradient*/
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      F(i,j)=defMap.F[q][i][j];
    }
  }
  for(unsigned int i=0;i<n_slip_systems;i++){
    crss[i]=history[q]->CRSS[i];
  }
  Fp=history[q]->Fp_previous;
  F_inv.invert(F);
  Fp_inv.invert(Fp);F_inv.invert(F);
  calculate_stress_strain<double, dim>(F,history[q]->Fp_previous,F_e,R_cauchy_green,E,SPK,Kirchhoff_stress, CauchyStress, Piola_stress,ElasticModulii);
  find_polar_decomposition<double, dim>(F_e, Rotation);  //polar decomposition of Fe
  calculate_schmidt_tensors<double, dim>(Rotation, slip_normal, slip_direction, SchmidtTensor, SchmidtTensorUnsym, OriginalSchmidt); 
  calculate_shear<double, dim>( SPK, SchmidtTensorUnsym, shearStress, crss, _ActiveSystems);
  ActiveSystems.reinit(_ActiveSystems.size());
  for(unsigned int I=0;I<_ActiveSystems.size();I++){
    ActiveSystems[I]=_ActiveSystems[I];
  }
  if(currentIteration==0 && q==0){
    std::cout<<"shear at begining\n";
    for(int i=0;i<n_slip_systems;i++)std::cout<<shearStress[i]<<" ";
    std::cout<<"\n\n";
  }
  unsigned int size=ActiveSystems.size(); 
  unsigned int cntr=0;
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int I=0;I<dim;I++)
      for(unsigned int j=0;j<dim;j++)
	for(unsigned int J=0;J<dim;J++)
	  for(unsigned int k=0;k<dim;k++)
	    for(unsigned int K=0;K<dim;K++)
	      for(unsigned int l=0;l<dim;l++)
		for(unsigned int L=0;L<dim;L++)
		  C[i][j][k][l]+=F(i,I)*F(j,J)*F(k,K)*F(l,L)*ElasticModulii[I][J][K][L];


 
  unsigned int whilecounter=0;
   Vector<double> gamma_total(n_slip_systems);
  if(size>0){
   
    FullMatrix<double> active_system_matrix, active_system_inv;
    Vector<double> active_system_rhs;
    Vector<double> gamma;
   
    std::vector<int> tempv;
    Table<2, double>B_mult_SPK(dim, dim), C_B_beta(dim, dim), RCB(dim, dim), temp(dim, dim), B_beta(dim, dim);
    Table<4,double>variation_Fp(dim, dim, dim, dim),variation_Fe(dim, dim, dim, dim), variation_E(dim, dim, dim, dim),variation_SPK(dim, dim, dim, dim), ElastoPlasticModulii(dim, dim, dim, dim);
    Table<3, double>variation_crss(n_slip_systems, dim, dim);
    for(unsigned int i=0;i<dim;i++)
      for(unsigned int j=0;j<dim;j++)
	for(unsigned int k=0;k<dim;k++)
	  for(unsigned int l=0;l<dim;l++){
	    variation_Fp[i][j][k][l]=0.;
	    variation_Fe[i][j][k][l]=0.;
	    variation_E[i][j][k][l]=0.;
	    variation_SPK[i][j][k][l]=0.;
	    ElastoPlasticModulii[i][j][k][l]=0.;
	  }
    for(int i=0;i<n_slip_systems;i++)
      for(int j=0;j<dim;j++)
	for(int k=0;k<dim;k++)
	  variation_crss[i][j][k]=0.;

    
    unsigned int alpha, beta, flag=1, flag1=1;
    // put a loop to check for yield surface criterion
    while(flag1==1) {
      flag=1;
      while(flag==1){//loop for calculating gamma---- removing negative slips
	active_system_matrix.reinit(size,size); active_system_inv.reinit(size, size);
	active_system_rhs.reinit(size); gamma.reinit(size);
	active_system_matrix=0.;active_system_inv=0.;
	active_system_rhs=0.; gamma=0.;
	for(unsigned int i=0;i<size;i++){
	  active_system_rhs[i]=0.;
	  gamma[i]=0.;
	  for(unsigned int j=0;j<size;j++){
	    active_system_matrix(i,j)=0.;
	    active_system_inv(i,j)=0.;
	  }
	}
	for(unsigned int I=0;I<size;I++){
	  alpha=ActiveSystems[I];
	  for(unsigned int J=0;J<size;J++){
	    beta=ActiveSystems[J];
	    
	    for(unsigned int i=0;i<dim;i++)
	      for(unsigned int j=0;j<dim;j++){
		B_beta[i][j]=0.;
		C_B_beta[i][j]=0.;
		RCB[i][j]=0.;
		temp[i][j]=0.;
		B_mult_SPK[i][j]=0.;
	      }
	    //calculate B_beta matrix
	    for(unsigned int i=0;i<dim;i++)
	      for(unsigned int j=0;j<dim;j++)
		for(unsigned int k=0;k<dim;k++)
		  B_beta[i][j]+=((R_cauchy_green(i,k)*OriginalSchmidt[beta][k][j])+(OriginalSchmidt[beta][k][i]*R_cauchy_green(k,j)))/2.0; //here schmidt tensor original or unsymmetric
	    for(unsigned int i=0;i<dim;i++)
	      for(unsigned int j=0;j<dim;j++)
		for(unsigned int k=0;k<dim;k++)
		B_mult_SPK[i][j]+=B_beta[i][k]*SPK(k,j);
	    for(unsigned int i=0;i<dim;i++)
	      for(unsigned int j=0;j<dim;j++)
		B_mult_SPK[i][j]*=2.0;
	    for(unsigned int i=0;i<dim;i++)
	      for(unsigned int j=0;j<dim;j++)
		for(unsigned int k=0;k<dim;k++)
		  for(unsigned int l=0;l<dim;l++)
		    C_B_beta[i][j]+=ElasticModulii[i][j][k][l]*B_beta[k][l]; 
	    for(unsigned int i=0;i<dim;i++)
	      for(unsigned int j=0;j<dim;j++)
		for(unsigned int k=0;k<dim;k++)
		  RCB[i][j]+=R_cauchy_green(i,k)*C_B_beta[k][j];
	    for(unsigned int i=0;i<dim;i++)
	      for(unsigned int j=0;j<dim;j++)
		temp[i][j]=RCB[i][j]+B_mult_SPK[i][j];
	    tensor_contraction<double, dim>(OriginalSchmidt[alpha], temp, active_system_matrix(I,J));
	    double hardening_alpha_beta=0.;
	    active_system_matrix(I,J)=active_system_matrix(I,J)*sign(shearStress[alpha])*sign(shearStress[beta]);
	    
	  }
	}
	
	
	/*if(currentIteration==4 && currentIncrement==2){
	  std::cout<<"\nactive system matrix";
	  for(int i=0;i<size;i++){
	    for(int j=0;j<size;j++){
	      std::cout<<active_system_matrix(i,j)<<" ";
	    }std::cout<<"\n";
	  }
	  std::cout<<"whilecounter"<<whilecounter<<"quadpoint"<<q;
	  }*/
	
	//setup active_system_rhs
	for(unsigned int I=0;I<size;I++){
	  int alpha=ActiveSystems[I];
	  active_system_rhs[I]=((std::abs(shearStress[alpha]))-crss[alpha]);
	  if(active_system_rhs[I]<0){std::cout<<"Error:: matrix rhs cant be negative";exit(-1);}
	}
	
	active_system_inv.invert(active_system_matrix);
        
	cntr=0;tempv.resize(0);
	for(unsigned int I=0;I<size;I++){
	  for(unsigned int J=0;J<size;J++){
	    gamma[I]+=active_system_inv(I,J)*active_system_rhs[J];
	  }
	  if(gamma[I]> 1e-6) {
	    cntr++;tempv.push_back(ActiveSystems[I]);
	    if(currentIteration==4/*gamma[I]>1.0e0*/) { std::cout << crss[I] << " " << shearStress[I] << " " << gamma[I] << " ";}
	    if(currentIteration==4/*gamma[I]>1.0e0*/) {for(unsigned int J=0;J<size;J++){ std::cout << active_system_inv(I,J) << " ";}}
	    if(currentIteration==4/*gamma[I]>1.0e0*/) {for(unsigned int J=0;J<size;J++){ std::cout << active_system_rhs(J) << " ";}}
	    if(currentIteration==4/*gamma[I]>1.0e0*/) {std::cout << "\n";}  
	  }//cntr counts no of +ve gamma values
	}
	
	if(cntr<size && cntr!=0){ //cntr==size-> no negative gamma values, and cntr=0-> all negative gamma
	  size=cntr;
	  ActiveSystems.reinit(0);
	  ActiveSystems.reinit(size);
	  for(unsigned int i=0;i<size;i++)ActiveSystems[i]=tempv[i];
	}
	else{flag=0;}
	
	//inner while loop end
      }
      

      
      //update Fp plastic deformation gradient (approach by Anand & Kothari) only if any slip system found active
      if(cntr>0){//if any slip system active
        
	Table<2, double> updateFp(dim, dim);
	for(unsigned int i=0;i<dim;i++)
	  for(unsigned int j=0;j<dim;j++)
	    updateFp[i][j]=0.;
	
	//add gamma to gamma old
	for(unsigned int i=0;i<n_slip_systems;i++){
	  for(unsigned int I=0;I<size;I++){
	    if(ActiveSystems[I]==i){gamma_total[i]+=gamma[I];break;}
	  }
	  
	  
	}
         Table<3,double>variation_gamma(size, dim, dim);
	 for(unsigned int i=0;i<size;i++)
	   for(unsigned int j=0;j<dim;j++)
	     for(unsigned int k=0;k<dim;k++)
	       variation_gamma[i][j][k]=0.;
	//calculate all variation Fp_inv, Fe and others
        
 	calculate_variation_Fe<dim>(variation_Fe, F, Fp, variation_Fp);
	calculate_variation_E<dim>(variation_E, variation_Fe, F_e);
	calculate_variation_SPK<dim>(variation_SPK, variation_E, ElasticModulii);
	
	double _gamma=0.;
	for(unsigned int I=0;I<n_slip_systems;I++){//previously loop over size of active systems
	  _gamma=0.;
	  _gamma=gamma_total[I]*sign(shearStress[alpha]);
	  for(unsigned int i=0;i<dim;i++)
	    for(unsigned int j=0;j<dim;j++){
	      updateFp[i][j]+=_gamma*OriginalSchmidt[I][i][j];
	    }
	}
	
	for(unsigned int i=0;i<dim;i++)
	  for(unsigned int j=0;j<dim;j++){
	    updateFp[i][j]+=(double)(i==j);
	  }
	
	calculate_variation_Fp<dim>(variation_Fp,variation_E,updateFp,  ActiveSystems, OriginalSchmidt, Fp, SPK, variation_SPK ,R_cauchy_green, gamma_total, active_system_matrix, active_system_rhs, shearStress, crss, ElasticModulii, variation_crss, variation_gamma);
	update_CRSS<dim>(crss, gamma_total,history[q]->CRSS,slip_normal);
	//calculate variational derivative of CRSS
	calculate_variation_crss<dim>(variation_crss, ActiveSystems, variation_gamma,  crss,gamma_total,slip_normal);
	Fp.reinit(dim, dim);
	for(unsigned int i=0;i<dim;i++){
	  for(unsigned int j=0;j<dim;j++){
	    for(unsigned int k=0;k<dim;k++){
	      Fp(i,j)+=updateFp[i][k]*(history[q]->Fp_previous(k,j));
	    }	
	  } 
	}
	
	
	double det_Fp;
	det_Fp=Fp.determinant();
	//to ensure that det(Fp)=1
	for(unsigned int i=0;i<dim;i++)
	  for(unsigned int j=0;j<dim;j++)
	    Fp(i,j)=Fp(i,j)*std::pow(det_Fp,-1./3.);
	//calculate Fe and repeat all procedures
	F_e.reinit(dim,dim); R_cauchy_green.reinit(dim,dim); E.reinit(dim, dim); SPK.reinit(dim, dim);Kirchhoff_stress.reinit(dim, dim);
	CauchyStress.reinit(dim, dim); Piola_stress.reinit(dim, dim);Rotation.reinit(dim, dim);
	Fp_inv.invert(Fp);
	calculate_stress_strain<double, dim>(F,Fp,F_e,R_cauchy_green,E,SPK,Kirchhoff_stress, CauchyStress, Piola_stress,ElasticModulii);
	find_polar_decomposition<double, dim>(F_e, Rotation);
	SchmidtTensor.resize(0);SchmidtTensorUnsym.resize(0); OriginalSchmidt.resize(0);
	calculate_schmidt_tensors<double, dim>(Rotation, slip_normal, slip_direction, SchmidtTensor, SchmidtTensorUnsym, OriginalSchmidt);
	_ActiveSystems.resize(0);

	calculate_shear<double, dim>(SPK, SchmidtTensorUnsym, shearStress, crss, _ActiveSystems);
	if(_ActiveSystems.size()>0){
	  ActiveSystems.reinit(_ActiveSystems.size());
	  for(unsigned int I=0;I<_ActiveSystems.size();I++){
	    ActiveSystems[I]=_ActiveSystems[I];
	  }
	  size=ActiveSystems.size();
	  
	}
	else{
	  flag1=0;
	  //std::cout<<"control here";exit(-1);
	  elastoplastic_tangent<dim>(F_e, Fp, F, SPK, variation_Fe, variation_Fp, variation_SPK, ElastoPlasticModulii,currentIncrement, currentIteration);
	  for(unsigned int i=0;i<dim;i++)
	    for(unsigned int j=0;j<dim;j++)
	      for(unsigned int k=0;k<dim;k++)
		for(unsigned int l=0;l<dim;l++)
		  C[i][j][k][l]=0.;
	  
	  for(unsigned int j=0;j<dim;j++)
	    for(unsigned int J=0;J<dim;J++)
	      for(unsigned int l=0;l<dim;l++)
		for(unsigned int L=0;L<dim;L++)
		  for(unsigned int i=0;i<dim;i++)
		    for(unsigned int k=0;k<dim;k++)
		      C[i][j][k][l]+=F(j,J)*F(l,L)*ElastoPlasticModulii[i][J][k][L];
	  
	}
	
      }
      
      
      //if all gamma negative then move so outer while ends
      else{flag1=0;}
      whilecounter++;
      
      
      //outer while    for consistency condition
    }
   
    //plasticity if ends
  }
  
  
  for(unsigned int i=0;i<n_slip_systems;i++){
    history[q]->CRSS_iteration[i]=crss[i];
    history[q]->Gamma_iteration[i]=gamma_total[i];
    history[q]->shear_Stress_iteration[i]=shearStress[i];
   
  }
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++){
      invF[i][j]=F_inv(i,j);
      history[q]->Fp_previousIteration(i,j)=Fp(i,j);
    }
}


//mechanics residual implementation
template <int dim>
void residualForMechanics(FEValues<dim>& fe_values, unsigned int DOF, Table<1, double >& ULocal, Table<1, double>& ULocalConv, deformationMap<double, dim>& defMap, unsigned int currentIteration,unsigned int currentIncrement,  std::vector<historyVariables<dim>* >& history, Vector<double>& RLocal, FullMatrix<double>& KLocal){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //quadrature loop
  for (unsigned int q=0; q<n_q_points; q++){
    
    //evaluate stress
    Table<2, double> Kirchhoff_stress(dim, dim), Piola_stress(dim, dim), invF(dim, dim);
    Table<4,double> C(dim,dim,dim,dim);
    evaluateStress<double, dim>(q, DOF, ULocal, defMap, currentIteration, currentIncrement, history, Kirchhoff_stress, Piola_stress, C, invF, fe_values);
    
    //evaluate Residual
    for (unsigned int A=0; A<dofs_per_cell; A++) {
      const unsigned int i = fe_values.get_fe().system_to_component_index(A).first - DOF;
      if (i>=0 && i<dim){
	// R = Grad(w)*P
	for (unsigned int j = 0; j<dim; j++){
	  for (unsigned int J = 0; J<dim; J++){  
	    RLocal(A) += -fe_values.shape_grad(A, q)[J]*invF[J][j]*(Piola_stress[i][J])*fe_values.JxW(q);
	  }
	  //evaluate tangent
	  for (unsigned int B=0; B<dofs_per_cell; B++) {   
	    const unsigned int k = fe_values.get_fe().system_to_component_index(B).first - DOF;
	    if (k>=0 && k<dim){
	      //K
	      for (unsigned int l = 0; l<dim; l++){
		for (unsigned int J = 0; J<dim; J++){
		  for (unsigned int L = 0; L<dim; L++){
		    KLocal(A,B)+=fe_values.shape_grad(A,q)[J]*invF[J][j]*C[i][j][k][l]*fe_values.shape_grad(B,q)[L]*invF[L][l]*fe_values.JxW(q);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    //add geomteric stiffness term here
    for(unsigned int A=0;A<dofs_per_cell;A++){
      for(unsigned int j=0;j<dim;j++){
	for(unsigned int J=0;J<dim;J++){
	  for(unsigned int B=0;B<dofs_per_cell;B++){
	    for(unsigned int k=0;k<dim;k++){
	      for(unsigned int K=0;K<dim;K++){
		KLocal(A,B)+=fe_values.shape_grad(A,q)[J]*invF[J][j]*Kirchhoff_stress[j][k]*fe_values.shape_grad(B,q)[K]*invF[K][k]*fe_values.JxW(q);
	      }
	    }
	  } 
	}
      }
    }
    
    
    //end quadrature loop
    
  }
}

#endif /* MECHANICS_H_ */


