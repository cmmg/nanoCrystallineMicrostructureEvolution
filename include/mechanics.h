//new
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2011
//authors: rudraa (2011, 2018)
//
#ifndef MECHANICS_H_
#define MECHANICS_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include "deformationMap.h"


template <int dim>
struct historyVariables{
public:
  historyVariables<dim>():
  alpha(0.0), beta (dim+1), betaIteration(dim+1), Ep(dim+1), EpIteration(dim+1), stress(0) ,elasStrain12(0.0){}

  //using std:::map to store time history variables
  double alpha, alphaIteration;
  double eqvStress, equvStrain;
  double elasStrain11, elasStrain12;
  dealii::Table<1, double > beta, betaIteration;
  double stress;
  dealii::Table<1, double > Ep, EpIteration;
};

template<int dim>
void assign_grain_id(Table<1, double>& angle,unsigned int & currentIncrement){
  std::srand(5.55);
  for(unsigned int i=0;i<dim;i++){
    if(i%2)angle[i]=0.0;
    else angle[i]=90.0;
  }
  
}

template<int dim>
void crystalRotation(Table<1, double>& phi, Table<1, double>& angle, std::vector<FullMatrix<double> >& rotationMatrices){
  double param=0.;
  FullMatrix<double> rotation(dim, dim);
  for(unsigned int I=0;I<n_diff_grains;I++){
    param=angle[I]*PI/180.0;
    rotation(0,0)=cos(param); rotation(0,1)=-sin(param); rotation(1,0)=sin(param);  rotation(1,1)=cos(param);
    if(dim==3){
      rotation(2,2)=1.0;
    }
    rotationMatrices.push_back(rotation);
    
  }
 
}

template<int dim>
void computeRotatedModulii(Table<4, double>&ElasticModulus, std::vector<FullMatrix<double> >&rotationMatrices, std::vector<Table<4, double> >&A_phi){

  for(unsigned int N=0;N<n_diff_grains;N++){
    Table<4, double>temp(dim, dim, dim, dim);
    FullMatrix<double> rotation(dim, dim);
    for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)rotation[i][j]=rotationMatrices[N][i][j];
    for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)for(unsigned int k=0;k<dim;k++)for(unsigned int l=0;l<dim;l++)temp[i][j][k][l]=0.;
    for(unsigned int I=0;I<dim;I++){
      for(unsigned int i=0;i<dim;i++){
	for(unsigned int J=0;J<dim;J++){
	  for(unsigned int j=0;j<dim;j++){
	    for(unsigned int K=0;K<dim;K++){
	      for(unsigned int k=0;k<dim;k++){
		for(unsigned int L=0;L<dim;L++){
		  for(unsigned int l=0;l<dim;l++){
		    temp[i][j][k][l]+=rotation[i][I]*rotation[j][J]*rotation[k][K]*rotation[l][L]*ElasticModulus[I][J][K][L];
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    A_phi.push_back(temp);
  }
}


//mechanics implementation
template <class T, int dim>
  void evaluateStress(unsigned int q, FEValues<dim>& fe_values, const unsigned int DOF, const Table<1, T>& ULocal,const Table<1, T>& ULocalConv, const deformationMap<T, dim>& defMap, unsigned int currentIteration,  std::vector<historyVariables<dim>*>& history,Table<1, double>&grainAngle,Table<1, double>&phi, Table<2, double>& secondPiola, Table<2, double>& piolaStress , Table<2, double>& cauchyStress , Table<4, double>&ElasticModulus ,Table<4, double>&ElasticTangentModulus, Table<2, double>& largeStrain,Table<2, double>& smallStrain ,FullMatrix<double>&Rotation,Table<2, double>&defGrad, double fractionalTime, double & freeEnergyMech,std::vector<double>& dF, std::vector<double>& dE, std::vector<double>& dF_dE, unsigned int & currentIncrement,Table<1, double>&h_phi, std::vector<Table<4, double> >&A_phi){ 


  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  //update history variables at the start of increment
  if (currentIteration==0){
    history[q]->alpha=history[q]->alphaIteration;
    history[q]->beta=history[q]->betaIteration;
    history[q]->Ep=history[q]->EpIteration;
  }
  //Table<1,double> grainOrientationAngle(n_diff_grains); phi(n_diff_grains);
  //calculate OrderParameter at quadrature point
  std::vector<FullMatrix<double> >rotationMatrices;
  
  for(unsigned int i=0;i<n_diff_grains;i++){grainAngle[i]=0.; phi[i]=0.;}
  for(unsigned int i=0;i<dofs_per_cell;i++){
    int ci=fe_values.get_fe().system_to_component_index(i).first-dim;
    if(ci>=0 && ci<n_diff_grains){
      phi[ci]+=fe_values.shape_value(i,q)*ULocal[i];
    }
  }

  
  double h_total=0;
  for(unsigned int i=0;i<n_diff_grains;i++)h_phi[i]=0.;
  for(unsigned int i=0;i<n_diff_grains;i++){
    h_phi[i]=(1.0-cos(PI*phi[i]))/2.0;
    h_total+=h_phi[i];
  }

  assign_grain_id<dim>(grainAngle,currentIncrement);
  crystalRotation<dim>(phi, grainAngle, rotationMatrices);
  
  //material properties
  Table<2,double> Fe (dim, dim);
  Table<2,double> E (dim, dim);
  Table<2,double> gradU (dim, dim);
  
  //Lame parameters
  

  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      gradU[i][j]=defMap.F[q][i][j] - (double)(i==j);
      defGrad[i][j]=defMap.F[q][i][j];
    }
  }
  
  FullMatrix<double> RCG(dim, dim);
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      for(unsigned int k=0;k<dim;k++){
	RCG(i,j)+=defGrad[k][i]*defGrad[k][j];
      }
    }
  }
  
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      largeStrain[i][j]=(0.5)*(RCG(i,j)-(double)(i==j));
      smallStrain[i][j]=(0.5)*(gradU[i][j]+gradU[j][i]);
      piolaStress[i][j]=0.0;
      secondPiola[i][j]=0.0;
      cauchyStress[i][j]=0.0;
    }
  }
  
  if(dim==2){
    double nu12=0.253; double nu21=0.253; double mu=500.0;
    ElasticModulus[0][0][0][0]=alpha1/(1-nu12*nu21);  ElasticModulus[0][0][1][1]=alpha1*nu21/(1-nu12*nu21);
    ElasticModulus[1][1][0][0]=beta1*nu12/(1-nu12*nu21);  ElasticModulus[1][1][1][1]=beta1/(1-nu12*nu21);
    ElasticModulus[0][1][0][1]= mu/2.;   ElasticModulus[0][1][1][0]=mu/2.;   ElasticModulus[1][0][0][1]=mu/2.;   ElasticModulus[1][0][1][0]=mu/2.;
  }

  if(dim==3){
    
    double Ex=1000.0; double Ey=2000.0; double Ez=1500.0;
    double nuxy=0.3; double nuyz=0.3; double nuzx=0.3;
    double nuyx=(Ey/Ex)*nuxy; double nuzy=(Ez/Ey)*nuyz; double nuxz=(Ex/Ez)*nuzx; 
    double Delta=(1.0-nuxy*nuyx-nuyz*nuzy-nuzx*nuxz-2.0*nuxy*nuyz*nuzx)/(Ex*Ey*Ez);
    double muxy=500.0; double muyz=600.0; double muzx=700.0;
    ElasticModulus[0][0][0][0]=(1.0-nuyz*nuzy)/(Ey*Ez*Delta);
    ElasticModulus[0][0][1][1]=(nuyx+nuzx*nuyz)/(Ey*Ez*Delta);
    ElasticModulus[0][0][2][2]=(nuzx+nuyx*nuzy)/(Ey*Ez*Delta);
    ElasticModulus[1][1][0][0]=(nuxy+nuxz*nuzy)/(Ez*Ex*Delta);
    ElasticModulus[1][1][1][1]=(1.0-nuzx*nuxz)/(Ez*Ex*Delta);
    ElasticModulus[1][1][2][2]=(nuzy+nuzx*nuxy)/(Ez*Ex*Delta);
    ElasticModulus[2][2][0][0]=(nuxz+nuxy*nuyz)/(Ex*Ey*Delta);
    ElasticModulus[2][2][1][1]=(nuyz+nuxz*nuyx)/(Ex*Ey*Delta);
    ElasticModulus[2][2][2][2]=(1.0-nuxy*nuyx)/(Ez*Ex*Delta);
    ElasticModulus[0][1][0][1]=muxy; ElasticModulus[0][1][1][0]=muxy;
    ElasticModulus[0][2][0][2]=muzx; ElasticModulus[0][2][2][0]=muzx;
    ElasticModulus[1][0][0][1]=muxy; ElasticModulus[1][0][1][0]=muxy;
    ElasticModulus[1][2][1][2]=muyz; ElasticModulus[1][2][2][1]=muyz;
    ElasticModulus[2][0][0][2]=muzx; ElasticModulus[2][0][2][0]=muzx;
    ElasticModulus[2][1][1][2]=muyz; ElasticModulus[2][1][2][1]=muyz;
  }
  computeRotatedModulii<dim>(ElasticModulus, rotationMatrices, A_phi);
  for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)for(unsigned int k=0;k<dim;k++)for(unsigned int l=0;l<dim;l++)ElasticTangentModulus[i][j][k][l]=0.;

  for(unsigned int N=0;N<n_diff_grains;N++){
    for(unsigned int i=0;i<dim;i++){
      for(unsigned int j=0;j<dim;j++){
	for(unsigned int k=0;k<dim;k++){
	  for(unsigned int l=0;l<dim;l++){
	    ElasticTangentModulus[i][j][k][l]+=h_phi[N]*A_phi[N][i][j][k][l]/h_total;
	  }
	}
      }
    }
  }

  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      for(unsigned int k=0;k<dim;k++){
	for(unsigned int l=0;l<dim;l++){
	  secondPiola[i][j]+=ElasticTangentModulus[i][j][k][l]*largeStrain[k][l];
	}
      }
    }
  }
  
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      for(unsigned int k=0;k<dim;k++){
	piolaStress[i][j]+=defGrad[i][k]*secondPiola[k][j];
      }
    }
  }

  // calculate Cauchy stress
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      for(unsigned int k=0;k<dim;k++){
	for(unsigned int l=0;l<dim;l++){
	  cauchyStress[i][j]+=ElasticTangentModulus[i][j][k][l]*smallStrain[k][l];
	}
      }
    }
  }
  
  //mechanical energy calculation
 
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++) {
	  if(!isFiniteStrain){
	  freeEnergyMech+=(0.5)*(smallStrain[i][j]*ElasticTangentModulus[i][j][k][l]*smallStrain[k][l])*fe_values.JxW(q);
	  }
	  else{
	    freeEnergyMech+=(0.5)*(largeStrain[i][j]*ElasticTangentModulus[i][j][k][l]*largeStrain[k][l])*fe_values.JxW(q);
	  }
	  
	}
  
  double phiSum=0.0; 
  for(unsigned int i=0;i<n_diff_grains;i++)phiSum+=phi[i]*phi[i];
  for(unsigned int i=0;i<n_diff_grains;i++){
    dF[i]+=(4./3.)*(-12.0*phi[i]*phi[i]+ 12 *phi[i]*phiSum)*fe_values.JxW(q);
  }
  for(unsigned int N=0;N<n_diff_grains;N++){
    for(unsigned int i=0;i<dim;i++){
      for(unsigned int j=0;j<dim;j++){
	for(unsigned int k=0;k<dim;k++){
	  for(unsigned int l=0;l<dim;l++){
	    if(isFiniteStrain){
	    dE[N]+=0.5*(sin(PI*phi[N])/h_total)*(PI/2.0)*(A_phi[N][i][j][k][k]-ElasticTangentModulus[i][j][k][l]/h_total)*largeStrain[i][j]*largeStrain[k][l]*fe_values.JxW(q);
	    }
	    else{
	      dE[N]+=0.5*(sin(PI*phi[N])/h_total)*(PI/2.0)*(A_phi[N][i][j][k][k]-ElasticTangentModulus[i][j][k][l]/h_total)*smallStrain[i][j]*smallStrain[k][l]*fe_values.JxW(q);
	    }
	  }
	}
      }
    }
  }
  
  for(unsigned int i=0;i<n_diff_grains;i++)dF_dE[i]=dF[i]-dE[i];
  history[q]->stress=secondPiola[0][0];
}


//mechanics residual implementation
template <int dim>
void residualForMechanics(FEValues<dim>& fe_values,FEFaceValues<dim> & fe_face_values,const typename DoFHandler<dim>::active_cell_iterator & cell, unsigned int DOF, Table<1, double >& ULocal, Table<1, double>& ULocalConv, deformationMap<double, dim>& defMap, unsigned int currentIteration,  std::vector<historyVariables<dim>* >& history, Vector<double>& RLocal, FullMatrix<double>& KLocal, double fractionalTime, double & freeEnergyMech, std::vector<double>& dF, std::vector<double>&dE, std::vector<double>& dF_dE ,unsigned int & currentIncrement){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
  //quadrature loop
  Table<1, double> traction(dim);
  traction[0]=100.0; traction[1]=0.0;
  
  //surface boundary condition
  //if(currentIncrement>=3 && currentIncrement<=9){
  /*for(unsigned int faceID=0;faceID<2*dim;faceID++){
    if(cell->face(faceID)->at_boundary()){
      //std::cout<<"control here1";
      if((std::abs(cell->face(faceID)->center()[0]-0.5)<1e-8)){
	//std::cout<<"control here2";
	fe_face_values.reinit(cell,faceID);
	for(unsigned int I=0;I<dofs_per_cell;I++){
	  unsigned int ci=fe_values.get_fe().system_to_component_index(I).first;
	  if(ci>=0 && ci<dim){
	    for(unsigned int q=0;q<n_face_q_points;q++){
	      //std:: cout<<"control here3"; 
	      double var=0;
	      var=fe_face_values.shape_value(I,q)*traction[ci]*fe_face_values.JxW(q);
	      RLocal[I]+=var;
	      //RLocal[I]+=fe_face_values.shape_value(I,q)*traction[ci]*fe_face_values.JxW(q);		
	      //if(var>1e-6){std::cout<<"ho gya bhai\t"<<var<<" "; exit(-1);}
	    }
	  }
	}
      }
    }
  }*/
    //}*/

  for (unsigned int q=0; q<n_q_points; q++){
    Table<1, double>phi(n_diff_grains), grainAngle(n_diff_grains);
    FullMatrix<double> Rotation(dim, dim);
    Table<1, double> h_phi(n_diff_grains); std::vector<Table<4, double> >A_phi;
    //evaluate stress
    Table<2, double> stress(dim, dim); Table<2, double> largeStrain(dim, dim);
    Table<4, double> ElasticModulus(dim, dim, dim, dim),ElasticTangentModulus(dim, dim, dim, dim);
    for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)for(unsigned int k=0;k<dim;k++)for(unsigned int l=0;l<dim;l++)ElasticModulus[i][j][k][l]=0.;
    Table<2, double>defGrad(dim,  dim); Table<2, double>secondPiola(dim, dim);
    Table<2, double>piolaStress(dim, dim); Table<2, double>smallStrain(dim, dim);
    Table<2, double>cauchyStress(dim, dim);
    evaluateStress<double, dim>(q,fe_values, DOF, ULocal, ULocalConv, defMap, currentIteration, history, grainAngle, phi,secondPiola, piolaStress, cauchyStress ,ElasticModulus,ElasticTangentModulus,largeStrain,smallStrain,Rotation,defGrad, fractionalTime, freeEnergyMech, dF, dE, dF_dE,currentIncrement,h_phi, A_phi);
    double h_total=0.;
    for(unsigned int N=0;N<n_diff_grains;N++)h_total+=h_phi[N];
   
 
    for(unsigned int I=0;I<dofs_per_cell;I++){
      const unsigned int ci = fe_values.get_fe().system_to_component_index(I).first - DOF;
      if (ci>=0 && ci<dim){

	for (unsigned int J = 0; J<dim; J++){
	  if(isFiniteStrain){
	    RLocal(I) += fe_values.shape_grad(I, q)[J]*(piolaStress[ci][J])*fe_values.JxW(q);
	  }
	  else{
	    RLocal[I] += fe_values.shape_grad(I, q)[J]* cauchyStress[ci][J]*fe_values.JxW(q);
	  }
	}
      }
      if(ci>=dim && ci<n_diff_grains+dim){
	Table<4, double>tempModulus(dim, dim, dim, dim);
	for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)for(unsigned int k=0;k<dim;k++)for(unsigned int l=0;l<dim;l++)tempModulus[i][j][k][l]=0.;
	
	for(unsigned int i=0;i<dim;i++){
	  for(unsigned int j=0;j<dim;j++){
	    for(unsigned int k=0;k<dim;k++){
	      for(unsigned int l=0;l<dim;l++){
		tempModulus[i][j][k][l]+=((sin(PI*phi[ci-dim])*PI/2.0)*A_phi[ci-dim][i][j][k][l]/h_total) - ElasticTangentModulus[i][j][k][l]*(PI*sin(PI*phi[ci-dim]))/(2.0*h_total);
	      }
	    }
	  }
	}
	for(unsigned int i=0;i<dim;i++){
	  for(unsigned int j=0;j<dim;j++){
	    for(unsigned int k=0;k<dim;k++){
	      for(unsigned int l=0;l<dim;l++){
		if(isFiniteStrain){
		  RLocal(I)+=Mobility_m*0.5*fe_values.shape_value(I,q)*largeStrain[i][j]*tempModulus[i][j][k][l]*largeStrain[k][l]*fe_values.JxW(q);
		}
		else{
		  RLocal[I]+=Mobility_m*0.5*fe_values.shape_value(I,q)*smallStrain[i][j]*tempModulus[i][j][k][l]*smallStrain[k][l]*fe_values.JxW(q);
		}
	      }
	    }
	  }
	}
	
		 
      }

      
    }

    // jacobian matrix (mechanics residual, material stiffness)
    if(isFiniteStrain){
      for(unsigned int A=0;A<dofs_per_cell;A++){
	unsigned int i=fe_values.get_fe().system_to_component_index(A).first;
	if(i>=0 && i<dim){
	  for(unsigned int J=0;J<dim;J++){
	    for(unsigned int B=0;B<dofs_per_cell;B++){
	      unsigned int l=fe_values.get_fe().system_to_component_index(B).first;
	      if(l>=0 && l<dim){
		for(unsigned int L=0;L<dim;L++){
		  for(unsigned int K=0;K<dim;K++){
		    for(unsigned int M=0;M<dim;M++){
		      KLocal(A,B)+=fe_values.shape_grad(A,q)[J]*fe_values.shape_grad(B,q)[L]*defGrad[i][K]*defGrad[l][M]*ElasticTangentModulus[K][J][L][M]*fe_values.JxW(q);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    if(!isFiniteStrain){
      for(unsigned int A=0;A<dofs_per_cell;A++){
	unsigned int i=fe_values.get_fe().system_to_component_index(A).first;
	if(i>=0 && i<dim){
	  for(unsigned int B=0;B<dofs_per_cell;B++){
	    unsigned int k=fe_values.get_fe().system_to_component_index(B).first;
	    if(k>=0 && k<dim){
	      for(unsigned int j=0;j<dim;j++){
		for(unsigned int l=0;l<dim;l++){
		  KLocal(A,B)+=fe_values.shape_grad(A,q)[j]*ElasticTangentModulus[i][j][k][l]*fe_values.shape_grad(B,q)[l]*fe_values.JxW(q);
		}
	      }
	    }
	  }
	}
      }
      
    }

    // jacobian matrix (geometric stiffness term)
    if(isFiniteStrain){
      for(unsigned int A=0;A<dofs_per_cell;A++){
	unsigned int i=fe_values.get_fe().system_to_component_index(A).first;
	if(i>=0 && i<dim){
	  for(unsigned int J=0;J<dim;J++){
	    for(unsigned int B=0;B<dofs_per_cell;B++){
	      unsigned int b=fe_values.get_fe().system_to_component_index(B).first;
	      if(b>=0 && b<dim){
		for(unsigned int K=0;K<dim;K++){
		  KLocal(A,B)+=fe_values.shape_grad(A,q)[K]*secondPiola[K][J]*fe_values.shape_grad(B,q)[J]*fe_values.JxW(q);
		}
	      }
	    }
	  } 
	}
      }
    }
    
    for(unsigned int A=0;A<dofs_per_cell;A++){
      int ca=fe_values.get_fe().system_to_component_index(A).first;
      if(ca>=0 && ca<dim){
	for(unsigned int j=0;j<dim;j++){
	  for(unsigned int B=0;B<dofs_per_cell;B++){
	    int cb=fe_values.get_fe().system_to_component_index(B).first;
	    if(cb>=dim && cb<dim+n_diff_grains){
	      Table<4, double>tempModulus(dim, dim, dim, dim);
	      for(unsigned int I=0;I<dim;I++)for(unsigned int J=0;J<dim;J++)for(unsigned int K=0;K<dim;K++)for(unsigned int L=0;L<dim;L++)tempModulus[I][J][K][L]=0.;
	      for(unsigned int I=0;I<dim;I++){
		for(unsigned int J=0;J<dim;J++){
		  for(unsigned int K=0;K<dim;K++){
		    for(unsigned int L=0;L<dim;L++){
		      tempModulus[I][J][K][L]+=(sin(PI*phi[cb-dim])*PI/2.0)*A_phi[cb-dim][I][J][K][L]/h_total-ElasticTangentModulus[I][J][K][L]*(sin(PI*phi[cb-dim])*PI/2.0)/h_total;
		    }
		  }
		}
	      }
	      if(!isFiniteStrain){
		for(unsigned int k=0;k<dim;k++){
		  for(unsigned int l=0;l<dim;l++){
		    KLocal(A,B)+=fe_values.shape_value(B,q)*tempModulus[ca][j][k][l]*smallStrain[k][l]*fe_values.shape_grad(A,q)[j]*fe_values.JxW(q);
		  }
		}
	      }
	      if(isFiniteStrain){
		for(unsigned int I=0;I<dim;I++){
		  for(unsigned int K=0;K<dim;K++){
		    for(unsigned int L=0;L<dim;L++){
		      KLocal(A,B)+=fe_values.shape_value(B,q)*defGrad[ca][I]*tempModulus[I][j][K][L]*largeStrain[K][L]*fe_values.shape_grad(A,q)[j]*fe_values.JxW(q);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    
    for(unsigned int A=0;A<dofs_per_cell;A++){
      int ca=fe_values.get_fe().system_to_component_index(A).first;
      if(ca>=dim && ca<dim+n_diff_grains){
	Table<4, double>tempModulus(dim, dim, dim, dim);
	for(unsigned int I=0;I<dim;I++)for(unsigned int J=0;J<dim;J++)for(unsigned int K=0;K<dim;K++)for(unsigned int L=0;L<dim;L++)tempModulus[I][J][K][L]=0.;
	for(unsigned int I=0;I<dim;I++){
	  for(unsigned int J=0;J<dim;J++){
	    for(unsigned int K=0;K<dim;K++){
	      for(unsigned int L=0;L<dim;L++){
		tempModulus[I][J][K][L]+=(sin(PI*phi[ca-dim])*PI/2.0)*A_phi[ca-dim][I][J][K][L]/h_total-ElasticTangentModulus[I][J][K][L]*(sin(PI*phi[ca-dim])*PI/2.0)/h_total;
	      }
	    }
	  }
	}
	
	if(!isFiniteStrain){
	  for(unsigned int k=0;k<dim;k++){
	    for(unsigned int l=0;l<dim;l++){
	      for(unsigned int B=0;B<dofs_per_cell;B++){
		int cb=fe_values.get_fe().system_to_component_index(B).first;
		if(cb>=0 && cb<dim){
		  for(unsigned int j=0;j<dim;j++){
		    KLocal(A,B)+=Mobility_m*fe_values.shape_grad(B,q)[j]*tempModulus[cb][j][k][l]*smallStrain[k][l]*fe_values.shape_value(A,q)*fe_values.JxW(q);
		  }
		}
	      }
	    }
	  }
	}
	if(isFiniteStrain){
	  for(unsigned int B=0;B<dofs_per_cell;B++){
	    unsigned int cb=fe_values.get_fe().system_to_component_index(B).first;
	    if(cb>=0 && cb<dim){
	      for(unsigned int I=0;I<dim;I++){
		for(unsigned int J=0;J<dim;J++){
		  for(unsigned int K=0;K<dim;K++){
		    for(unsigned int L=0;L<dim;L++){
		      for(unsigned int M=0;M<dim;M++){
			KLocal(A,B)+=Mobility_m* fe_values.shape_value(A,q)*(0.5)*tempModulus[I][J][K][L]*largeStrain[I][J]*( (double)(M==K) *defGrad[cb][L] + defGrad[cb][K]* (double)(M==L) )*fe_values.shape_grad(B,q)[M]*fe_values.JxW(q);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    for(unsigned int A=0;A<dofs_per_cell;A++){
      int ca=fe_values.get_fe().system_to_component_index(A).first;
      if(ca>=dim && ca<n_diff_grains+dim){
	for(unsigned int B=0;B<dofs_per_cell;B++){
	  int cb=fe_values.get_fe().system_to_component_index(B).first;
	  if(cb>=dim && cb<n_diff_grains+dim){
	    
	    Table<4, double> tempModulus(dim, dim, dim, dim);
	    for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)for(unsigned int k=0;k<dim;k++)for(unsigned int l=0;l<dim;l++)tempModulus[i][j][k][l]=0.;
	    
	    if(ca==cb){
	      for(unsigned int I=0;I<dim;I++){
		for(unsigned int J=0;J<dim;J++){
		  for(unsigned int K=0;K<dim;K++){
		    for(unsigned int L=0;L<dim;L++){
		      tempModulus[I][J][K][L]+=(PI*PI*cos(PI*phi[ca-dim])/2.0)*A_phi[ca-dim][I][J][K][L]/h_total - (PI*PI*sin(PI*phi[ca-dim])*sin(PI*phi[ca-dim])/2.0)*A_phi[ca-dim][I][J][K][L]/(h_total*h_total) - ElasticTangentModulus[I][J][K][L]*(PI*PI*cos(PI*phi[ca-dim])/2.0)/h_total + ElasticTangentModulus[I][J][K][L]*
			(PI*PI*sin(PI*phi[ca-dim])*sin(PI*phi[ca-dim])/2.0)/(h_total*h_total);
		    }
		  }
		}
	      }
	    }
	    else{
	      for(unsigned int I=0;I<dim;I++){
		for(unsigned int J=0;J<dim;J++){
		  for(unsigned int K=0;K<dim;K++){
		    for(unsigned int L=0;L<dim;L++){
		      tempModulus[I][J][K][L]+=-(PI*PI*sin(PI*phi[ca-dim])*sin(PI*phi[cb-dim])/4.0)*(A_phi[ca-dim][I][J][K][L]+A_phi[cb-dim][I][J][K][L])/(h_total*h_total) + (PI*PI*sin(PI*phi[ca-dim])*sin(PI*phi[cb-dim])/2.0)*ElasticTangentModulus[I][J][K][L]/(h_total*h_total); 
		    }
		  }
		}
	      }
	    }
	      
	    for(unsigned int i=0;i<dim;i++){
	      for(unsigned int j=0;j<dim;j++){
		for(unsigned int k=0;k<dim;k++){
		  for(unsigned int l=0;l<dim;l++){
		    if(!isFiniteStrain){
		      KLocal(A,B)+=Mobility_m*0.5*fe_values.shape_value(A,q)*smallStrain[i][j]*tempModulus[i][j][k][l]*smallStrain[k][l]*fe_values.shape_value(B,q)*fe_values.JxW(q);
		    }
		    if(isFiniteStrain){
		      KLocal(A,B)+=Mobility_m*0.5*fe_values.shape_value(A,q)*largeStrain[i][j]*tempModulus[i][j][k][l]*largeStrain[k][l]*fe_values.shape_value(B,q)*fe_values.JxW(q);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    
    
  }
}


#endif /* MECHANICS_H_ */
