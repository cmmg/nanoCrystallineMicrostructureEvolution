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
  alpha(0.0), beta (dim+1), betaIteration(dim+1), Ep(dim+1), EpIteration(dim+1), stress(0.0), orientation(0.0) ,elasStrain12(0.0){}

  //using std:::map to store time history variables
  double alpha, alphaIteration;
  double eqvStress, equvStrain;
  double elasStrain11, elasStrain12;
  dealii::Table<1, double > beta, betaIteration;
  double stress, orientation;
  dealii::Table<1, double > Ep, EpIteration;
};

template<int dim>
void assign_grain_id(Table<1, double>& angle,unsigned int & currentIncrement){
  std::srand(5.55);
  for(unsigned int i=0;i<n_diff_grains;i++){
    if(i<4)angle[i]=90.0;
    else angle[i]=0.0;
  }

  //for(unsigned int i=0;i<n_diff_grains;i++){
      //angle[i]=(double)(std::rand()%180);
  //  if(i%3==0) angle[i]=0.0;
  //  if(i%3==1) angle[i]=45.0;
  //  if(i%3==2) angle[i]=90.0;
  // }
    //angle[0]=0.0; angle[1]=30.0; angle[2]=60.0; angle[3]=90.0; angle[4]=0.0; angle[5]=30.0; angle[6]=60.0; angle[7]=90.0;
  //if(currentIncrement<10)angle[0]=0.0;
  //else{
  //angle[0]=(((double)(currentIncrement)-10)/100.0)*90.0; 
  //}
  
  //  angle[1]=90.0;
  //exit(-1);
}

template<int dim>
void crystalRotation(Table<1, double>& phi, Table<1, double>& angle, std::vector<FullMatrix<double> >& rotationMatrices){
  double param=0.;
  FullMatrix<double> rotation(dim, dim);//(dim+1, dim+1);
  for(unsigned int I=0;I<n_diff_grains;I++){
    param=angle[I]*PI/180.0;
    rotation(0,0)=cos(param); rotation(0,1)=-sin(param); rotation(1,0)=sin(param);  rotation(1,1)=cos(param);
    rotationMatrices.push_back(rotation);
    
  }
  //exit(-1);
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
  void evaluateStress(unsigned int q, FEValues<dim>& fe_values, const unsigned int DOF, const Table<1, T>& ULocal,const Table<1, T>& ULocalConv, const deformationMap<T, dim>& defMap, unsigned int currentIteration,  std::vector<historyVariables<dim>*>& history,Table<1, double>&grainAngle,Table<1, double>&phi, Table<2, double>& stress, FullMatrix<double>& C_consistent,Table<4, double>&ElasticModulus ,Table<4, double>&ElasticTangentModulus, Table<2, double>& strain,FullMatrix<double>&Rotation,Table<2, double>&defGrad, double fractionalTime, double & freeEnergyMech, unsigned int & currentIncrement,Table<1, double>&h_phi, std::vector<Table<4, double> >&A_phi){ 


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
  history[q]->orientation=0.0;
  for(unsigned int i=0;i<n_diff_grains;i++){
    history[q]->orientation+=phi[i]*grainAngle[i];
  }
  //material properties
  Table<2,double> Fe (dim, dim);
  Table<2,double> E (dim, dim);
  Table<2,double> gradU (dim, dim);
  double Y=elasticModulus, nu=PoissonsRatio;
  //Lame parameters
  double lambda=(nu*Y)/((1+nu)*(1-2*nu)), mu=Y/(2*(1+nu));

  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      // Fe[i][j]=defMap.F[q][i][j];
      gradU[i][j]=defMap.F[q][i][j] - (double)(i==j);
    }
  }
  
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      E[i][j] = 0.5*(gradU[i][j]+gradU[j][i]);// strain      
      strain[i][j]=E[i][j];
    }
  }
 
  double nu12=0.253; double nu21=0.253;
 
  ElasticModulus[0][0][0][0]=alpha1/(1-nu12*nu21);  ElasticModulus[0][0][1][1]=alpha1*nu21/(1-nu12*nu21);
  ElasticModulus[1][1][0][0]=beta1*nu12/(1-nu12*nu21);  ElasticModulus[1][1][1][1]=beta1/(1-nu12*nu21);
  ElasticModulus[0][1][0][1]= mu/2.;   ElasticModulus[0][1][1][0]=mu/2.;   ElasticModulus[1][0][0][1]=mu/2.;   ElasticModulus[1][0][1][0]=mu/2.;

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

  
  
  //mechanical energy calculation
  //for(unsigned int I=0;I<n_diff_grains;I++){
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++) {
	  freeEnergyMech+=(0.5)*(E[i][j]*ElasticTangentModulus[i][j][k][l]*E[k][l])*fe_values.JxW(q);
	  
	}
  
  
  for(int i=0;i<dim;i++) for(int j=0;j<dim;j++) stress[i][j]=0.;

  for(unsigned int i=0;i<dim;i++) for(unsigned int j=0;j<dim;j++) for(unsigned int k=0;k<dim;k++) for(unsigned int l=0;l<dim;l++){stress[i][j]+=ElasticTangentModulus[i][j][k][l]*E[k][l];}
  
  //stress[0][0]=trial_stress[0];stress[1][1]=trial_stress[1];stress[0][1]=trial_stress[2];stress[1][0]=trial_stress[2];
  /*history[q]->EpIteration=vec_Ep;
  history[q]->betaIteration=beta;
  history[q]->alphaIteration=alpha;
  history[q]->elasStrain11=stress[0][0];
  history[q]->elasStrain12=gamma;*/
  history[q]->stress=stress[0][0];
  
}

template<int dim>
void calculate_diff_rotation(int I, FullMatrix<double>& diffRotation,double fractionalTime, unsigned int & currentIncrement){
  Table<1, double>angle(n_diff_grains);
  assign_grain_id<dim>(angle, currentIncrement);
  double param=angle[I]*PI/180.0;
  diffRotation(0,0)=cos(param);   diffRotation(0,1)=-sin(param);   diffRotation(1,0)=sin(param);   diffRotation(1,1)=cos(param); 

}

//mechanics residual implementation
template <int dim>
void residualForMechanics(FEValues<dim>& fe_values,FEFaceValues<dim> & fe_face_values,const typename DoFHandler<dim>::active_cell_iterator & cell, unsigned int DOF, Table<1, double >& ULocal, Table<1, double>& ULocalConv, deformationMap<double, dim>& defMap, unsigned int currentIteration,  std::vector<historyVariables<dim>* >& history, Vector<double>& RLocal, FullMatrix<double>& KLocal, double fractionalTime, double & freeEnergyMech, unsigned int & currentIncrement){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
  //quadrature loop
  Table<1, double> traction(dim);
  traction[0]=100.0; traction[1]=0.0;
  
  
	      

  for (unsigned int q=0; q<n_q_points; q++){
    Table<1, double>phi(n_diff_grains), grainAngle(n_diff_grains);
    FullMatrix<double> Rotation(dim, dim);
    Table<1, double> h_phi(n_diff_grains); std::vector<Table<4, double> >A_phi;
    //evaluate stress
    Table<2, double> stress(dim, dim);
    Table<2, double> strain(dim, dim);
    FullMatrix<double> C(dim+1,dim+1);
    Table<2,double> B(dofs_per_cell,dim+1);
    Table<1,double> _B(dofs_per_cell);
    Table<4, double> ElasticModulus(dim, dim, dim, dim),ElasticTangentModulus(dim, dim, dim, dim);
    for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)for(unsigned int k=0;k<dim;k++)for(unsigned int l=0;l<dim;l++)ElasticModulus[i][j][k][l]=0.;
    Table<2, double>defGrad(dim, dim);
    
    evaluateStress<double, dim>(q,fe_values, DOF, ULocal, ULocalConv, defMap, currentIteration, history, grainAngle, phi,stress, C, ElasticModulus,ElasticTangentModulus,strain,Rotation,defGrad, fractionalTime, freeEnergyMech, currentIncrement,h_phi, A_phi);
    double h_total=0.;
    for(unsigned int N=0;N<n_diff_grains;N++)h_total+=h_phi[N];
   
  //Lame parameters
 
    for(unsigned int I=0;I<dofs_per_cell;I++){
      const unsigned int ci = fe_values.get_fe().system_to_component_index(I).first - DOF;
      if (ci>=0 && ci<dim){

	for (unsigned int di = 0; di<dim; di++){
	  RLocal(I) += fe_values.shape_grad(I, q)[di]*(stress[ci][di])*fe_values.JxW(q);
	   
	}
      }
      if(ci>=dim && ci<n_diff_grains+dim){
	FullMatrix<double> diffRotation(dim, dim);
	calculate_diff_rotation<dim>(ci-dim, diffRotation,fractionalTime, currentIncrement);

	Table<4, double>tempModulus(dim, dim, dim, dim);// incorporating all differential terms// \partial C / \partial \phi_I
	for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)for(unsigned int k=0;k<dim;k++)for(unsigned int l=0;l<dim;l++)tempModulus[i][j][k][l]=0.;
	// put \frac{\partial C}{\partial \phi} here
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
		RLocal(I)+=Mobility*0.5*fe_values.shape_value(I,q)*strain[i][j]*tempModulus[i][j][k][l]*strain[k][l]*fe_values.JxW(q);
	      }
	    }
	  }
	}
	
	 // if ci<n_diff_grains condition end here	 
      }

      
    }

    
    for (unsigned int i=0;i<dofs_per_cell;i++){
      const unsigned int ci = fe_values.get_fe().system_to_component_index(i).first - DOF;
      if(ci<dim){
	for(unsigned int j=0;j<dim;j++){
	  for(unsigned int k=0;k<dofs_per_cell;k++){
	    const unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
	    if(ck<dim){
	      for(unsigned int l=0;l<dim;l++){
		KLocal(i,k)+=fe_values.shape_grad(i,q)[j]*ElasticTangentModulus[ci][j][ck][l]*fe_values.shape_grad(k,q)[l]*fe_values.JxW(q);
	      }
	    }
	  }
	}
      }
    }
    
    
    // if(currentIncrement>=8){
    for(unsigned int A=0;A<dofs_per_cell;A++){
      int ca=fe_values.get_fe().system_to_component_index(A).first;
      if(ca>=0 && ca<dim){
	for(unsigned int j=0;j<dim;j++){
	  for(unsigned int B=0;B<dofs_per_cell;B++){
	    int cb=fe_values.get_fe().system_to_component_index(B).first;
	    if(cb>=dim && cb<dim+n_diff_grains){
	      FullMatrix<double>diffRotation(dim, dim);
	      calculate_diff_rotation<dim>(cb-dim,diffRotation,fractionalTime, currentIncrement);
	      Table<4, double>tempModulus(dim, dim, dim, dim);// incorporating all differential terms// \partial C / \partial \phi_I
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
	      for(unsigned int k=0;k<dim;k++){
		for(unsigned int l=0;l<dim;l++){
		  KLocal(A,B)+=fe_values.shape_value(B,q)*tempModulus[ca][j][k][l]*strain[k][l]*fe_values.shape_grad(A,q)[j]*fe_values.JxW(q);
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
	FullMatrix<double>diffRotation(dim, dim);
	calculate_diff_rotation<dim>(ca-dim,diffRotation,fractionalTime, currentIncrement);
	Table<4, double>tempModulus(dim, dim, dim, dim);// incorporating all differential terms// \partial C / \partial \phi_I
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
	
	
	for(unsigned int k=0;k<dim;k++){
	  for(unsigned int l=0;l<dim;l++){
	    for(unsigned int B=0;B<dofs_per_cell;B++){
	      int cb=fe_values.get_fe().system_to_component_index(B).first;
	      if(cb>=0 && cb<dim){
		for(unsigned int j=0;j<dim;j++){
		  KLocal(A,B)+=Mobility*fe_values.shape_grad(B,q)[j]*tempModulus[cb][j][k][l]*strain[k][l]*fe_values.shape_value(A,q)*fe_values.JxW(q);
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
	FullMatrix<double>diffRotation1(dim, dim);
	calculate_diff_rotation<dim>(ca-dim, diffRotation1,fractionalTime, currentIncrement);
	for(unsigned int B=0;B<dofs_per_cell;B++){
	  int cb=fe_values.get_fe().system_to_component_index(B).first;
	  if(cb>=dim && cb<n_diff_grains+dim){
	    FullMatrix<double>diffRotation2(dim, dim);
	    calculate_diff_rotation<dim>(cb-dim, diffRotation2,fractionalTime, currentIncrement);
	    Table<4, double> tempModulus(dim, dim, dim, dim);
	    for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)for(unsigned int k=0;k<dim;k++)for(unsigned int l=0;l<dim;l++)tempModulus[i][j][k][l]=0.;
	    
	    if(ca==cb){
	      for(unsigned int I=0;I<dim;I++){
		for(unsigned int J=0;J<dim;J++){
		  for(unsigned int K=0;K<dim;K++){
		    for(unsigned int L=0;L<dim;L++){
		      tempModulus[I][J][K][L]+=(PI*PI*cos(PI*phi[ca-dim])/2.0)*A_phi[ca-dim][I][J][K][L]/h_total - (PI*PI*sin(PI*phi[ca-dim])*sin(PI*phi[ca-dim])/2.0)*A_phi[ca-dim][I][J][K][L]/(h_total*h_total) - ElasticTangentModulus[I][J][K][L]*(PI*PI*cos(PI*phi[ca-dim])/2.0)/h_total + ElasticTangentModulus[I][J][K][L]*(PI*PI*sin(PI*phi[ca-dim])*sin(PI*phi[ca-dim])/2.0)/(h_total*h_total);
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
		    KLocal(A,B)+=Mobility*0.5*fe_values.shape_value(A,q)*strain[i][j]*tempModulus[i][j][k][l]*strain[k][l]*fe_values.shape_value(B,q)*fe_values.JxW(q);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    
    //}
    
    
  }
}


#endif /* MECHANICS_H_ */
