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
  alpha(0.0), beta (dim+1), betaIteration(dim+1), Ep(dim+1), EpIteration(dim+1){}

  //using std:::map to store time history variables
  double alpha, alphaIteration;
  double eqvStress, equvStrain;
  double elasStrain11, elasStrain12;
  dealii::Table<1, double > beta, betaIteration;
  double stress;
  dealii::Table<1, double > Ep, EpIteration;
};

template<int dim>
void assign_grain_id(Table<1, double>& angle){
  std::srand(5.55);
  for(unsigned int i=0;i<n_diff_grains;i++){
    angle[i]=(double)(std::rand()%180);
    // std::cout<<angle[i]<<" ";
  }
  //exit(-1);
}

template<int dim>
void crystalRotation(Table<1, double>& phi, Table<1, double>& angle, FullMatrix<double>& crystalRotation){
  double param=0.;
  FullMatrix<double> rotation(dim, dim);//(dim+1, dim+1);
  for(unsigned int I=0;I<n_diff_grains;I++){
    param=angle[I]*PI/180.0;
    // rotation=IdentityMatrix(dim, dim);
    rotation(0,0)=cos(param); rotation(0,1)=-sin(param); rotation(1,0)=sin(param);  rotation(1,1)=cos(param);
    //for(int i=0;i<dim;i++)for(int j=0;j<dim;j++)std::cout<<rotation(i,j)<<" ";std::cout<<"\n";
    for(unsigned int i=0;i<dim;i++){
      for(unsigned int j=0;j<dim;j++){
	crystalRotation(i,j)+=phi[I]*rotation(i,j);
      }
    }
  }
  //exit(-1);
}

//mechanics implementation
template <class T, int dim>
  void evaluateStress(unsigned int q, FEValues<dim>& fe_values, const unsigned int DOF, const Table<1, T>& ULocal, const deformationMap<T, dim>& defMap, unsigned int currentIteration,  std::vector<historyVariables<dim>*>& history, Table<2, double>& stress, FullMatrix<double>& C_consistent, Table<4, double>&ElasticModulus, Table<2, double>&defGrad, double fractionalTime, double &freeEnergyMech){ 


  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  //update history variables at the start of increment
  if (currentIteration==0){
    history[q]->alpha=history[q]->alphaIteration;
    history[q]->beta=history[q]->betaIteration;
    history[q]->Ep=history[q]->EpIteration;
  }
  Table<1,double> grainOrientationAngle(n_diff_grains), phi(n_diff_grains);
  //calculate OrderParameter at quadrature point
  for(unsigned int i=0;i<n_diff_grains;i++){grainOrientationAngle[i]=0.; phi[i]=0.;}
  for(unsigned int i=0;i<dofs_per_cell;i++){
    int ci=fe_values.get_fe().system_to_component_index(i).first-dim;
    if(ci>=0 && ci<n_diff_grains){
      phi[ci]+=fe_values.shape_value(i,q)*ULocal[i];
    }
  }
  assign_grain_id<dim>(grainOrientationAngle);
  FullMatrix<double>Rotation(dim, dim);
  crystalRotation<dim>(phi, grainOrientationAngle, Rotation);
  
  //material properties
  Table<2,double> Fe (dim, dim);
  Table<2,double> E (dim, dim);
  Table<1,double>vec_E(dim+1),vec_Ep(dim+1);

  
  FullMatrix<double> C(dim+1,dim+1), C_inv(dim+1,dim+1), temp_C(dim+1,dim+1), C_algorithmic(dim+1,dim+1),C1(dim+1,dim+1);
  Table<2,double> projection(dim+1,dim+1), Q(dim+1,dim+1);
  Table<2,double> gradU (dim, dim);
  double Y=elasticModulus, nu=PoissonsRatio;
  //Lame parameters
  double lambda=(nu*Y)/((1+nu)*(1-2*nu)), mu=Y/(2*(1+nu));
  //double kappa= lambda+(2.0/3.0)*mu;
  
  //double H=0.0;
  
  //double D_g=0.0;
  //Fe
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      // Fe[i][j]=defMap.F[q][i][j];
      gradU[i][j]=defMap.F[q][i][j] - (double)(i==j);
    }
  }
  
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      E[i][j] = 0.5*(gradU[i][j]+gradU[j][i]);// strain      
    }
  }
  //strain due to rotation
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      E[i][j]=E[i][j]-fractionalTime*lambda1*Rotation(i,j);
    }
  }
  for(unsigned int i=0;i<dim+1;i++){
    for(unsigned int j=0;j<dim+1;j++){
      C(i,j)=0.0;
      temp_C(i,j)=0.0;
      C_algorithmic(i,j)=0;
      projection[i][j]=0.0;
      Q[i][j]=0.0;  
    }
  }
  
  Table<1,double> basis1(dim),basis2(dim), rotbasis1(dim), rotbasis2(dim);
  Table<2, double> matBasis1(dim, dim), matBasis2(dim, dim);
  for(unsigned int i=0;i<dim;i++) for(unsigned int j=0;j<dim;j++){matBasis1[i][j]=0.; matBasis2[i][j]=0.;}
  basis1[0]=1.; basis1[1]=0.;  basis2[0]=0.; basis2[1]=1.;  rotbasis1[0]=0.; rotbasis1[1]=0.;  rotbasis2[0]=0.; rotbasis2[1]=0.;
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      rotbasis1[i]+=Rotation(i,j)*basis1[j];
      rotbasis2[i]+=Rotation(i,j)*basis2[j];
    }
  }
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      matBasis1[i][j]+=rotbasis1[i]*rotbasis1[j];
      matBasis2[i][j]+=rotbasis2[i]*rotbasis2[j];
    }
  }
  
  
  projection[0][0]=projection[1][1]=2./3.; projection[0][1]=projection[1][0]=-1./3.;  projection[2][2]=2;
  C(0,0)=C(1,1)=(Y/(1-nu*nu));   C(0,1)=C(1,0)=nu*Y/(1-nu*nu);  C(2,2)=Y/(2*(1+nu));
  //Q[0][0]=Q[1][1]=Q[0][1]=std::sqrt(1./2.);   Q[1][0]=(-1)*std::sqrt(1./2.);     Q[2][2]=1;
  //C_inv.invert(C);
  /* ElasticModulus[0][0][0][0]=Y/(1-nu*nu);  ElasticModulus[1][1][1][1]=Y/(1-nu*nu);
  ElasticModulus[0][0][1][1]=nu*Y/(1-nu*nu);   ElasticModulus[1][1][0][0]=Y*nu/(1-nu*nu);
  ElasticModulus[0][1][0][1]= mu/2.;   ElasticModulus[0][1][1][0]=mu/2.;   ElasticModulus[1][0][0][1]=mu/2.;   ElasticModulus[1][0][1][0]=mu/2.;*/
  double nu12=0.32; double nu21=0.45;
  ElasticModulus[0][0][0][0]=alpha1/(1-nu12*nu21);  ElasticModulus[0][0][1][1]=alpha1*nu21/(1-nu12*nu21);
  ElasticModulus[1][1][0][0]=beta1*nu12/(1-nu12*nu21);  ElasticModulus[1][1][1][1]=beta1/(1-nu12*nu21);
  ElasticModulus[0][1][0][1]= mu/2.;   ElasticModulus[0][1][1][0]=mu/2.;   ElasticModulus[1][0][0][1]=mu/2.;   ElasticModulus[1][0][1][0]=mu/2.;

  //calculate mechanical energy

  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++) {
	  freeEnergyMech+=(0.5)*(E[i][j]*ElasticModulus[i][j][k][l]*E[k][l])*fe_values.JxW(q);
	}
  
  for(int i=0;i<dim;i++) for(int j=0;j<dim;j++) stress[i][j]=0.;

  for(unsigned int i=0;i<dim;i++) for(unsigned int j=0;j<dim;j++) for(unsigned int k=0;k<dim;k++) for(unsigned int l=0;l<dim;l++){stress[i][j]+=ElasticModulus[i][j][k][l]*E[k][l];}
  
  //stress[0][0]=trial_stress[0];stress[1][1]=trial_stress[1];stress[0][1]=trial_stress[2];stress[1][0]=trial_stress[2];
  /*history[q]->EpIteration=vec_Ep;
  history[q]->betaIteration=beta;
  history[q]->alphaIteration=alpha;
  history[q]->elasStrain11=stress[0][0];
  history[q]->elasStrain12=gamma;
  history[q]->stress=stress[0][0];*/
}

template<int dim>
void calculate_diff_rotation(int I, FullMatrix<double>& diffRotation,double fractionalTime){
  Table<1, double>angle(n_diff_grains);
  assign_grain_id<dim>(angle);
  double param=angle[I]*PI/180.0;
  diffRotation(0,0)=cos(param);   diffRotation(0,0)=-sin(param);   diffRotation(0,0)=sin(param);   diffRotation(0,0)=cos(param); 
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      diffRotation(i,j)=-fractionalTime*lambda1*diffRotation(i,j);
    }
  }
}

//mechanics residual implementation
template <int dim>
void residualForMechanics(FEValues<dim>& fe_values, unsigned int DOF, Table<1, double >& ULocal, Table<1, double>& ULocalConv, deformationMap<double, dim>& defMap, unsigned int currentIteration,  std::vector<historyVariables<dim>* >& history, Vector<double>& RLocal, FullMatrix<double>& KLocal, double fractionalTime, double & freeEnergyMech){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //quadrature loop
  for (unsigned int q=0; q<n_q_points; q++){
    
    //evaluate stress
    Table<2, double> stress(dim, dim);
    FullMatrix<double> C(dim+1,dim+1);
    Table<2,double> B(dofs_per_cell,dim+1);
    Table<1,double> _B(dofs_per_cell);
    Table<4, double> ElasticModulus(dim, dim, dim, dim);
    for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)for(unsigned int k=0;k<dim;k++)for(unsigned int l=0;l<dim;l++)ElasticModulus[i][j][k][l]=0.;
    Table<2, double>defGrad(dim, dim);
    evaluateStress<double, dim>(q,fe_values, DOF, ULocal, defMap, currentIteration, history, stress, C, ElasticModulus,defGrad, fractionalTime,freeEnergyMech);
 
  //Lame parameters
 
    for(unsigned int i=0;i<dofs_per_cell;i++){
      const unsigned int ci = fe_values.get_fe().system_to_component_index(i).first - DOF;
      if (ci>=0 && ci<dim){

	for (unsigned int di = 0; di<dim; di++){
	  RLocal(i) += fe_values.shape_grad(i, q)[di]*(stress[ci][di])*fe_values.JxW(q);
	   
	}
      }
      if(ci>=dim && ci<n_diff_grains+dim){
	FullMatrix<double> diffRotation(dim, dim);
	calculate_diff_rotation<dim>(ci-dim, diffRotation,fractionalTime);
	for(unsigned int i=0;i<dim;i++){
	  for(unsigned int j=0;j<dim;j++){
	    RLocal(i)+=fe_values.shape_value(i,q)*stress[i][j]*diffRotation(i,j)*fe_values.JxW(q);
	  }
	}
      }

      
    }
    
   
    
    for(unsigned int i=0;i<dofs_per_cell;i++){
      const unsigned int ci = fe_values.get_fe().system_to_component_index(i).first - DOF;
      if(ci<dim){
	for(unsigned int j=0;j<dim;j++){
	  for(unsigned int k=0;k<dofs_per_cell;k++){
	    const unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
	    if(ck<dim){
	      for(unsigned int l=0;l<dim;l++){
		KLocal(i,k)+=fe_values.shape_grad(i,q)[j]*ElasticModulus[ci][j][ck][l]*fe_values.shape_grad(k,q)[l]*fe_values.JxW(q);
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
	      FullMatrix<double>diffRotation(dim, dim);
	      calculate_diff_rotation<dim>(cb-dim,diffRotation,fractionalTime);
	      for(unsigned int k=0;k<dim;k++){
		for(unsigned int l=0;l<dim;l++){
		  KLocal(A,B)+=fe_values.shape_value(B,q)*ElasticModulus[ca][j][k][l]*diffRotation(k,l)*fe_values.shape_grad(A,q)[j]*fe_values.JxW(q);
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
	calculate_diff_rotation<dim>(ca-dim,diffRotation,fractionalTime);
	for(unsigned int k=0;k<dim;k++){
	  for(unsigned int l=0;l<dim;l++){
	    for(unsigned int B=0;B<dofs_per_cell;B++){
	      int cb=fe_values.get_fe().system_to_component_index(B).first;
	      if(cb>=0 && cb<dim){
		for(unsigned int j=0;j<dim;j++){
		  KLocal(A,B)+=fe_values.shape_grad(B,q)[j]*ElasticModulus[cb][j][k][l]*diffRotation(k,l)*fe_values.shape_value(A,q)*fe_values.JxW(q);
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
	calculate_diff_rotation<dim>(ca-dim, diffRotation1,fractionalTime);
	for(unsigned int B=0;B<dofs_per_cell;B++){
	  int cb=fe_values.get_fe().system_to_component_index(B).first;
	  if(cb>=dim && cb<n_diff_grains+dim){
	    FullMatrix<double>diffRotation2(dim, dim);
	    calculate_diff_rotation<dim>(cb-dim, diffRotation2,fractionalTime);
	    for(unsigned int i=0;i<dim;i++){
	      for(unsigned int j=0;j<dim;j++){
		for(unsigned int k=0;k<dim;k++){
		  for(unsigned int l=0;l<dim;l++){
		    KLocal(A,B)+=fe_values.shape_value(A,q)*diffRotation1(i,j)*ElasticModulus[i][j][k][l]*diffRotation2(k,l)*fe_values.shape_value(B,q)*fe_values.JxW(q);
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
