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


template <int dim>
struct historyVariables{
public:
  historyVariables<dim>():
  
  alpha(0.0),  bLeftcauchy(dim, dim), bLeftcauchyIteration(dim, dim), backstress(dim), backstressIteration(dim),defgrad(dim, dim), defgradIteration(dim, dim){}

  //using std:::map to store time history variables
  double alpha, alphaIteration;
  double eqvstress, equvStrain;
  double elasStrain11, elasStrain22;
  dealii::Table<2, double > bLeftcauchy, bLeftcauchyIteration;
  dealii::FullMatrix<double> defgrad, defgradIteration;
  dealii::Table<1, double > backstress, backstressIteration;
};
//mechanics implementation
template <class T, int dim>
  void evaluateStress(unsigned int q, const unsigned int DOF, const Table<1, T>& ULocal, const deformationMap<T, dim>& defMap, unsigned int currentIteration,unsigned int currentIncrement,  std::vector<historyVariables<dim>*>& history, Table<2, double>& Kirchhoff_stress, Table<2, double>& Piola_stress, Table<4,double>& C,  Table<2,double>& invF, FEValues<dim>& fe_values){ 
  
  //update history variables at the start of increment
  if (currentIteration==0){
    history[q]->alpha=history[q]->alphaIteration;
    history[q]->defgrad=history[q]->defgradIteration;
    history[q]->bLeftcauchy=history[q]->bLeftcauchyIteration;
    history[q]->backstress=history[q]->backstressIteration;
  }
  
  Table<2, double> be(dim,dim), ft(dim, dim), P1(dim, dim), P2(dim, dim), P3(dim, dim), a_elastoplastic(dim, dim),df_dbeta2(dim, dim);
  //Table<2, double> Fn(dim,dim);
  FullMatrix<double> _be(dim,dim),_beT(dim, dim), defgradInv(dim, dim), C_cg(dim, dim), C_cgInv(dim, dim), Ft(dim, dim),FtInv(dim, dim), a_elastic(dim, dim),a_elasticInv(dim, dim);
  FullMatrix<double> h1(dim, dim),h1Inv(dim, dim), h2(dim, dim), h3(dim, dim),h3Inv(dim, dim);
  Table<4, double> II(dim, dim, dim, dim), IIbe(dim, dim, dim, dim),_D(dim, dim, dim, dim), _Dr(dim, dim, dim, dim), _Drr(dim, dim, dim, dim), C_ep(dim, dim, dim, dim), Cr(dim, dim, dim, dim), Crr(dim, dim, dim, dim) ; 
  Table<5, double> C_tr(dim, dim, dim, dim, dim);
  Table<1, double> principal_stress(dim), dev_principal_stress(dim), eigval(dim), _backstress(dim),  _eigval(dim);
  LAPACKFullMatrix<double> b_eTReig(dim, dim), Ident(dim, dim);
  std::vector< Vector<double> > eigvec(dim);
  Table<3, double> nxn(dim, dim, dim);
  
  // double counter=0;
  double nu=PoissonsRatio, Y=elasticModulus;
  double lambda=(nu*Y)/((1+nu)*(1-2*nu)), mu=Y/(2*(1+nu));
  double kappa= lambda+(2.0/3.0)*mu;
  double yieldcriteria=0,norm=0;
  double K=isotropic_hardening;
  double H=kinematic_hardening;
  double tau_y=Yield_stress;
  double gamma=0.0, alpha=history[q]->alpha, function=0.,d_function=0.;
  double tol=1e-6;
  double h4=0.,denominator=0.;
  unsigned int multiplicity;double trace_be=0, det_be=0., det_Ccg=0.;
  double matrix_norm=0.;
  //calculate Ft current deformmation gradient

  _backstress=history[q]->backstress;
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      Ft(i,j)=defMap.F[q][i][j];
      invF(i,j)=defMap.invF[q][i][j];
      Ident(i,j)=(double)(i==j);
    }
  }
  
  defgradInv.invert(history[q]->defgrad);
  FtInv.invert(Ft);
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++){
	ft[i][j]+=Ft(i,k)*defgradInv(k,j);
	C_cg(i,j)+=Ft(k,i)*Ft(k,j);
      }
  
  C_cgInv.invert(C_cg);
  
  det_Ccg=C_cg.determinant();
 
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++){
	  be[i][l]+=ft[i][j]*(history[q]->bLeftcauchy[j][k])*ft[l][k];
	}
  
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      _be(i,j)=be[i][j];
      _beT(i,j)=be[j][i];
      trace_be+=(i==j)*be[i][j];
    }
  }
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      matrix_norm+=(_be(i,j)-_beT(i,j))*(_be(i,j)-_beT(i,j));
    }
  }

  if(matrix_norm>1e-10){std::cout<<"matrix not symmetric";exit(-1);}
  
  b_eTReig=_be;
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++){
	  II[i][j][k][l]=(0.5)*( (i==k)*(j==l) +(i==l)*(j==k));
	  IIbe[i][j][k][l]=0.5*(be[i][k]*be[j][l] + be[i][l]*be[j][k] );
	}
  
  det_be=_be.determinant();
  
  b_eTReig.compute_generalized_eigenvalues_symmetric(Ident, eigvec);
  
  for(unsigned int i=0;i<dim;i++){
    _eigval[i]=abs(b_eTReig.eigenvalue(i));
    
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	nxn[i][j][k]=eigvec[i][j]*eigvec[i][k];
  }
  
  
  if((_eigval[0]<1e-10 && _eigval[0]>-(1e-10)) || (_eigval[1]<1e-10 && _eigval[1]>-1e-10) || (_eigval[2]<1e-10 && _eigval[2]>-1e-10)){std::cout<<"****alert 0 eigenvalue detected";std::cout<<"_eigvalues"<<_eigval[0]<<" "<<_eigval[1]<<" "<<_eigval[2];exit(-1);}

  //eigen values calculated
  
  
  if((_eigval[0]-_eigval[1])<tol && (_eigval[0]-_eigval[1])>-tol)
    {
      if((_eigval[2]-_eigval[1])<tol && (_eigval[2]-_eigval[1])>-tol)
	multiplicity=3;
      else
	multiplicity=2;
    }
  else
    {
      if((_eigval[2]-_eigval[1])<tol && (_eigval[2]-_eigval[1])>-tol)
	multiplicity=2;
      else
	multiplicity=1;
    }
  
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++){
      double s1=(0.5)*(std::log(_eigval[0]*_eigval[1]*_eigval[2]));
      double s2=0,s3=0;
      for(unsigned int a=0;a<dim;a++){
	s2+=(i==a)*(a==j)*(1-0.5*log(_eigval[a]))/_eigval[a];
	s3+=(i==a)*(0.5*log(_eigval[a]))/_eigval[a];
      }
      a_elastic(i,j)=lambda-(lambda*s1*_eigval[j]*(i==j)/_eigval[i])+(2*mu)*s2*_eigval[i]*_eigval[j] +  lambda*s1*(i==j) + (2*mu)* s3* _eigval[i]*(i==j); 
    }
  
  
  //define principal stresses
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
	principal_stress[i]+=(0.5)*(a_elastic(i,j))*std::log(_eigval[j]);
      }
  }
  
  for(unsigned int i=0;i<dim;i++)dev_principal_stress[i]=principal_stress[i]-(1./3.)*(principal_stress[0]+principal_stress[1]+principal_stress[2]);
   
  
  for(unsigned int i=0;i<dim;i++) norm+=(dev_principal_stress[i]-_backstress[i])*(dev_principal_stress[i]-_backstress[i]);
  
  yieldcriteria=std::sqrt(norm)-std::sqrt(2./3.)*(tau_y+K*history[q]->alpha);
 
  if(yieldcriteria<tol)
    {
     
      //elastic
     
      //make _Dr symmetric reference from CS Jog book
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++)
	  for(unsigned int k=0;k<dim;k++)
	    for(unsigned int l=0;l<dim;l++)
	      _Drr[i][j][k][l]=C_cgInv(i,k)*C_cgInv(j,l);

      
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++)
	  for(unsigned int k=0;k<dim;k++)
	    for(unsigned int l=0;l<dim;l++)
	      _Dr[i][j][k][l]+=(1./4.)*(_Drr[i][j][k][l] + _Drr[j][i][k][l] + _Drr[i][j][l][k] + _Drr[j][i][l][k]   );
      
      
      
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++)
	  for(unsigned int k=0;k<dim;k++)
	    for(unsigned int l=0;l<dim;l++)
	      _D[i][j][k][l]+=lambda*C_cgInv(i,j)*C_cgInv(k,l)-lambda*std::log(det_Ccg)*_Dr[i][j][k][l];
      
       if(multiplicity==1){
	
	//calculate P1, P2, P3
	for(unsigned int i=0;i<dim;i++)
	  for(unsigned int j=0;j<dim; j++)
	    for(unsigned int k=0;k<dim;k++){
	      P1[i][j]+=(C_cg(i,k)-(_eigval[2])*(i==k))*(C_cg(k,j)-(_eigval[1])*(k==j))/(((_eigval[0])-(_eigval[1]))*((_eigval[0])-(_eigval[2])));
	      P2[i][j]+=(C_cg(i,k)-(_eigval[0])*(i==k))*(C_cg(k,j)-(_eigval[2])*(k==j))/(((_eigval[1])-(_eigval[2]))*((_eigval[1])-(_eigval[0])));
	      P3[i][j]+=(C_cg(i,k)-(_eigval[1])*(i==k))*(C_cg(k,j)-(_eigval[0])*(k==j))/(((_eigval[2])-(_eigval[0]))*((_eigval[2])-(_eigval[1])));
	    }//for loop end
	
	for(unsigned int i=0;i<dim;i++)
	  for(unsigned int j=0;j<dim;j++)
	    for(unsigned int k=0;k<dim;k++)
	      for(unsigned int l=0;l<dim;l++)
		Crr[i][j][k][l]+=( ((1-log(_eigval[0]))*P1[i][k]*P1[j][l]/(_eigval[0]*_eigval[0])) +  ((1-log(_eigval[1]))*P2[i][k]*P2[j][l]/(_eigval[1]*_eigval[1])) + ((1-log(_eigval[2]))*P3[i][k]*P3[j][l]/(_eigval[2]*_eigval[2])) + (log(_eigval[0])/_eigval[0] - log(_eigval[1])/_eigval[1])*(P1[i][k]*P2[j][l]+P2[i][k]*P1[j][l])/(_eigval[0]-_eigval[1]) +  (log(_eigval[1])/_eigval[1] - log(_eigval[2])/_eigval[2])*(P2[i][k]*P3[j][l]+P3[i][k]*P2[j][l])/(_eigval[1]-_eigval[2]) +  (log(_eigval[2])/_eigval[2] - log(_eigval[0])/_eigval[0])*(P3[i][k]*P1[j][l]+P1[i][k]*P3[j][l])/(_eigval[2]-_eigval[0]));
	
	//make Cr symmetric
	for(unsigned int i=0;i<dim;i++)
	  for(unsigned int j=0;j<dim;j++)
	    for(unsigned int k=0;k<dim;k++)
	      for(unsigned int l=0;l<dim;l++)
	      Cr[i][j][k][l]+=(1./4.)*(Crr[i][j][k][l]+Crr[i][j][l][k]+Crr[j][i][k][l]+Crr[j][i][l][k]);
	
	for(unsigned int i=0;i<dim;i++)
	  for(unsigned int j=0;j<dim;j++)
	    for(unsigned int k=0;k<dim;k++)
	      for(unsigned int l=0;l<dim;l++)
		C_ep[i][j][k][l]=_D[i][j][k][l] + 2*mu*Cr[i][j][k][l];
	for(unsigned int i=0;i<dim;i++)
	  for(unsigned int j=0;j<dim;j++)
	    for(unsigned int k=0;k<dim;k++)
	      for(unsigned int l=0;l<dim;l++){
		if((C_ep[i][j][k][l]-C_ep[i][j][l][k])>tol || (C_ep[i][j][k][l]-C_ep[i][j][l][k])<-tol){std::cout<<"C_not_symmetric_minor_symmetry1_m=1"; exit(-1);}
		if((C_ep[i][j][k][l]-C_ep[j][i][k][l])>tol || (C_ep[i][j][k][l]-C_ep[j][i][k][l])<-tol){std::cout<<"C_not_symmetric_minor_symmetry2_m=1"; exit(-1);}
		if((C_ep[i][j][k][l]-C_ep[j][i][l][k])>tol || (C_ep[i][j][k][l]-C_ep[j][i][l][k])<-tol){std::cout<<"C_not_symmetric_minor_symmetry3_m=1"; exit(-1);}
		if((C_ep[i][j][k][l]-C_ep[k][l][i][j])>tol || (C_ep[i][j][k][l]-C_ep[k][l][i][j])<-tol){std::cout<<"C_not_symmetric_major_symmetry3_m=1"; exit(-1);}	      
	      }
	
	
       }//if multiplicity==1 end
       
       if(multiplicity==3)
	 {
	   for(unsigned int i=0;i<dim;i++)
	     for(unsigned int j=0;j<dim;j++)
	       for(unsigned int k=0;k<dim;k++)
		 for(unsigned int l=0;l<dim;l++){
		   C_ep[i][j][k][l]+=_D[i][j][k][l]+2*mu*(1-std::log(_eigval[0]))*II[i][j][k][l]/(_eigval[0]*_eigval[0]);
		 }
	   for(unsigned int i=0;i<dim;i++)
	     for(unsigned int j=0;j<dim;j++)
	       for(unsigned int k=0;k<dim;k++)
		 for(unsigned int l=0;l<dim;l++){
		   if((C_ep[i][j][k][l]-C_ep[i][j][l][k])>tol || (C_ep[i][j][k][l]-C_ep[i][j][l][k])<-tol){std::cout<<"C_not_symmetric_minor_symmetry1_m=3"; exit(-1);}
		   if((C_ep[i][j][k][l]-C_ep[j][i][k][l])>tol || (C_ep[i][j][k][l]-C_ep[j][i][k][l])<-tol){std::cout<<"C_not_symmetric_minor_symmetry2_m=3"; exit(-1);}
		   if((C_ep[i][j][k][l]-C_ep[j][i][l][k])>tol || (C_ep[i][j][k][l]-C_ep[j][i][l][k])<-tol){std::cout<<"C_not_symmetric_minor_symmetry3_m=3"; exit(-1);}
		   if((C_ep[i][j][k][l]-C_ep[k][l][i][j])>tol || (C_ep[i][j][k][l]-C_ep[k][l][i][j])<-tol){std::cout<<"C_not_symmetric_major_symmetry3_m=3"; exit(-1);}	      
		 }
	   
	   
	   
	 }
       
       //check C for symmetry
       
       for(unsigned int i=0;i<dim;i++)
	 for(unsigned int j=0;j<dim;j++)
	   for(unsigned int k=0;k<dim;k++)
	     for(unsigned int l=0;l<dim;l++)
	       for(unsigned int I=0;I<dim;I++)
		 for(unsigned int J=0;J<dim;J++)
		   for(unsigned int K=0;K<dim;K++)
		     for(unsigned int L=0;L<dim;L++)
		       C[i][j][k][l]+=Ft(i,I)*Ft(j,J)*Ft(k,K)*Ft(l,L)*C_ep[I][J][K][L];
       
       for(unsigned int i=0;i<dim;i++)
	 for(unsigned int j=0;j<dim;j++)
	   for(unsigned int k=0;k<dim;k++)
	     for(unsigned int l=0;l<dim;l++){
	       if((C[i][j][k][l]-C[i][j][l][k])>tol || (C[i][j][k][l]-C[i][j][l][k])<-tol){std::cout<<"C_not_symmetric_minor_symmetry1_m=1"; exit(-1);}
	       if((C[i][j][k][l]-C[j][i][k][l])>tol || (C[i][j][k][l]-C[j][i][k][l])<-tol){std::cout<<"C_not_symmetric_minor_symmetry2_m=1"; exit(-1);}
	       if((C[i][j][k][l]-C[j][i][l][k])>tol || (C[i][j][k][l]-C[j][i][l][k])<-tol){std::cout<<"C_not_symmetric_minor_symmetry3_m=1"; exit(-1);}
	       if((C[i][j][k][l]-C[k][l][i][j])>tol || (C[i][j][k][l]-C[k][l][i][j])<-tol){std::cout<<"C_not_symmetric_major_symmetry3_m=1"; exit(-1);}	      
	     }
       
       
       
       
       
       //to avoid division by zero at first increment C =st.venont kirchhoff model
       if(currentIncrement==1 && currentIteration==0)
	 {
	   for(unsigned int i=0;i<dim;i++)
	     for(unsigned int j=0;j<dim;j++)
	       for(unsigned int k=0;k<dim;k++)
		 for(unsigned int l=0;l<dim;l++)
		   C[i][j][k][l]=lambda*(i==j)*(k==l) +mu*((i==k)*(j==l) + (i==l)*(j==k));
	 }
       
    }
  //elastic part ends
  
  else
    {
      //plastic
      //calculate gamma
      //{std::cout<<"plasticity triggered";exit(-1);}
      unsigned int ctr=0;
      function=yieldcriteria-(2*mu*gamma)-std::sqrt(2./3.)*K*(alpha-history[q]->alpha);
      double mod_f=function;
      if(function<0)mod_f=-1.0*function;
      while(mod_f>tol)
	{
	  d_function=-2*mu-(2./3.)*gamma;
	  gamma=gamma-(function/d_function);
	  alpha=history[q]->alpha + std::sqrt(2./3.)*gamma;
	  function=yieldcriteria-(2*mu*gamma)-std::sqrt(2./3.)*K*(alpha-history[q]->alpha);
	  mod_f=function;
	  if(function<0)mod_f=-1.0*function;
	  ctr++;
	  if(ctr>30){std::cout<<"max itr for gamma";exit(-1);}
	  }

      //  gamma=yieldcriteria/(2*mu);
      //alpha=history[q]->alpha + std::sqrt(2./3.)*gamma;
      //GAMMA AND ALPHA CALCULATED
      
      //update eigen values i.e. principal stretches
      
      for(unsigned int i=0;i<dim;i++){
	principal_stress[i]=principal_stress[i]-2*mu*gamma*(dev_principal_stress[i]-history[q]->backstress[i])/std::sqrt(norm);
	_backstress[i]=(history[q]->backstress[i])+(2./3.)*gamma*H*(dev_principal_stress[i]-history[q]->backstress[i])/std::sqrt(norm);
      } //update eigen values
      norm=0.;
      for(unsigned int i=0;i<dim;i++)dev_principal_stress[i]=principal_stress[i]-(1./3.)*(principal_stress[0]+principal_stress[1]+principal_stress[2]);
      for(unsigned int i=0;i<dim;i++) norm+=(dev_principal_stress[i]-_backstress[i])*(dev_principal_stress[i]-_backstress[i]);
      
      for(unsigned int A=0;A<dim;A++){
	unsigned int B=(A+1)%3;
	unsigned int C=(B+1)%3;
	eigval[A]=_eigval[A]*std::exp( (-2*gamma/(3.*std::sqrt(norm)))* (2*(principal_stress[A]-_backstress[A])-(principal_stress[B]-_backstress[B])-(principal_stress[C]-_backstress[C]))  );
      }


      //calculate C_tr
      for(unsigned int A=0;A<dim;A++)
	for(unsigned int i=0;i<dim;i++)
	  for(unsigned int j=0;j<dim;j++)
	    for(unsigned int k=0;k<dim;k++)
	      for(unsigned int l=0;l<dim;l++){
		unsigned int B=(A+1)%3;
		unsigned int C=(B+1)%3;
		C_tr[A][i][j][k][l]+=(1/((_eigval[A]-_eigval[B])*(_eigval[A]*_eigval[C])))*(IIbe[i][j][k][l]- be[i][j]*be[k][l]- (det_be/eigval[A])*(II[i][j][k][l]-((double)(i==j)-nxn[A][i][j])*((double)(k==l)-nxn[A][k][l])) + eigval[A]*(be[i][j]*nxn[A][k][l] + nxn[A][i][j]*be[k][l]+ (trace_be-4*eigval[A])* nxn[A][i][j]*nxn[A][k][l])   ) ;
	      }
      
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++){
	  double s1=(0.5)*(std::log(eigval[0]*eigval[1]*eigval[2]));
	  double s2=0,s3=0;
	  for(unsigned int a=0;a<dim;a++){
	    s2+=(i==a)*(a==j)*(1-0.5*log(eigval[a]))/eigval[a];
	    s3+=(i==a)*(0.5*log(eigval[a]))/eigval[a];
	  }
	  a_elastic(i,j)=lambda-(lambda*s1*eigval[j]*(i==j)/eigval[i])+(2*mu)*s2*eigval[i]*eigval[j] +  lambda*s1*(i==j) + (2*mu)* s3* eigval[i]*(i==j);
	}
      a_elasticInv.invert(a_elastic);
     
      //calculate h1 h2 h3 to calculate C_ep
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++){
	  df_dbeta2[i][j]+=-(principal_stress[i]-_backstress[i])*(principal_stress[j]-_backstress[j])/(std::pow(norm,3./2.)) + ((double)(i==j)-(1./3.))/(std::sqrt(norm));
	}
      
      //h1
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++){
	  h1Inv(i,j)=(i==j);// + (2./3.)*H*gamma*df_dbeta2[i][j];
	}
      h1.invert(h1Inv);
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++){
	  h2(i,j)+=(double)(i==j) ;
	  for(unsigned int k=0;k<dim;k++){
	    h2(i,j)+=0.;// -(2./3.)*gamma*H*h1(i,k)*df_dbeta2[k][j];
	  }
	}
      
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++){
	  h3Inv(i,j)+=a_elasticInv(i,j);
	  for(unsigned int k=0;k<dim;k++){
	    h3Inv(i,j)+=gamma*df_dbeta2[i][k]*h2(k,j);
	  }
	}
      h3.invert(h3Inv);
      h4=1;//-gamma*(K)*(0); //note that \partial^2{f}/ \partial^2{q}=0
      denominator=0;
      //calculate a_elastoplastic
      
      for(unsigned int m=0;m<dim;m++)
	for(unsigned int n=0;n<dim;n++)
	  for(unsigned int p=0;p<dim;p++)
	    for(unsigned int r=0;r<dim;r++)
	      denominator+=(dev_principal_stress[m]-_backstress[m])*(h2(m,n)*h3(n,p)*h2(p,r)+(2./3.)*H*h1(m,r))*(dev_principal_stress[r]-_backstress[r])/norm;
      
      denominator=h4*denominator+(2./3.)*K;

      Table<1,double> a_epL(dim), a_epR( dim);

      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++)
	  for(unsigned int k=0;k<dim;k++){
	    a_epR[i]+=h3(i,j)*h2(j,k)*(dev_principal_stress[k]-_backstress[k])/sqrt(norm);
	    a_epL[i]+=h3(i,j)*h2(j,k)*(dev_principal_stress[k]-_backstress[k])/sqrt(norm);
	  }
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++)
	  a_elastoplastic[i][j]=h4*a_epL[i]*a_epR[j];
      
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++)
	  a_elastoplastic[i][j]=h3(i,j)-a_elastoplastic[i][j]/denominator;
      //a_elastoplastic calculated
      
      //calculate C_elastplastic algorithmic
      
      for(unsigned int A=0;A<dim;A++)
	for(unsigned int B=0;B<dim;B++)
	  for(unsigned int i=0;i<dim;i++)
	    for(unsigned int j=0;j<dim;j++)
	      for(unsigned int k=0;k<dim;k++)
		for(unsigned int l=0;l<dim;l++){
		  C_ep[i][j][k][l]+=a_elastoplastic[A][B]*nxn[A][i][j]*nxn[B][k][l];
		}
      for(unsigned int A=0;A<dim;A++)
	for(unsigned int i=0;i<dim;i++)
	  for(unsigned int j=0;j<dim;j++)
	    for(unsigned int k=0;k<dim;k++)
	      for(unsigned int l=0;l<dim;l++){
		C_ep[i][j][k][l]+=2*principal_stress[A]*C_tr[A][i][j][k][l];
	      }
       
     
      //C_ep calculated
      C=C_ep;
      for(unsigned int i=0;i<dim;i++)_eigval[i]=eigval[i];
    }
  //plastic condition ends

  
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++){
      history[q]->bLeftcauchyIteration[i][j]=_eigval[0]*nxn[0][i][j]+_eigval[1]*nxn[1][i][j]+_eigval[2]*nxn[2][i][j];
      
      history[q]->defgradIteration(i,j)=Ft(i,j); 
    }
  //evaluate stress=kirchhoff stress
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int A=0;A<dim;A++)
	Kirchhoff_stress[i][j]+=principal_stress[A]*nxn[A][i][j];
  
  //calculate piola stress
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	Piola_stress[i][j]+=Kirchhoff_stress[i][k]*FtInv(j,k);
  
  history[q]->alphaIteration=alpha;
  history[q]->backstressIteration=_backstress;
  history[q]->elasStrain11=gamma;
 
  
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
    
       
    //end quadrature loop
    
  }
}

#endif /* MECHANICS_H_ */
