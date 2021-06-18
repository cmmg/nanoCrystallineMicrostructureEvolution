//new
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2012
//authors: rudraa (2012, 2018)
//

#ifndef CHEMO_H_
#define CHEMO_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"



template<int dim>
void evaluateFieldAtQuadraturePoint(unsigned int q, Table<1, double >&phi, Table<1, double >&phi_conv, Table<2, double >&phi_j, double&sol, double&sol_conv, Table<1, double >&sol_j,double&mu ,Table<1, double >&mu_j,FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,double >& ULocal, dealii::Table<1, double>& ULocalConv, Vector<double>& R,FullMatrix<double>&local_matrix,unsigned int currentIncrement,unsigned int currentIteration ,std::vector<historyVariables<dim>* >& history, double &freeEnergyChemBulk, double & freeEnergyChemGB) {
  int var=0; if(isMechanics)var=dim;
  unsigned int dofs_per_cell=fe_values.dofs_per_cell;
  for(unsigned int i=0;i<dofs_per_cell;i++){
    
    int ci=fe_values.get_fe().system_to_component_index(i).first - var;
    if(ci>=0 && ci<n_diff_grains){
      phi[ci]+=fe_values.shape_value(i,q)*ULocal[i];
      phi_conv[ci]+=fe_values.shape_value(i,q)*ULocalConv[i];
      for(unsigned int j=0;j<dim;j++){
	phi_j[ci][j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];
      }
    }
    if(ci==n_diff_grains){
      sol+=fe_values.shape_value(i,q)*ULocal[i];
      sol_conv+=fe_values.shape_value(i,q)*ULocalConv[i];
      for(unsigned int j=0;j<dim;j++){
	sol_j[j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];
      }
    }
    if(ci==n_diff_grains+1){
      mu+=fe_values.shape_value(i,q)*ULocal[i];
      for(unsigned int j=0;j<dim;j++){
	mu_j[j]+=fe_values.shape_grad(i,q)[j]*ULocal[i];
      }
    }

    //dof loop ends
  }

  //std::cout<<"\nsolute concentration\t"<<sol<<"\n";
  //std::cout<<"chemical potential\t"<<mu<<"\n";
  //exit(-1);
  double phiSQ=0.;
  for(unsigned int N=0;N<n_diff_grains;N++){
    phiSQ+=phi[N]*phi[N];
    freeEnergyChemBulk+=(4./3.)*(1.0-4.0* phi[N]*phi[N]*phi[N])*fe_values.JxW(q);
  }
  freeEnergyChemBulk+=4.0* phiSQ*phiSQ*fe_values.JxW(q);

  for(unsigned int N=0;N<n_diff_grains;N++){
    for(unsigned int j=0;j<dim;j++){
      freeEnergyChemGB+=(InterfaceEnergyParameter/2.)*(phi_j[N][j]*phi_j[N][j])*fe_values.JxW(q);   
    }
  }


  // definition ends
}

template<int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1,double >& ULocal, dealii::Table<1, double>& ULocalConv, Vector<double>& R,FullMatrix<double>&local_matrix,FullMatrix<double>& jacobian,unsigned int currentIncrement,  unsigned int currentIteration ,std::vector<historyVariables<dim>* >& history, double & freeEnergyChemBulk, double & freeEnergyChemGB){
 
  unsigned int dofs_per_cell=fe_values.dofs_per_cell;
  unsigned int n_q_points=fe_values.n_quadrature_points;
  double epsilon=InterfaceEnergyParameter;
  double M_sol=M_alpha, M_phi=Mobility;
  for(unsigned int q=0;q<n_q_points;q++){
    Table<1, double >phi(n_diff_grains), phi_conv(n_diff_grains), sol_j(dim), mu_j(dim);
    Table<2,double >phi_j(n_diff_grains, dim);
    
    for(unsigned int i=0;i<n_diff_grains;i++){phi[i]=0.; phi_conv[i]=0.;}
    for(unsigned int i=0;i<dim;i++){sol_j[i]=0.; mu_j[i]=0.;}
    for(unsigned int i=0;i<n_diff_grains;i++){
      for(unsigned int j=0;j<dim;j++){
	phi_j[i][j]=0.;
      }
    }
   
    double mu=0.0, sol=0.0, sol_conv=0.0;
    evaluateFieldAtQuadraturePoint<dim>( q, phi, phi_conv, phi_j, sol, sol_conv, sol_j, mu , mu_j, fe_values,  DOF,fe_face_values,cell, dt,ULocal,  ULocalConv,  R,local_matrix,  currentIncrement,currentIteration, history, freeEnergyChemBulk, freeEnergyChemGB) ;
    if(!isSoluteDrag){
      //if(currentIncrement>=50 && currentIncrement<150){M_phi=0.0;}
      /*std::cout<<"phi values";
      for(unsigned int i=0;i<n_diff_grains;i++){
	std::cout<<phi[i]<<"\t";
      }exit(-1);*/
      for(unsigned int i=0;i<dofs_per_cell;i++){
	int ci=fe_values.get_fe().system_to_component_index(i).first - DOF;
	
	if(ci>=0 && ci<n_diff_grains){
	  double phi2_sum=0.;
	  double W_sol=0.;
	  W_sol=1.0 ;
	  for(unsigned int I=0;I<n_diff_grains;I++){
	    phi2_sum+=pow(phi[I],2);
	  }
	  
	  R[i]+=(1/dt)*fe_values.shape_value(i,q)*(phi[ci]-phi_conv[ci])*fe_values.JxW(q);
	  std::cout<<R[i]<<"\t";
	  R[i]+=fe_values.shape_value(i,q)*(4./3.)*(M_phi)*(-12.0*phi[ci]*phi[ci]+12.0*phi[ci]*phi2_sum)*fe_values.JxW(q);
	  std::cout<<R[i]<<"\t";
	  for(unsigned int j=0;j<dim;j++){
	    R[i]+=(M_phi)*epsilon*fe_values.shape_grad(i,q)[j]*phi_j[ci][j]*fe_values.JxW(q);
	    std::cout<<R[i]<<"\t";
	  }
	  
	}
	
      }
      
      double g_phi=0., phi2_sum=0., phi3_sum=0.;
      for(unsigned int N=0;N<n_diff_grains;N++){
	phi2_sum+=pow(phi[N],2);
	phi3_sum+=pow(phi[N],3);
      }
      
      for(unsigned int A=0;A<dofs_per_cell;A++){
	for(unsigned int B=0;B<dofs_per_cell;B++){
	  int ca=fe_values.get_fe().system_to_component_index(A).first- DOF;
	  int cb=fe_values.get_fe().system_to_component_index(B).first - DOF;
	  if(ca>=0 && ca<n_diff_grains ){
	    if(cb>=0 && cb<n_diff_grains){
	      if(ca==cb){
		
		local_matrix(A,B)+=fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*(1.0/dt)*fe_values.JxW(q);
		std::cout<<"\n\n"<<local_matrix(A,B)<<"\t";
		std::cout<<fe_values.shape_value(A,q)<<" "<<fe_values.shape_value(B,q)<<" "<<1.0/dt<<" "<<fe_values.JxW(q)<<"\t\t";
		local_matrix(A,B)+=M_phi*fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*(4./3.)*(-24.0*phi[ca]+12.0*phi2_sum+24.0*phi[ca]*phi[ca])*fe_values.JxW(q);
		std::cout<<local_matrix(A,B)<<"\t";
		for(unsigned int i=0;i<dim;i++){
		  local_matrix(A,B)+=epsilon*M_phi*fe_values.shape_grad(A,q)[i]*fe_values.shape_grad(B,q)[i]*fe_values.JxW(q);
		  std::cout<<local_matrix(A,B)<<"\t";
		}
		
	      }
	      else{
		local_matrix(A,B)+=M_phi*fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*(96.0/3.0)*(phi[ca]*phi[cb])*fe_values.JxW(q);
		std::cout<<local_matrix(A,B)<<"\t";
		exit(-1);
	      }
	      
	    }//if cj ends
	    
	  }//if ci ends
	  
	  
	  
	}//dofs B ends
	
	
      }//dofs A ends
      
      for(unsigned int i=0;i<dofs_per_cell;i++){
	for(unsigned int j=0;j<dofs_per_cell;j++){
	  std::cout<<local_matrix(i,j)<<" ";
	}std::cout<<"\n";
      }std::cout<<"\n\n";
      for(unsigned int i=0;i<dofs_per_cell;i++){
	std::cout<<R[i]<<" ";
      }
      exit(-1);

    } // !isSoluteDrag statement ends here

    if(isSoluteDrag){
      //std::cout<<"control here in !isSoluteDrag";exit(-1);
      double W_sol=0.;
      W_sol=WA*(1.-sol) + WB*(sol);
      if(currentIncrement<=40){W_sol=1.0; M_sol=0.0;}
      //if(currentIncrement>12500){M_sol=0.0;}
      if(currentIncrement>40 && currentIncrement<=2500){M_phi=0.0;}
      //if(currentIncrement>=5 && currentIncrement<150){M_phi=0.0;}
      for(unsigned int A=0;A<dofs_per_cell;A++){
	unsigned int ca=fe_values.get_fe().system_to_component_index(A).first- DOF;
	if(ca>=0 && ca<n_diff_grains){
	  double phi2sum=0.0;
	  for(unsigned int I=0;I<n_diff_grains;I++){
	    phi2sum+= phi[I]*phi[I];
	  }
	  R[A]+=(1./dt)* fe_values.shape_value(A,q)*(phi[ca]-phi_conv[ca])*fe_values.JxW(q);
	  R[A]+=M_phi*W_sol* fe_values.shape_value(A,q)*(4./3.)*(-12.0*phi[ca]*phi[ca]+12.0* phi2sum*phi[ca])*fe_values.JxW(q);
	  for(unsigned int j=0;j<dim;j++){
	    R[A]+=M_phi*W_sol*epsilon*fe_values.shape_grad(A,q)[j]*phi_j[ca][j]*fe_values.JxW(q);
	  }
	}// ca>=0 && ca<n_diff_grains block ends here
	
	if(ca==n_diff_grains){
	  R[A]+=(1./dt)* (sol-sol_conv)*fe_values.shape_value(A,q)*fe_values.JxW(q);
	  for(unsigned int j=0;j<dim;j++){
	    R[A]+=M_sol*fe_values.shape_grad(A,q)[j]*mu_j[j]*fe_values.JxW(q);
	  }
	}//ca == n_diff_grains ends here

	if(ca==n_diff_grains+1){
	  double g_phi=0., phi2Sum=0., phi3Sum=0.;
	  for(unsigned int i=0;i<n_diff_grains;i++){
	    phi2Sum+=phi[i]*phi[i];
	    phi3Sum+=phi[i]*phi[i]*phi[i];
	  }
	  g_phi=(4./3.)*(1.0- 4.0* phi3Sum + 3.0 * phi2Sum*phi2Sum);
	  //solute energy double well G
	  double dG_dsol= 800*sol*(sol-1.0)*(sol-0.5)+ g_phi*(WB-WA);//(-4.0*sol-5);
	  R[A]+=fe_values.shape_value(A,q)* (mu-dG_dsol)*fe_values.JxW(q);
	  for(unsigned int j=0;j<dim;j++){
	    R[A]-=kappa*fe_values.shape_grad(A,q)[j]*sol_j[j]*fe_values.JxW(q);
	  }

	}//ca==n_diff_grains+1 ends here

      }// dofs A loop ends here

      // jacobian matrix calculation starts here
      for(unsigned int A=0;A<dofs_per_cell;A++){
	unsigned int ca=fe_values.get_fe().system_to_component_index(A).first- DOF;
	for(unsigned int B=0;B<dofs_per_cell;B++){
	  unsigned int cb=fe_values.get_fe().system_to_component_index(B).first- DOF;
	  
	  if(ca>=0 && ca<n_diff_grains){
	    if(cb>=0 && cb<n_diff_grains){
	      if(ca==cb){
		local_matrix(A,B)+=(1.0/dt)*fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*fe_values.JxW(q);
		double fValue=0.0, phi2sum=0.0;
		for (unsigned int i=0;i<n_diff_grains;i++){
		  phi2sum+=phi[i]*phi[i];//here                                                     
		}
		fValue=(4./3.)* (-24.0*phi[ca]+24.0*phi[ca]*phi[ca] + 12.0* phi2sum  );
		local_matrix(A,B)+=M_phi* W_sol * fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*fValue*fe_values.JxW(q);
		for(unsigned int j=0;j<dim;j++){
		  local_matrix(A,B)+=M_phi * epsilon * W_sol * fe_values.shape_grad(A,q)[j]*fe_values.shape_grad(B,q)[j]*fe_values.JxW(q);
		}
	      }//ca==cb ends here
	      else{
		local_matrix(A,B)+=M_phi*W_sol * fe_values.shape_value(A,q)* fe_values.shape_value(B,q)* (32.0)* phi[ca]*phi[cb]*fe_values.JxW(q);
	      }//ca!=cb ends here
	    }//cb<n_diff_grains ends here
	    
	    if(cb==n_diff_grains){
	      double phi2sum=0.0;
	      for(unsigned int i=0;i<n_diff_grains;i++){
		phi2sum+=phi[i]*phi[i];
	      }
	      local_matrix(A,B)+=M_phi * (WB - WA) * fe_values.shape_value(A,q)* fe_values.shape_value(B,q)*(4.0/3.0)* (-12.0 *phi[ca]*phi[ca]+ 12.0 * phi[ca]*phi2sum )*fe_values.JxW(q);
	      for(unsigned int j=0;j<dim;j++){
		local_matrix(A,B)+=M_phi * epsilon *fe_values.shape_value(B,q)* (WB-WA) *fe_values.shape_grad(A,q)[j]*phi_j[ca][j]*fe_values.JxW(q) ;
	      }

	    }//cb==n_diff_grains
	    
	    if(cb==n_diff_grains+1){
	      local_matrix(A,B)+=0.0;
	    }//cb==n_diff_grains+1
	    
	  }//ca<n_diff_grains ends here
	  
	  if(ca==n_diff_grains){
	    if(cb>=0 && cb<n_diff_grains){
	      local_matrix(A,B)+=0.0;
	    }//cb>=0 && cb<n_diff_grains ends here
	    if(cb==n_diff_grains){
	      local_matrix(A,B)+=(1.0/dt)* fe_values.shape_value(A,q) * fe_values.shape_value(B,q) * fe_values.JxW(q);
	    }//cb==n_diff_grains ends here
	    if(cb==n_diff_grains+1){
	      for(unsigned int j=0;j<dim;j++){
		local_matrix(A,B)+=M_sol * fe_values.shape_grad(A,q)[j]* fe_values.shape_grad(B,q)[j]*fe_values.JxW(q) ;
	      }
	    }//cb==n_diff_grains+1
	  
	  }//ca==n_diff_grains ends here

	  if(ca==n_diff_grains+1){
	    if(cb>=0 && cb<n_diff_grains){
	      double phi2sum=0.0, value=0.0;
	      for(unsigned int i=0;i<n_diff_grains;i++){
		phi2sum+=phi[i]*phi[i];
	      }
	      value=(4.0/3.0) * (-12.0 * phi[cb]*phi[cb] + 12.0 * phi[cb]*phi2sum);
	      local_matrix(A,B)-=fe_values.shape_value(A,q)*fe_values.shape_value(B,q) * value * (WB - WA)*fe_values.JxW(q);//here jacobian changed from + to -
	    }//cb<n_diff_grains ends here
	    if(cb==n_diff_grains){
	      double value=0.0;
	      value=800 * (sol * (sol-0.5) + sol*(sol-1.0) + (sol-0.5)*(sol-1.0));
	      local_matrix(A,B)+= (-1.0)* fe_values.shape_value(A,q)*fe_values.shape_value(B,q) * value * fe_values.JxW(q);
	      for(unsigned int j=0;j<dim;j++){
		local_matrix(A,B)+=(-1.0)*kappa * fe_values.shape_grad(A,q)[j] * fe_values.shape_grad(B,q)[j]*fe_values.JxW(q);
	      }
	    }//cb==n_diff_grains ends here
	    
	    if(cb==n_diff_grains+1){
	      local_matrix(A,B)+=fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*fe_values.JxW(q);
	    }//cb==n_diff_grains
	    
	  }//ca==n_diff_grains+1 ends here
	  
	  
	  
	}//dofs B loop ends here
      }// dofs A loop ends here

    }// isSoluterag block ends here


  }//qudrature loop ends here
    
}

#endif /* CHEMO_H_ */

