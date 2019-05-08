/*template<int dim>
void myTable():Table<dim>(){
  myTable(){
    for i, j myTabLe(i,j)=0.0;
  }
  };*/

template<int dim>
void calculate_variation_Fe(Table<4, double>&variation_Fe, FullMatrix<double>&F, FullMatrix<double>&Fp, Table<4, double>&variation_Fp ){
  FullMatrix<double>Fp_inv(dim, dim);
  Table<4, double>varFe2(dim, dim, dim, dim), variation_Fpinv(dim, dim, dim, dim),varFe1(dim, dim, dim, dim) ;
  
  Fp_inv.invert(Fp);
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++){
	  variation_Fe[i][j][k][l]=0.; 
	  varFe2[i][j][k][l]=0.;
	  variation_Fpinv[i][j][k][l]=0.;
	  varFe1[i][j][k][l]=0.;
	}
  
  for(unsigned int J=0;J<dim;J++)
    for(unsigned int A=0;A<dim;A++)
      for(unsigned int I=0;I<dim;I++)
	for(unsigned int l=0;l<dim;l++)
	  for(unsigned int L=0;L<dim;L++)
	    for(unsigned int B=0;B<dim;B++)
	      variation_Fpinv[J][B][l][L]+=Fp_inv(J,A)*variation_Fp[A][I][l][L]*Fp_inv(I,B);
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int A=0;A<dim;A++)
      for(unsigned int l=0;l<dim;l++)
	for(unsigned int L=0;L<dim;L++)
	  varFe1[i][A][l][L]+=(double)(i==l)*Fp_inv(L,A);
 
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int I=0;I<dim;I++)
      for(unsigned int A=0;A<dim;A++)
	for(unsigned int l=0;l<dim;l++)
	  for(unsigned int L=0;L<dim;L++)
	    varFe2[i][A][l][L]+=F(i,I)*variation_Fpinv[I][A][l][L];
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  variation_Fe[i][j][k][l]=varFe1[i][j][k][l]-varFe2[i][j][k][l]; 
  
}


template<int dim>
void calculate_variation_E(Table<4, double>& variation_E, Table<4,double>& variation_Fe,FullMatrix<double>& Fe){
  //variation of green-lagrange strain
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  variation_E[i][j][k][l]=0.;
 
  
  for(unsigned int I=0;I<dim;I++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int l=0;l<dim;l++)
	for(unsigned int L=0;L<dim;L++)
	  for(unsigned int J=0;J<dim;J++)
	    variation_E[I][J][l][L]+=0.5*(variation_Fe[I][j][l][L]*Fe(j,J)+Fe(I,j)*variation_Fe[j][J][l][L]);
}

template<int dim>
void calculate_variation_SPK(Table<4,double>& variation_SPK, Table<4,double>& variation_E, Table<4,double>& ElasticModulii){
  //variation of second piola kirchhoff stress
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  variation_SPK[i][j][k][l]=0.;
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  for(unsigned int m=0;m<dim;m++)
	    for(unsigned int n=0;n<dim;n++)
	      variation_SPK[i][j][m][n]+=ElasticModulii[i][j][k][l]*variation_E[k][l][m][n];
  
}

template<int dim>
void calculate_variation_Bmatrix(Table<2,double>& OriginalSchmidt, Table<4,double>&variation_E, Table<4,double>& variation_B){
  //variation of B=(0.5)*(SchmidtTensor(trandpose)*R_cauchy_green+R-cauchy-green*SchmidtTensor)
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  variation_B[i][j][k][l]=0.;
    
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int m=0;m<dim;m++)
	  for(unsigned int n=0;n<dim;n++)
	    variation_B[i][j][m][n]+=0.5*(OriginalSchmidt[k][i]*variation_E[k][j][m][n]+variation_E[i][k][m][n]*OriginalSchmidt[k][j]);
  //check formulation if this doesn't work
  
}

template<int dim>
void calculate_B_matrix(Table<2,double>&SchmidtTensor, Table<2,double>&R_cauchy_green, Table<2,double>&B){

  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      B[i][j]=0.;
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	B[i][j]+=0.5*(SchmidtTensor[k][i]*R_cauchy_green[k][j]+ R_cauchy_green[i][k]* SchmidtTensor[k][j]);

}


template<int dim>
void calculate_variation_system_matrix(Vector<int>& ActiveSystems,std::vector<Table<2,double> >&OriginalSchmidt,FullMatrix<double>&active_system_matrix,Table<4,double>&variation_system_matrix,Table<4,double>&variation_E,Table<1,double>&shearStress,Table<2,double>&R_cauchy_green,Table<2,double>&SPK,Table<4,double>&ElasticModulii){
  unsigned int size=ActiveSystems.size();
  Table<4,double>variation_SPK(dim, dim, dim, dim);
  calculate_variation_SPK<dim>(variation_SPK,variation_E,ElasticModulii);//variation of second piola kirchhoff stress
 
  for(unsigned int i=0;i<size;i++)
    for(unsigned int j=0;j<size;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++){
	  variation_system_matrix[i][j][k][l]=0.;
	  
	}

  for(unsigned int I=0;I<size;I++){
    unsigned int alpha=ActiveSystems[I];
    for(unsigned int J=0;J<size;J++){
      unsigned int beta=ActiveSystems[J];
      Table<4,double>varRCB1(dim,dim, dim, dim),varRCB2(dim, dim,dim,dim),varBT1(dim, dim, dim, dim),varBT2(dim, dim, dim, dim),variation_B(dim, dim, dim, dim);
      Table<2,double>B(dim, dim);
      for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)B[i][j]=0.;
       for(unsigned int i=0;i<size;i++)
	 for(unsigned int j=0;j<size;j++)
	   for(unsigned int k=0;k<dim;k++)
	     for(unsigned int l=0;l<dim;l++){
		 varRCB1[i][j][k][l]=0.;
		 varRCB2[i][j][k][l]=0.;
		 varBT1[i][j][k][l]=0.;
		 varBT2[i][j][k][l]=0.;
		 variation_B[i][j][k][l]=0.;
	       }
	       
      calculate_B_matrix<dim>(OriginalSchmidt[beta],R_cauchy_green,B);
      calculate_variation_Bmatrix<dim>(OriginalSchmidt[beta],variation_E,variation_B);
      
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int k=0;k<dim;k++)
	  for(unsigned int m=0;m<dim;m++)
	    for(unsigned int n=0;n<dim;n++)
	      for(unsigned int j=0;j<dim;j++)
		for(unsigned int p=0;p<dim;p++)
		  for(unsigned int q=0;q<dim;q++){
		    varRCB1[i][j][m][n]+=2.0*variation_E[i][k][m][n]*ElasticModulii[k][j][p][q]*B[p][q];
		    varRCB2[i][j][m][n]+=R_cauchy_green[i][k]*ElasticModulii[k][j][p][q]*variation_B[p][q][m][n];
		  }
    
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++)
	  for(unsigned int k=0;k<dim;k++)
	    for(unsigned int m=0;m<dim;m++)
	      for(unsigned int n=0;n<dim;n++){
		varBT1[i][j][m][n]+=2.0*(variation_B[i][k][m][n]*SPK[k][j]);
		varBT2[i][j][m][n]+=2.0*(B[i][k]*variation_SPK[k][j][m][n]);
	      }
      
      for(unsigned int i=0;i<dim;i++)
	for(unsigned int j=0;j<dim;j++)
	  for(unsigned int m=0;m<dim;m++)
	    for(unsigned int n=0;n<dim;n++)
	      variation_system_matrix[I][J][m][n]+=(std::abs(shearStress[alpha])/shearStress[alpha])*(std::abs(shearStress[beta])/shearStress[beta])*(varRCB1[i][j][m][n]+varRCB2[i][j][m][n]+varBT1[i][j][m][n]+varBT2[i][j][m][n])*OriginalSchmidt[alpha][i][j];
      //beta for loop end
      
    }
    //alpha for loop end
  }
  
  //function definition end
}

template<int dim>
void calculate_variation_rhs(Table<3, double>&variation_rhs,Table<4, double>&variation_E,Table<4, double>&variation_SPK,Table<2,double>&SPK,Table<2, double>&R_cauchy_green,Vector<int>&ActiveSystems ,std::vector<Table<2, double> >&OriginalSchmidt,Table<1, double>&shearStress){

  unsigned int size=ActiveSystems.size();
  Table<4, double>temp(dim, dim, dim, dim);
  for(unsigned int i=0;i<size;i++)
    for(unsigned int j=0;j<size;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  temp[i][j][k][l]=0.;
  for(unsigned int i=0;i<size;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	variation_rhs[i][j][k]=0.;
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int m=0;m<dim;m++)
	  for(unsigned int n=0;n<dim;n++)
	    temp[i][j][m][n]+=(2.0*variation_E[i][k][m][n]*SPK[k][j]+R_cauchy_green[i][k]*variation_SPK[k][j][m][n]);
  for(unsigned int I=0;I<size;I++){
    unsigned int alpha=ActiveSystems[I];
    for(unsigned int i=0;i<dim;i++)
      for(unsigned int j=0;j<dim;j++)
	for(unsigned int m=0;m<dim;m++)
	  for(unsigned int n=0;n<dim;n++)
    variation_rhs[I][m][n]+=(std::abs(shearStress[alpha])/shearStress[alpha])*temp[i][j][m][n]*OriginalSchmidt[alpha][i][j];
  }
  
}


template<int dim>
void calculate_variation_gamma(Table<3, double>&variation_gamma,Table<2, double>&R_cauchy_green,FullMatrix<double>&active_system_matrix,Table<4,double>&variation_system_matrix,std::vector<Table<2, double> >&OriginalSchmidt, Vector<int>ActiveSystems,Vector<double>gamma_total, Table<2, double>&SPK,Table<4,double>&variation_E, Table<4, double>&variation_SPK, Table<1,double>&shearStress, Table<1,double>&crss){
  
  unsigned int size=ActiveSystems.size();
  FullMatrix<double>active_system_inv(size, size);active_system_inv=0.0;
  for(unsigned int i=0;i<size;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	variation_gamma[i][j][k]=0.;
  active_system_inv.left_invert(active_system_matrix);
  //active_system_inv=active_system_matrix.invert();
  Table<3,double>variation_rhs(size, dim, dim);
  Table<3, double>temp(size, dim, dim),temp1(size, dim, dim);

  for(unsigned int i=0;i<size;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++){
	variation_rhs[i][j][k]=0.;
	temp[i][j][k]=0.;
	temp1[i][j][k]=0.;
      }
  calculate_variation_rhs<dim>(variation_rhs, variation_E, variation_SPK, SPK, R_cauchy_green, ActiveSystems ,OriginalSchmidt,shearStress);
  
  for(unsigned int I=0;I<size;I++){
    //unsigned int alpha=ActiveSystems[I];
    for(unsigned int J=0;J<size;J++){
      unsigned int beta=ActiveSystems[J];
      for(unsigned int m=0;m<dim;m++){
	for(unsigned int n=0;n<dim;n++){
	  temp[I][m][n]+=variation_system_matrix[I][J][m][n]*gamma_total[beta];//*(std::abs(shearStress[beta])/shearStress[beta]);
	  
	}
      }
    }
  }
  for(unsigned int i=0;i<size;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	temp1[i][j][k]=variation_rhs[i][j][k]-temp[i][j][k];
  for(unsigned int I=0;I<size;I++){
    unsigned int beta=ActiveSystems[I];
    for(unsigned int J=0;J<size;J++){
      unsigned int alpha=ActiveSystems[J];
      for(unsigned int m=0;m<dim;m++)
	for(unsigned int n=0;n<dim;n++)
	  variation_gamma[I][m][n]+=active_system_inv(J,I)*temp1[J][m][n];
    }
  }
 
}


template<int dim>
void calculate_variation_Fp(Table<4, double>&variation_Fp,Table<4, double>&variation_E, Table<2, double>&updateFp, Vector<int>& ActiveSystems, std::vector< Table<2,double> >&OriginalSchmidt, FullMatrix<double>&Fp,Table<2, double>&SPK, Table<4, double>&variation_SPK,Table<2, double>&R_cauchy_green ,Vector<double>&gamma_total, FullMatrix<double>&active_system_matrix, Vector<double>&active_system_rhs,Table<1, double>&shearStress, Table<1, double>&crss, Table<4, double>&ElasticModulii){
  Table<4, double> variationtempFp(dim, dim, dim, dim);
  unsigned int size=ActiveSystems.size();
  Table<4,double>variation_system_matrix(size, size, dim, dim),temp2(dim, dim, dim, dim),temp3(dim, dim, dim, dim),temp4(dim, dim, dim, dim);
  Table<3,double>variation_gamma(size, dim, dim);
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++){
	  variationtempFp[i][j][k][l]=0.;
	  temp2[i][j][k][l]=0.;
	  temp3[i][j][k][l]=0.;
	  temp4[i][j][k][l]=0.;
	}
  for(unsigned int i=0;i<size;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	variation_gamma[i][j][k]=0.;
  for(unsigned int i=0;i<size;i++)
    for(unsigned int j=0;j<size;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  variation_system_matrix[i][j][k][l]=0.;
  
  
  calculate_variation_system_matrix<dim>( ActiveSystems, OriginalSchmidt, active_system_matrix, variation_system_matrix, variation_E, shearStress,R_cauchy_green, SPK, ElasticModulii);
  calculate_variation_gamma<dim>(variation_gamma,R_cauchy_green, active_system_matrix, variation_system_matrix, OriginalSchmidt, ActiveSystems, gamma_total, SPK, variation_E, variation_SPK, shearStress, crss);
 
  for(unsigned int I=0;I<size;I++){
    unsigned int alpha=ActiveSystems[I];
    Table<4, double>temp1(dim, dim, dim, dim);
    for(unsigned int i=0;i<dim;i++)
      for(unsigned int j=0;j<dim;j++)
	for(unsigned int k=0;k<dim;k++)
	  for(unsigned int l=0;l<dim;l++)
	    temp1[i][j][k][l]=0.;
    for(unsigned int i=0;i<dim;i++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int m=0;m<dim;m++)
	  for(unsigned int n=0;n<dim;n++)
	    temp1[i][k][m][n]+=(std::abs(shearStress[alpha])/shearStress[alpha])*variation_gamma[I][m][n]*OriginalSchmidt[alpha][i][k];
    
    for(unsigned int i=0;i<dim;i++)
      for(unsigned int j=0;j<dim;j++)
	for(unsigned int k=0;k<dim;k++)
	  for(unsigned int l=0;l<dim;l++)
	    temp2[i][j][k][l]+=temp1[i][j][k][l];
    
  }

  for(unsigned int i=0;i<dim;i++)
    for(unsigned int k=0;k<dim;k++)
      for(unsigned int m=0;m<dim;m++)
	for(unsigned int n=0;n<dim;n++)
	  for(unsigned int j=0;j<dim;j++)
	    temp3[i][j][m][n]+=temp2[i][k][m][n]*Fp(k,j);
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int k=0;k<dim;k++)
      for(unsigned int j=0;j<dim;j++)
	for(unsigned int m=0;m<dim;m++)
	  for(unsigned int n=0;n<dim;n++)
	    temp4[i][j][m][n]+=updateFp[i][k]*variation_Fp[k][j][m][n];
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  variation_Fp[i][j][k][l]=0.;
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  variation_Fp[i][j][k][l]+=temp3[i][j][k][l]+temp4[i][j][k][l];
  
   
  
}


