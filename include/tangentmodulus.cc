//update elastoplastic tangent modulus
template<int dim>
void elastoplastic_tangent(FullMatrix<double>&Fe, FullMatrix<double>&Fp, FullMatrix<double>&F, Table<2, double>&SPK, Table<4, double>& variation_Fe, Table<4, double>&variation_Fp, Table<4, double>&variation_SPK,Table<4, double>&ElastoPlasticModulii, unsigned int currentIncrement, unsigned int currentIteration){
  FullMatrix<double> Fp_inv(dim, dim);
  FullMatrix<double> F_inv(dim, dim);
  F_inv.invert(F);
  Fp_inv.invert(Fp);
  Table<4, double>EM1(dim, dim, dim, dim), EM2(dim, dim, dim, dim),EM3(dim, dim, dim, dim),EM4(dim, dim, dim, dim);
  Table<4, double>variation_F_inv_T(dim, dim, dim, dim);
  //EM1=0.; EM2=0.; EM3=0.; EM4=0.;
   for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++){
	  ElastoPlasticModulii[i][j][k][l]=0.;
	  EM1[i][j][k][l]=0.;
	  EM2[i][j][k][l]=0.;
	  EM3[i][j][k][l]=0.;
	  EM4[i][j][k][l]=0.;
	  variation_F_inv_T[i][j][k][l]=0.;
	}
  
  //compute part 1 of tangent modulus variation
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int A=0;A<dim;A++)
      for(unsigned int l=0;l<dim;l++)
	for(unsigned int L=0;L<dim;L++)
	  for(unsigned int B=0;B<dim;B++)
	    for(unsigned int j=0;j<dim;j++)
	      for(unsigned int J=0;J<dim;J++)
		EM1[i][J][l][L]+=variation_Fe[i][A][l][L]*SPK[A][B]*Fe(B,j)*F_inv(j,J);
  	  

  

  //compute part 1 of tangent modulus variation
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int A=0;A<dim;A++)
      for(unsigned int B=0;B<dim;B++)
	for(unsigned int l=0;l<dim;l++)
	  for(unsigned int L=0;L<dim;L++)
	    for(unsigned int j=0;j<dim;j++)
	      for(unsigned int J=0;J<dim;J++)
		EM2[i][J][l][L]+=Fe(i,A)*variation_SPK[A][B][l][L]*Fe(B,j)*F_inv(j,J);
  /*if(currentIncrement==13 && currentIteration==1){
    //for(unsigned int i=0;i<dim;i++){
    //for(unsigned int j=0;j<dim;j++){
	for(unsigned int k=0;k<dim;k++){
	  for(unsigned int l=0;l<dim;l++){
	    std::cout<<Fe/[k][l]<<" ";
	  }std::cout<<"\n";
	}
	//  }
  //}std::cout<<"\n";
    
  }*/
  
  
  //compute part 1 of tangent modulus variation
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int A=0;A<dim;A++)
      for(unsigned int B=0;B<dim;B++)
	for(unsigned int j=0;j<dim;j++)
	  for(unsigned int l=0;l<dim;l++)
	    for(unsigned int L=0;L<dim;L++)
	      for(unsigned int J=0;J<dim;J++)
		EM3[i][J][l][L]+=Fe(i,A)*SPK[A][B]*variation_Fe[B][j][l][L]*F_inv(j,J);
		  
  //calculate variation of F_inv_transpose
  
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int J=0;J<dim;J++)      
      for(unsigned int l=0;l<dim;l++)
	for(unsigned int L=0;L<dim;L++)
	  variation_F_inv_T[i][J][l][L]+=(-1.0)*F_inv(i,l)*F_inv(J,L);
  
  
  //compute part 1 of tangent modulus variation
  for(unsigned int i=0;i<dim;i++)
    for(unsigned int A=0;A<dim;A++)
      for(unsigned int B=0;B<dim;B++)
	for(unsigned int j=0;j<dim;j++)
	  for(unsigned int l=0;l<dim;l++)
	    for(unsigned int L=0;L<dim;L++)
	      for(unsigned int J=0;J<dim;J++)
		EM4[i][J][l][L]+=Fe(i,A)*SPK[A][B]*Fe(B,j)*variation_F_inv_T[j][J][l][L];

  for(unsigned int i=0;i<dim;i++)
    for(unsigned int j=0;j<dim;j++)
      for(unsigned int k=0;k<dim;k++)
	for(unsigned int l=0;l<dim;l++)
	  ElastoPlasticModulii[i][j][k][l]=EM1[i][j][k][l]+EM2[i][j][k][l]+EM3[i][j][k][l]+EM4[i][j][k][l];

   
}
