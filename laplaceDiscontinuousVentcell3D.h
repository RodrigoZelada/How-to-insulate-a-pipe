void laplaceDiscontinuousVentcell3D(Mat A, Vec b, mesh3D &Th, vector<double> &kappa, vector<vector<double>> &u, map<int,int> Gamma, mesh3D &Th2, double (*f)(double, double, double), double (*TD)(double, double, double), DirichletBC &DBC1, DirichletBC &DBC2)
{
	int N = Th.ntet;
	int Nt = Th.nt;

	int D = Th.d; //dimension 2

	//Integers
	int n;
	int i, j, ii, jj, jj_tilde;

	double volume, area;

	vector<double> coefii(4), coefjj(4), coefjj_tilde(4);
	vector<double> normal(3),normalOmegai(3),gradii(3),gradjj(3),gradjj_tilde(3);

	int n_tilde, j_tilde, itet, jtet;
	vector<vector<double>> coef = baseSpace(Th);
	vector<vector<double>> coef2 = baseSpace(Th2);
	double kRegion, sign; 
	int region;

	PetscInt ii_global, jj_global, jj_global_tilde;
	PetscReal val;

	int Nv1 = DBC1.nvNotDBC;
	int iibelongDBC,jjbelongDBC,jjbelongDBC_tilde;

	//Iterate each tetrahedra: T_n
	for (n = 0; n < N; n++){
		volume = Th.volume[n];
		region = Th.tetrahedra[n][D+1]; 
		kRegion = k1;
		if (region == region2){
			kRegion = k2;
		}

		//volume terms
		for (i = 0; i <= D; i++){

			ii = Th.tetrahedra[n][i];
			for (int k=0;k<=D;k++){
				coefii[k] = coef[n][4*i+k];
				if (k<Th.d){gradii[k] = coefii[k];}
			}

			if (region == region1){
				iibelongDBC = DBC1.vertices[ii][1]; 
				ii_global = DBC1.vertices[ii][0];
			}
			else { //region2
				iibelongDBC = DBC2.vertices[ii][1];
				ii_global = DBC2.vertices[ii][0]+ Nv1; 
			}

			//Iterate a pair of vertex:
			for (j = 0; j <= D; j++){
				jj = Th.tetrahedra[n][j];
				for (int k=0;k<=D;k++){
					coefjj[k] = coef[n][4*j+k];
					if (k<Th.d){gradjj[k] = coefjj[k];}
				}

				if (region == region1){
					jjbelongDBC = DBC1.vertices[jj][1];
					jj_global = DBC1.vertices[jj][0]; 
				}
				else { //region2
					jjbelongDBC = DBC2.vertices[jj][1];
					jj_global = DBC2.vertices[jj][0] + Nv1; 
				}

				val = 0.;
				//a(u,v) = grad(u)*grad(v),
				val += volume*kRegion*dot(gradii,gradjj);

				//Exact for P2
				if (region == region1){
					for (int p=0; p<=D; p++){
						int pp=Th.tetrahedra[n][p];
						//xp1 = point pp + component 1 = {xpp,ypp,zpp,1.};
						vector<double> xp1 = {Th.vertices[pp][0],Th.vertices[pp][1],Th.vertices[pp][2],1.};
						vector<double> up = {u[pp][0],u[pp][1],u[pp][2]};
						val -= rhocp*(volume/20.)*dot(coefii,xp1)*dot(gradjj,up);


						for (int q=0; q<p; q++){
							int qq=Th.tetrahedra[n][q];
							vector<double> xq1 = {Th.vertices[qq][0],Th.vertices[qq][1],Th.vertices[qq][2],1.};
							vector<double> uq = {u[qq][0],u[qq][1],u[qq][2]};
							vector<double> xpq1 = {0.5*xp1[0] + 0.5*xq1[0],0.5*xp1[1] + 0.5*xq1[1],0.5*xp1[2] + 0.5*xq1[2],1.};
							vector<double> upq = {0.5*up[0] + 0.5*uq[0],0.5*up[1] + 0.5*uq[1],0.5*up[2] + 0.5*uq[2]};
							val += rhocp*(volume/5.)*dot(coefii,xpq1)*dot(gradjj, upq);
						}
					}
				}
				if ( (iibelongDBC == 0) && (jjbelongDBC==0) ){
					MatSetValues(A,1,&ii_global,1,&jj_global,&val,ADD_VALUES);
				}
				
				if ( (iibelongDBC == 0) && (jjbelongDBC==1) ){ //if j is Dirichlet
					val=-val*TD(Th.vertices[jj][0],Th.vertices[jj][1],Th.vertices[jj][2]);
					VecSetValues(b,1,&ii_global,&val,ADD_VALUES);
				}

			}
		}
	}
    for (int t=0; t<Nt; t++){
        int label_triangle = Th.triangles[t][D];
        n = Th.triangleToTetra[t];
		area = Th.area[t];

		for (int i=0; i<D; i++){
			normalOmegai[i] = Th.normal[t][i];
		}

		switch (label_triangle)
		{

			case 10: //it belongs to region2
				for (i=0; i<D; i++){
					ii = Th.triangles[t][i];
		
					for (int p=0; p<=D; p++){
						if (ii == Th.tetrahedra[n][p]){itet = p;}
					}
		
					for (int k=0;k<=D;k++){
						coefii[k] = coef[n][4*itet+k];
						if (k<Th.d){gradii[k] = coefii[k];}
					}

					ii_global = DBC2.vertices[ii][0] + Nv1; iibelongDBC = DBC2.vertices[ii][1];
					val = 0.;
					for (int p=0; p<D; p++){
						int pp = Th.triangles[t][p];
						vector<double>  xp1 = {Th.vertices[pp][0], Th.vertices[pp][1],Th.vertices[pp][2],1.};
						for (int q=0; q<p; q++){
							int qq=Th.triangles[t][q];
							vector<double>  xq1 = {Th.vertices[qq][0], Th.vertices[qq][1], Th.vertices[qq][2],1.};
							vector<double>  xpq1 = {0.5*xp1[0] + 0.5*xq1[0], 0.5*xp1[1] + 0.5*xq1[1], 0.5*xp1[2] + 0.5*xq1[2],1.};
							val += (area/3.)*alpha*dot(coefii,xpq1)*f(xpq1[0],xpq1[1],xpq1[2]);
						}
					}
					if ( iibelongDBC == 0 ){
						VecSetValues(b,1,&ii_global,&val,ADD_VALUES); 
					}

					for (j=0; j<D; j++){
						jj = Th.triangles[t][j];
						jj_global = DBC2.vertices[jj][0] + Nv1; jjbelongDBC = DBC2.vertices[jj][1];

						for (int p=0; p<=D; p++){
							if (jj == Th.tetrahedra[n][p]){jtet = p;}
						}

						for (int k=0;k<=D;k++){
							coefjj[k] = coef[n][4*jtet+k];
							if (k<Th.d){gradjj[k] = coefjj[k];}
						}

						val = 0.;
						for (int p=0; p<D; p++){
							int pp = Th.triangles[t][p];
							vector<double> xp1 = {Th.vertices[pp][0], Th.vertices[pp][1], Th.vertices[pp][2],1.};

							for (int q=0; q<p; q++){
								int qq = Th.triangles[t][q];
								vector<double> xq1 = {Th.vertices[qq][0], Th.vertices[qq][1], Th.vertices[qq][2],1.};
								vector<double> xpq1 = {0.5*xq1[0] + 0.5*xp1[0],0.5*xq1[1] + 0.5*xp1[1],0.5*xq1[2] + 0.5*xp1[2],1.};
								val += (area/3.)*alpha*dot(coefii,xpq1)*dot(coefjj,xpq1);
							}
						}
						if ( (iibelongDBC == 0) && (jjbelongDBC==0) ){
							MatSetValues(A,1,&ii_global,1,&jj_global,&val,ADD_VALUES);
						}
					}
				}
				break;
			
			case 120: //it belongs to region2
				for (i=0; i<D; i++){
					ii = Th.triangles[t][i];
		
					for (int p=0; p<=D; p++){
						if (ii == Th.tetrahedra[n][p]){itet = p;}
					}
		
					for (int k=0;k<=D;k++){
						coefii[k] = coef[n][4*itet+k];
						if (k<Th.d){gradii[k] = coefii[k];}
					}

					ii_global = DBC2.vertices[ii][0] + Nv1; iibelongDBC = DBC2.vertices[ii][1];

					val = 0.;
					for (int p=0; p<D; p++){
						int pp = Th.triangles[t][p];
						vector<double> xp1 = {Th.vertices[pp][0], Th.vertices[pp][1],Th.vertices[pp][2],1.};
						for (int q=0; q<p; q++){
							int qq=Th.triangles[t][q];
							vector<double> xq1 = {Th.vertices[qq][0], Th.vertices[qq][1], Th.vertices[qq][2],1.};
							vector<double> xpq1 = {0.5*xp1[0] + 0.5*xq1[0],0.5*xp1[1] + 0.5*xq1[1],0.5*xp1[2] + 0.5*xq1[2],1.};
							val += (area/3.)*alpha*dot(coefii,xpq1)*f(xpq1[0],xpq1[1],xpq1[2]);
						}
					}
					if ( iibelongDBC == 0){
						VecSetValues(b,1,&ii_global,&val,ADD_VALUES); 
					}

					for (j=0; j<D; j++){
						jj = Th.triangles[t][j];
						jj_global = DBC2.vertices[jj][0] + Nv1; jjbelongDBC = DBC2.vertices[jj][1];

						for (int p=0; p<=D; p++){
							if (jj == Th.tetrahedra[n][p]){jtet = p;}
						}

						for (int k=0;k<=D;k++){
							coefjj[k] = coef[n][4*jtet+k];
							if (k<Th.d){gradjj[k] = coefjj[k];}
						}

						val = 0.;
						for (int p=0; p<D; p++){
							int pp = Th.triangles[t][p];
							vector<double> xp1 = {Th.vertices[pp][0], Th.vertices[pp][1], Th.vertices[pp][2],1.};

							for (int q=0; q<p; q++){
								int qq = Th.triangles[t][q];
								vector<double> xq1 = {Th.vertices[qq][0], Th.vertices[qq][1], Th.vertices[qq][2],1.};
								vector<double> xpq1 = {0.5*xq1[0] + 0.5*xp1[0],0.5*xq1[1] + 0.5*xp1[1],0.5*xq1[2] + 0.5*xp1[2],1.};
								val += (area/3.)*alpha*dot(coefii,xpq1)*dot(coefjj,xpq1);
							}
						}
						if ( (iibelongDBC == 0) && (jjbelongDBC==0) ){
							MatSetValues(A,1,&ii_global,1,&jj_global,&val,ADD_VALUES);
						}
						if ( (iibelongDBC == 0) && (jjbelongDBC==1) ){ //if j is Dirichlet
							val=-val*TD(Th.vertices[jj][0],Th.vertices[jj][1],Th.vertices[jj][2]);
							VecSetValues(b,1,&ii_global,&val,ADD_VALUES);
						}
					}
				}
				break;

			case 130: //it belongs to region2
				for (i=0; i<D; i++){
					ii = Th.triangles[t][i];
		
					for (int p=0; p<=D; p++){
						if (ii == Th.tetrahedra[n][p]){itet = p;}
					}
		
					for (int k=0;k<=D;k++){
						coefii[k] = coef[n][4*itet+k];
						if (k<Th.d){gradii[k] = coefii[k];}
					}

					ii_global = DBC2.vertices[ii][0] + Nv1; iibelongDBC = DBC2.vertices[ii][1];
					val = 0.; 
					for (int p=0; p<D; p++){
						int pp = Th.triangles[t][p];
						vector<double> xp1 = {Th.vertices[pp][0], Th.vertices[pp][1],Th.vertices[pp][2],1.};
						for (int q=0; q<p; q++){
							int qq=Th.triangles[t][q];
							vector<double> xq1 = {Th.vertices[qq][0], Th.vertices[qq][1], Th.vertices[qq][2],1.};
							vector<double> xpq1 = {0.5*xp1[0] + 0.5*xq1[0],0.5*xp1[1] + 0.5*xq1[1],0.5*xp1[2] + 0.5*xq1[2],1.};
							val += (area/3.)*alpha*dot(coefii,xpq1)*f(xpq1[0],xpq1[1],xpq1[2]);
						}
					}
					if ( iibelongDBC == 0 ){
						VecSetValues(b,1,&ii_global,&val,ADD_VALUES); 
					}

					for (j=0; j<D; j++){
						jj = Th.triangles[t][j];
						jj_global = DBC2.vertices[jj][0] + Nv1; jjbelongDBC = DBC2.vertices[jj][1];

						for (int p=0; p<=D; p++){
							if (jj == Th.tetrahedra[n][p]){jtet = p;}
						}

						for (int k=0;k<=D;k++){
							coefjj[k] = coef[n][4*jtet+k];
							if (k<Th.d){gradjj[k] = coefjj[k];}
						}

						val = 0.;
						for (int p=0; p<D; p++){
							int pp = Th.triangles[t][p];
							vector<double> xp1 = {Th.vertices[pp][0], Th.vertices[pp][1], Th.vertices[pp][2],1.};

							for (int q=0; q<p; q++){
								int qq = Th.triangles[t][q];
								vector<double> xq1 = {Th.vertices[qq][0], Th.vertices[qq][1], Th.vertices[qq][2],1.};
								vector<double> xpq1 = {0.5*xq1[0] + 0.5*xp1[0],0.5*xq1[1] + 0.5*xp1[1],0.5*xq1[2] + 0.5*xp1[2],1.};
								val += (area/3.)*alpha*dot(coefii,xpq1)*dot(coefjj,xpq1);
							}
						}
						if ( (iibelongDBC == 0) && (jjbelongDBC==0) ){
							MatSetValues(A,1,&ii_global,1,&jj_global,&val,ADD_VALUES);
						}
						if ( (iibelongDBC == 0) && (jjbelongDBC==1) ){ //if j is Dirichlet
							val=-val*TD(Th.vertices[jj][0],Th.vertices[jj][1],Th.vertices[jj][2]);
							VecSetValues(b,1,&ii_global,&val,ADD_VALUES);
						}
					}
				}
				break;
			
			case 100: //label_interface
				if (region == region1){ 
					normal[0]=normalOmegai[0];normal[1]=normalOmegai[1];normal[2]=normalOmegai[2];
					sign = 1.;
				}
				else{
					normal[0]=-normalOmegai[0];normal[1]=-normalOmegai[1];normal[2]=-normalOmegai[2];
					sign = -1.;
				}
				
				for (i=0; i<D; i++){
					ii = Th.triangles[t][i];

					for (int p=0; p<=D; p++){
						if (ii == Th.tetrahedra[n][p]){itet = p;}
					}

					for (int k=0;k<=D;k++){
						coefii[k] = coef[n][4*itet+k];
						if (k<D){gradii[k] = coefii[k];}
					}

					int t_tilde = Gamma.at(t);
					n_tilde = Th2.triangleToTetra[t_tilde];

					if (region == region1){ 
						iibelongDBC = DBC1.vertices[ii][1]; 
						ii_global = DBC1.vertices[ii][0]; 
					}
					else{
						iibelongDBC = DBC2.vertices[ii][1]; 
						ii_global = DBC2.vertices[ii][0]+Nv1; 
					}

					for (j=0; j<D; j++){
						jj = Th.triangles[t][j];
						jj_tilde = Th2.triangles[t_tilde][j];

						if (region == region1){ 
							jjbelongDBC = DBC1.vertices[jj][1];  jjbelongDBC_tilde = DBC2.vertices[jj_tilde][1]; 
							jj_global = DBC1.vertices[jj][0]; jj_global_tilde = DBC2.vertices[jj_tilde][0] + Nv1;
						}
						else{
							jjbelongDBC = DBC2.vertices[jj][1]; jjbelongDBC_tilde = DBC1.vertices[jj_tilde][1];
							jj_global = DBC2.vertices[jj][0]+Nv1;jj_global_tilde =  DBC1.vertices[jj_tilde][0];
						}

						for (int p=0; p<=D; p++){
							if (jj == Th.tetrahedra[n][p]){jtet = p;}
							if (jj_tilde == Th2.tetrahedra[n_tilde][p]){j_tilde = p;}
						}

						for (int k=0;k<=D;k++){
							coefjj[k] = coef[n][4*jtet+k];
							coefjj_tilde[k] = coef2[n_tilde][4*j_tilde+k];
							if (k<D){gradjj[k] = coefjj[k]; gradjj_tilde[k] = coefjj_tilde[k];}
						}

						//Exact for P2
						val = 0.;
						for (int p=0; p<D; p++){
							int pp = Th.triangles[t][p];
							vector<double> xp1 = {Th.vertices[pp][0], Th.vertices[pp][1], Th.vertices[pp][2],1.};

							for (int q=0; q<p; q++){
								int qq = Th.triangles[t][q];
								vector<double> xq1 = {Th.vertices[qq][0], Th.vertices[qq][1], Th.vertices[qq][2],1.};
								vector<double> xpq1 = {0.5*xq1[0] + 0.5*xp1[0],0.5*xq1[1] + 0.5*xp1[1],0.5*xq1[2] + 0.5*xp1[2],1.};
			
								//ks*[T]*[S] = T1*S1 + T2*S2 - T2*S1 - T1*S2 
								val += (ks/epsilon)*(area/3.)*dot(coefii,xpq1)*dot(coefjj,xpq1);

								//[T]*<S> = 0.5*[T1-T2]*(S1+S2) = S1*(T1-T2) + S2*(T1-T2) -> S<->i, T<->j
								double Hpq = 0.5*(kappa[pp] + kappa[qq]);
								val += sign*0.5*(area/3.)*ks*Hpq*dot(coefii,xpq1)*dot(coefjj,xpq1);
							}
						}
						
						val += 0.25*ks*epsilon*area*(dot(gradii, gradjj)  - dot(gradii, normal)*dot(gradjj, normal) );

						if ( (iibelongDBC == 0) && (jjbelongDBC == 0) ){
							MatSetValues(A,1,&ii_global,1,&jj_global,&val,ADD_VALUES);
						}
						if ( (iibelongDBC == 0) && (jjbelongDBC==1) ){ //if j is Dirichlet
							val=-val*TD(Th.vertices[jj][0],Th.vertices[jj][1],Th.vertices[jj][2]); //T2=0
							VecSetValues(b,1,&ii_global,&val,ADD_VALUES);
						}
						
						//Exact for P2
						val = 0.;
						for (int p=0; p<D; p++){
							int pp = Th.triangles[t][p];
							vector<double> xp1 = {Th.vertices[pp][0], Th.vertices[pp][1], Th.vertices[pp][2],1.};

							for (int q=0; q<p; q++){
								int qq = Th.triangles[t][q];
								vector<double> xq1 = {Th.vertices[qq][0], Th.vertices[qq][1], Th.vertices[qq][2],1.};
								vector<double> xpq1 = {0.5*xq1[0] + 0.5*xp1[0],0.5*xq1[1] + 0.5*xp1[1],0.5*xq1[2] + 0.5*xp1[2],1.};
			
								//ks*[T]*[S] = T1*S1 + T2*S2 - T2*S1 - T1*S2 
								val -= (ks/epsilon)*(area/3.)*dot(coefii,xpq1)*dot(coefjj_tilde,xpq1);

								//[T]*<S> = 0.5*[T1-T2]*(S1+S2) , T= SUM_j
								double Hpq = 0.5*(kappa[pp] + kappa[qq]);
								val -= sign*0.5*(area/3.)*ks*Hpq*dot(coefii,xpq1)*dot(coefjj_tilde,xpq1);
							}
						}
						
						val += 0.25*ks*epsilon*area*(dot(gradii, gradjj_tilde)  - dot(gradii, normal)*dot(gradjj_tilde, normal) );

						if ( (iibelongDBC == 0) && (jjbelongDBC_tilde == 0) ){
							MatSetValues(A,1,&ii_global,1,&jj_global_tilde,&val,ADD_VALUES);
						}
						if ( (iibelongDBC == 0) && (jjbelongDBC_tilde==1) ){ //if j is Dirichlet
							val=-val*TD(Th.vertices[jj][0],Th.vertices[jj][1],Th.vertices[jj][2]); //T2=0
							VecSetValues(b,1,&ii_global,&val,ADD_VALUES);
						}
						

					}
					

					
				}

				break;
		}
        
    }

}
