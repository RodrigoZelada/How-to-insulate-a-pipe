void nitscheDiscontinuousVentcell3D(Mat A, Vec b, mesh3D &Th, vector<double> &kappa, vector<vector<double>> &u, map<int,int> Gamma, mesh3D &Th2, double (*f)(double, double, double), double (*TD)(double, double, double), DirichletBC &DBC1, DirichletBC &DBC2)
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
	int iibelongDBC,jjbelongDBC;

	double ae,be,ce,length,eta,kmin,gamma;
	kmin = min(k1,k2);
	gamma = lambda*kmin;

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
			
			case 100: //label_interface
				if (region == region1){ 
					normal[0]=normalOmegai[0];normal[1]=normalOmegai[1];normal[2]=normalOmegai[2];
					sign = 1.; kRegion=k1; 
				}
				else{
					normal[0]=-normalOmegai[0];normal[1]=-normalOmegai[1];normal[2]=-normalOmegai[2];
					sign = -1.; kRegion=k2; 
				}
				vector<double> v0={Th.vertices[Th.triangles[t][0]][0],Th.vertices[Th.triangles[t][0]][1],Th.vertices[Th.triangles[t][0]][2]};
				vector<double> v1={Th.vertices[Th.triangles[t][1]][0],Th.vertices[Th.triangles[t][1]][1],Th.vertices[Th.triangles[t][1]][2]};
				vector<double> v2={Th.vertices[Th.triangles[t][2]][0],Th.vertices[Th.triangles[t][2]][1],Th.vertices[Th.triangles[t][2]][2]};
				vector<double> v01 = {v0[0]-v1[0],v0[1]-v1[1],v0[2]-v1[2]};
				vector<double> v02 = {v0[0]-v2[0],v0[1]-v2[1],v0[2]-v2[2]};
				vector<double> v21 = {v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]};
				length = max(max(norm(v01),norm(v02)),norm(v21));
				eta = epsilon/ks;
				ae = 1./(eta + gamma*length); be = (gamma*length*eta)/(eta+gamma*length); ce = (gamma*length)/(eta+gamma*length);
				
				int t_tilde = Gamma.at(t);
				n_tilde = Th2.triangleToTetra[t_tilde];

				for (i=0; i<=D; i++){ //antes era i<D
					ii = Th.tetrahedra[n][i]; itet=i;

					for (int k=0;k<=D;k++){
						coefii[k] = coef[n][4*itet+k];
						if (k<D){gradii[k] = coefii[k];}
					}

					if (region == region1){ 
						iibelongDBC = DBC1.vertices[ii][1]; 
						ii_global = DBC1.vertices[ii][0]; 
					}
					else{
						iibelongDBC = DBC2.vertices[ii][1]; 
						ii_global = DBC2.vertices[ii][0]+Nv1; 
					}

					for (j=0; j<=D; j++){
						jj = Th.tetrahedra[n][j]; jtet=j;

						if (region == region1){ 
							jjbelongDBC = DBC1.vertices[jj][1]; jj_global = DBC1.vertices[jj][0];
						}
						else{
							jjbelongDBC = DBC2.vertices[jj][1];jj_global = DBC2.vertices[jj][0]+Nv1;
						}

						for (int k=0;k<=D;k++){
							coefjj[k] = coef[n][4*jtet+k];
							if (k<D){gradjj[k] = coefjj[k]; }
						}

						val = 0.;
						val += 0.25*ks*epsilon*area*(dot(gradii, gradjj)  - dot(gradii, normal)*dot(gradjj, normal) );
						//- int_F b * <k grad(u)*n> <k grad (v) *n> = - int_F 0.25*b*(k1 grad(u1)*n + k2 grad(u2)*n)*(k1 grad(v1)*n + k2 grad(v2)*n) 
						val -= 0.25*be*area*(kRegion*dot(gradii, normal)*kRegion*dot(gradjj,normal) );

						//Exact for P2
						for (int p=0; p<D; p++){
							int pp = Th.triangles[t][p];
							vector<double> xp1 = {Th.vertices[pp][0], Th.vertices[pp][1], Th.vertices[pp][2],1.};

							// -int1d(c*(<kgradN(Th)>[Sh] + <kgradN(Sh)>[Th]))
							//val -= sign*(0.5*ce*area*(kRegion*dot(gradjj,normal) *  dot(coefii,xp1) ) + 0.5*ce*area*(kRegion*dot(gradii,normal) *  dot(coefjj,xp1) ));
							//-int_F ce*[kgradN(T)*S] = -int_F ce*(S1*k1*gradN(T1)-S2*k2*gradN(T2))
							val -= sign*ce*area*kRegion*dot(gradjj,normal) *  dot(coefii,xp1);

							for (int q=0; q<p; q++){
								int qq = Th.triangles[t][q];
								vector<double> xq1 = {Th.vertices[qq][0], Th.vertices[qq][1], Th.vertices[qq][2],1.};
								vector<double> xpq1 = {0.5*xq1[0] + 0.5*xp1[0],0.5*xq1[1] + 0.5*xp1[1],0.5*xq1[2] + 0.5*xp1[2],1.};
			
								//ks*[T]*[S] = T1*S1 + T2*S2 - T2*S1 - T1*S2 
								val += ae*(area/3.)*dot(coefii,xpq1)*dot(coefjj,xpq1);

								//[T]*<S> = 0.5*[T1-T2]*(S1+S2) = S1*(T1-T2) + S2*(T1-T2) -> S<->i, T<->j
								double Hpq = 0.5*(kappa[pp] + kappa[qq]);
								val += sign*0.5*(area/3.)*ks*Hpq*dot(coefii,xpq1)*dot(coefjj,xpq1);
							}
						}
						
						
						if ( (iibelongDBC == 0) && (jjbelongDBC == 0) ){
							MatSetValues(A,1,&ii_global,1,&jj_global,&val,ADD_VALUES);
						}
						if ( (iibelongDBC == 0) && (jjbelongDBC==1) ){ //if j is Dirichlet
							val=-val*T1(Th.vertices[jj][0],Th.vertices[jj][1],Th.vertices[jj][2]); //T2=0
							VecSetValues(b,1,&ii_global,&val,ADD_VALUES);
						}
						

					}

					//jj_tilde
					for (j_tilde=0;j_tilde<=D;j_tilde++){
						jj_tilde = Th2.tetrahedra[n_tilde][j_tilde];

						if (region == region1){ 
							jjbelongDBC = DBC2.vertices[jj_tilde][1]; jj_global_tilde = DBC2.vertices[jj_tilde][0] + Nv1;
						}
						else{
							jjbelongDBC = DBC1.vertices[jj_tilde][1]; jj_global_tilde =  DBC1.vertices[jj_tilde][0];
						}

						for (int l=0;l<=D;l++){
							coefjj_tilde[l] = coef2[n_tilde][4*j_tilde+l];
							if (l<Th.d){gradjj_tilde[l] = coefjj_tilde[l];}
						}

						//Exact for P2
						val = 0.;
						val += 0.25*ks*epsilon*area*(dot(gradii, gradjj_tilde)  - dot(gradii, normal)*dot(gradjj_tilde, normal) );
						//- int_F b * <k grad(u)*n> <k grad (v) *n> = - int_F 0.25*b*(k1 grad(u1)*n + k2 grad(u2)*n)*(k1 grad(v1)*n + k2 grad(v2)*n) 
						val -= 0.25*be*area*(k1*dot(gradii, normal)*k2*dot(gradjj_tilde,normal) );

						for (int p=0; p<D; p++){
							int pp = Th.triangles[t][p];
							vector<double> xp1 = {Th.vertices[pp][0], Th.vertices[pp][1], Th.vertices[pp][2],1.};

							// -int1d(c*(<kgradN(Th)>[Sh] + <kgradN(Sh)>[Th]))
							//val -= sign*(0.5*ce*area*(kOtherRegion*dot(gradjj_tilde,normal) *  dot(coefii,xp1) ) - 0.5*ce*area*(kRegion*dot(gradii,normal) *  dot(coefjj_tilde,xp1) ));

							for (int q=0; q<p; q++){
								int qq = Th.triangles[t][q];
								vector<double> xq1 = {Th.vertices[qq][0], Th.vertices[qq][1], Th.vertices[qq][2],1.};
								vector<double> xpq1 = {0.5*xq1[0] + 0.5*xp1[0],0.5*xq1[1] + 0.5*xp1[1],0.5*xq1[2] + 0.5*xp1[2],1.};
			
								//ks*[T]*[S] = T1*S1 + T2*S2 - T2*S1 - T1*S2 
								val -= ae*(area/3.)*dot(coefii,xpq1)*dot(coefjj_tilde,xpq1);

								//[T]*<S> = 0.5*[T1-T2]*(S1+S2) , T= SUM_j
								double Hpq = 0.5*(kappa[pp] + kappa[qq]);
								val -= sign*0.5*(area/3.)*ks*Hpq*dot(coefii,xpq1)*dot(coefjj_tilde,xpq1);
							}
						}

						
						if ( (iibelongDBC == 0) && (jjbelongDBC == 0) ){
							MatSetValues(A,1,&ii_global,1,&jj_global_tilde,&val,ADD_VALUES);
						}
						if ( (iibelongDBC == 0) && (jjbelongDBC==1) ){ //if j is Dirichlet
							val=-val*T1(Th.vertices[jj][0],Th.vertices[jj][1],Th.vertices[jj][2]); //T2=0
							VecSetValues(b,1,&ii_global,&val,ADD_VALUES);
						}
					
				}

					
				}

				break;
		
        }
    }

}
