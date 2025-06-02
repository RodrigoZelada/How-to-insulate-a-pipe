//define P1 given functions
vector<double> P1Function(mesh3D &Th, double (*f)(double, double, double)){
  vector<double> TExact(Th.nv);
  for (int i=0; i < Th.nv; i++){
    TExact[i] = f(Th.vertices[i][0],Th.vertices[i][1],Th.vertices[i][2]);
  }
  return TExact;
}

vector<vector<double>> GradP1Function(mesh3D &Th, vector<double> &thetaP1, vector<vector<double>> &coef){
  vector<vector<double>> GradthetaP1(Th.ntet,vector<double>(Th.d));

  for (int n = 0; n < Th.ntet; n++){ //gradient is constant by tetra
		for (int k=0; k<Th.d; k++){ //dimensions x,y,z
			GradthetaP1[n][k] = 0.;
		}
		for (int k=0;k<Th.d;k++){ //dimensions x,y,z
			for (int i=0; i<=Th.d; i++){ //each vertex of the tetrahedron
				int ii = Th.tetrahedra[n][i];
        GradthetaP1[n][k] += thetaP1[ii]*coef[n][4*i+k];
			}
    }
  }
    return GradthetaP1;
}

vector<vector<double>> P1Vector(mesh3D &Th, double (*f1)(double, double, double), double (*f2)(double, double, double),double (*f3)(double, double, double) ){
  vector<vector<double>> TExact(Th.nv,vector<double>(Th.d));
  for (int i=0; i < Th.nv; i++){
    TExact[i][0] = f1(Th.vertices[i][0],Th.vertices[i][1],Th.vertices[i][2]);
    TExact[i][1] = f2(Th.vertices[i][0],Th.vertices[i][1],Th.vertices[i][2]);
    TExact[i][2] = f3(Th.vertices[i][0],Th.vertices[i][1],Th.vertices[i][2]);
  }
  return TExact;
}

vector<vector<double>> GradP1Vector(mesh3D &Th, vector<vector<double>>  &thetaP1, vector<vector<double>> &coef){
  vector<vector<double>> GradthetaP1(Th.ntet,vector<double>(Th.d*(Th.d+1)));

  for (int n = 0; n < Th.ntet; n++){ //tetrahedra
		//gradient
		for (int k=0; k<Th.d; k++){ //dimensions x,y,z
      for (int l=0 ; l<Th.d; l++){ //each function thetax, thetay, thetaz
			  GradthetaP1[n][Th.d*l+k] = 0.;
      }
		}
		for (int k=0;k<Th.d;k++){  //dimensions x,y,z
      for (int l=0; l<Th.d; l++){  //each function thetax, thetay, thetaz
        for (int i=0; i<=Th.d; i++){ //each vertex of the tetrahedron
          int ii = Th.tetrahedra[n][i];
          GradthetaP1[n][Th.d*l+k] += thetaP1[ii][l]*coef[n][4*i+k]; 
        }
      }

    }
  }
    return GradthetaP1;
}