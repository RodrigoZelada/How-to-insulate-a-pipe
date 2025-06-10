//g++ MatrixBinaryPetsc.cpp -o MatrixBinaryPetsc -std=c++14 -larmadillo
//./MatrixBinaryPetsc output/case3Vol
//time ./MatrixBinaryPetsc output/case3Vol -ksp_view_mat ascii:A.txt -mat_view draw -draw_save
//time mpiexec -n 8 ./MatrixBinaryPetsc output/case3Vol -pc_type lu -pc_factor_mat_solver_type superlu_dist -ksp_monitor
//time mpiexec -n 8 ./MatrixBinaryPetsc output/case3Vol -ksp_type gmres -pc_type bjacobi -sub_pc_type ilu -ksp_monitor
//time mpiexec -n 8 ./MatrixBinaryPetsc output/case3Vol -ksp_type gmres -pc_type asm -sub_pc_type ilu -ksp_monitor

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <time.h>
#include <string>
#include <vector>
#include <map>
#include <set>

#include <petsc.h>

using namespace std;

string output;
double tau;

const int region1 = 20; //fluid
const int region2 = 3; //isolant
const int label_interface = 100;
const int label_robin = 10;
const double k1=1.5; //0.15*1e-6; //thermal diffusivity D or alpha of water;
const double k2=1.; //0.1*1e-6; //thermal diffusivity D or alpha of rubber;
const double rhocp = 1e2; //1e6instead of 1e-6 in diffusivity

const double epsilon = 1e-2;
const double ks = 5.;//111. 111.*1e-6; //thermal diffusivity D or alpha of cupper;
const double alpha = 1e0;

double T1(double x, double y, double z) { return 40.;} //40.
double T2(double x, double y, double z) { return 40.;} //40.
double F1(double x, double y, double z) { return 0.;}
double F2(double x, double y, double z) { return 0.;}
double Text(double x, double y, double z) { return 0.;}
double T00(double x, double y, double z) { return 0.;}

double nx(double x,double y, double z){return 0.;}
double ny(double x,double y, double z){return (y-0.5)/(sqrt(pow(y-0.5,2.)+pow(z-0.5,2.))+1e-5);}
double nz(double x,double y, double z){return (z-0.5)/(sqrt(pow(y-0.5,2.)+pow(z-0.5,2.))+1e-5);}

double thetax(double x,double y, double z){return -tau*nx(x,y,z);}
double thetay(double x,double y, double z){return -tau*ny(x,y,z);}
double thetaz(double x,double y, double z){return -tau*nz(x,y,z);}

double dxthetax(double x,double y, double z){return 0.;}
double dythetax(double x,double y, double z){return 0.;}
double dzthetax(double x,double y, double z){return 0.;}

double dxthetay(double x,double y, double z){return 0.;}
double dythetay(double x,double y, double z){return -tau*pow(z-0.5,2.)/(pow(pow(y-0.5,2) + pow(z-0.5,2),1.5)+1e-5);}
double dzthetay(double x,double y, double z){return tau*((y-0.5)*(z-0.5))/(pow(pow(y-0.5,2.) + pow(z-0.5,2.),1.5)+1e-5);}

double dxthetaz(double x,double y, double z){return 0.;}
double dythetaz(double x,double y, double z){return tau*((y-0.5)*(z-0.5))/(pow(pow(y-0.5,2) + pow(z-0.5,2.),1.5)+1e-5);}
double dzthetaz(double x,double y, double z){return -tau*pow(y-0.5,2.)/(pow(pow(y-0.5,2) + pow(z-0.5,2),1.5)+1e-5);}

#include "linalg.h"
#include "mesh3D.h"
#include "trunc.h"
#include "medit.h"
#include "fespace.h"
#include "data.h"
#include "newIndicesWithoutBoundaryConditions.h"
#include "laplaceDiscontinuousVentcell3DTP2.h"
#include "adjointDiscontinuousVentcell3DExt.h"

int main(int argc, char ** argv) {
	output = argv[1];
	tau = atof(argv[2]);
	cout << "tau = " << tau << ", thetay(1,1,1) = " << thetay(1,1,1) << endl;
	/*Duplicate degrees at interface*/
	mesh3D Th(output+"/Th.mesh");
	map<int,int> mapGlobalTetrahedra1;
	map<int,int> mapGlobalTriangles1;
	map<int,int> mapGlobalVertices1;
	map<int,int> mapGlobalTetrahedra2;
	map<int,int> mapGlobalTriangles2;
	map<int,int> mapGlobalVertices2;
	vector<int> localToGlobalVertices1;
	vector<int> localToGlobalVertices2;
	
	mesh3D Th1 = trunc(Th,region1,mapGlobalTetrahedra1,mapGlobalTriangles1,mapGlobalVertices1,localToGlobalVertices1);
	mesh3D Th2 = trunc(Th,region2,mapGlobalTetrahedra2,mapGlobalTriangles2,mapGlobalVertices2,localToGlobalVertices2);

	int lDBC1[1] = {1};
	int nDBC1 = 1;
	DirichletBC DBC1(Th1,lDBC1,nDBC1);
	DirichletBC DBC2(Th2);

	map<int, int> Gamma1;
	map<int, int> Gamma2;

	int t1,t2;
	for (int t_global=0; t_global<Th.nt; t_global++){
		if (Th.triangles[t_global][Th.d]==label_interface){
			t1 = mapGlobalTriangles1.at(t_global);
			t2 = mapGlobalTriangles2.at(t_global);
			Gamma1.insert({t1, t2});
			Gamma2.insert({t2, t1});
		}
	}

	vector<double> kappa = readSolution(output+"/kappa.sol");
	vector<double>  kappa1 = truncP1(kappa,mapGlobalVertices1);
	vector<double>  kappa2 = truncP1(kappa,mapGlobalVertices2);
	vector<vector<double>>  u = readSolutionMat(output+"/ux.sol",output+"/uy.sol",output+"/uz.sol");

	int Nv1 = DBC1.nvNotDBC;
	int Nv2 = DBC2.nvNotDBC;
	int N = Nv1+Nv2;

	Mat A;
    Vec x,b;
    KSP ksp;
	//PetscInt Istart,Iend;
	PetscCall(PetscInitialize(&argc,&argv,NULL,"Solve laplacian\n"));

    PetscCall(VecCreate(PETSC_COMM_WORLD,&b));
    PetscCall(VecSetSizes(b,PETSC_DECIDE,N));
    PetscCall(VecSetFromOptions(b));
    //PetscCall(VecGetOwnershipRange(b,&Istart,&Iend));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
    PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
    PetscCall(MatSetFromOptions(A));
	PetscCall(MatSetUp(A)); //PetscCall(MatSeqAIJSetPreallocation(A, 30, NULL)); 
    //PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
	
	laplaceDiscontinuousVentcell3DTP2(A, b, Th1, kappa1, u, Gamma1, Th2, &Text, &T1, DBC1, DBC2);
	laplaceDiscontinuousVentcell3DTP2(A, b, Th2, kappa2, u, Gamma2, Th1, &Text, &T1, DBC1, DBC2);

	PetscCall(VecAssemblyBegin(b));
    PetscCall(VecAssemblyEnd(b));
    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

	PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
    PetscCall(KSPSetOperators(ksp,A,A));
    PetscCall(KSPSetFromOptions(ksp));

	PetscCall(VecDuplicate(b,&x));
    PetscCall(VecSet(x,0.0));
	PetscCall(KSPSolve(ksp,b,x));
    PetscCall(VecViewFromOptions(x,NULL,"-vec_view"));

	/*Save solution*/

	PetscScalar *T;
    PetscCall(VecGetArray(x,&T));
	PetscMPIInt    rank;
	PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
	vector<double> uh1 = P1Function(Th1, &T1); vector<double> uh2 = P1Function(Th2, &T2);
	if (rank==0){
		for (int i=0;i<Th1.nv;i++){
			if (DBC1.vertices[i][0] >= 0){
				uh1[i] = T[DBC1.vertices[i][0]];
			}
		}
		for (int i=0;i<Th2.nv;i++){
			if (DBC2.vertices[i][0] >= 0){
				uh2[i] = T[Nv1+DBC2.vertices[i][0]]; 
			}
		}
		savemesh(Th2,output+"/Th2.mesh");
		saveSolution(output+"/T1.sol",uh1);saveSolution(output+"/T2.sol",uh2);
	}

	PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&A));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&b));

	/*Compute J*/

	vector<double> T0 = P1Function(Th2,&Text);
	vector<double> diffTExt(Th2.nv); for (int i=0; i<Th2.nv;i++){diffTExt[i] = uh2[i]-T0[i];}

	vector<double> normalOmegai(Th2.d);
	vector<vector<double>> coef2 = baseSpace(Th2);
	
	vector<vector<double>> Gradtheta(Th2.d,vector<double>(Th2.d)); 
	vector<vector<double>> IGradtheta(Th2.d,vector<double>(Th2.d));  
	
	vector<vector<double>> NormalP1 = P1Vector(Th2, &nx, &ny, &nz);
	vector<vector<double>> GradThetaxP1 = P1Vector(Th2, &dxthetax, &dythetax, &dzthetax);
	vector<vector<double>> GradThetayP1 = P1Vector(Th2, &dxthetay, &dythetay, &dzthetay);
	vector<vector<double>> GradThetazP1 = P1Vector(Th2, &dxthetaz, &dythetaz, &dzthetaz);

	double IGT[3][3];

	double J=0.;
    for (int t=0; t<Th2.nt; t++){

        int label_triangle = Th2.triangles[t][Th2.d];
		double area = Th2.area[t];
		if (label_triangle == label_robin){	
			for (int p=0; p<Th2.d; p++){
				int pp = Th2.triangles[t][p];
				vector<double> xp1 = {Th2.vertices[pp][0], Th2.vertices[pp][1], Th2.vertices[pp][2],1.};

				for (int q=0; q<p; q++){
					int qq=Th2.triangles[t][q];
					vector<double> xq1 = {Th2.vertices[qq][0], Th2.vertices[qq][1], Th2.vertices[qq][2],1.};
					vector<double> xpq1 = {0.5*xq1[0] + 0.5*xp1[0],0.5*xq1[1] + 0.5*xp1[1],0.5*xq1[2] + 0.5*xp1[2],1.};
					
					double val = 0.5*diffTExt[pp] + 0.5*diffTExt[qq];
					normalOmegai[0]=0.5*NormalP1[pp][0]+0.5*NormalP1[pp][0];normalOmegai[1]=0.5*NormalP1[pp][1]+0.5*NormalP1[qq][1];normalOmegai[2]=0.5*NormalP1[pp][2]+0.5*NormalP1[qq][2];

					//(I+grad Theta)^⁻t = (I-grad Theta^t)
					IGradtheta[0][0] = 1.-0.5*GradThetaxP1[pp][0]-0.5*GradThetaxP1[qq][0];   IGradtheta[0][1] =   -0.5*GradThetayP1[pp][0]-0.5*GradThetayP1[qq][0]; IGradtheta[0][2] =   -0.5*GradThetazP1[pp][0]-0.5*GradThetazP1[qq][0];
					IGradtheta[1][0] =   -0.5*GradThetaxP1[pp][1]-0.5*GradThetaxP1[qq][1];   IGradtheta[1][1] = 1.-0.5*GradThetayP1[pp][1]-0.5*GradThetayP1[qq][1]; IGradtheta[1][2] =   -0.5*GradThetazP1[pp][1]-0.5*GradThetazP1[qq][1];
					IGradtheta[2][0] =   -0.5*GradThetaxP1[pp][2]-0.5*GradThetaxP1[qq][2];   IGradtheta[2][1] =   -0.5*GradThetayP1[pp][2]-0.5*GradThetayP1[qq][2]; IGradtheta[2][2] = 1.-0.5*GradThetazP1[pp][2]-0.5*GradThetazP1[qq][2];

					//det(I+grad Theta)
					IGT[0][0] = 1.+0.5*GradThetaxP1[pp][0]+0.5*GradThetaxP1[qq][0]; IGT[0][1] =    0.5*GradThetaxP1[pp][1]+0.5*GradThetaxP1[qq][1]; IGT[0][2] =    0.5*GradThetaxP1[pp][2]+0.5*GradThetaxP1[qq][2];
					IGT[1][0] =    0.5*GradThetayP1[pp][0]+0.5*GradThetayP1[qq][0]; IGT[1][1] = 1.+0.5*GradThetayP1[pp][1]+0.5*GradThetayP1[qq][1]; IGT[1][2] =    0.5*GradThetayP1[pp][2]+0.5*GradThetayP1[qq][2];
					IGT[2][0] =    0.5*GradThetazP1[pp][0]+0.5*GradThetazP1[qq][0]; IGT[2][1] =    0.5*GradThetazP1[pp][1]+0.5*GradThetazP1[qq][1]; IGT[2][2] = 1.+0.5*GradThetazP1[pp][2]+0.5*GradThetazP1[qq][2];

					J += (area/3.)*pow(alpha,2.)*pow(val,2.)*det3x3(IGT)*norm(MatrixByVector(IGradtheta,normalOmegai));
				}
			}
		}
	}
	cout << "Text = " << Text(0.,0.,0.) << ", int2d = " << J << endl;
	ofstream file(output+"/J.gp");
	file.precision(8);
	file << J << endl;
	file.close();

	/*Adjoint*/
	Mat AA;
    Vec xA,bA;
    KSP kspA;

    PetscCall(VecCreate(PETSC_COMM_WORLD,&bA));
    PetscCall(VecSetSizes(bA,PETSC_DECIDE,N));
    PetscCall(VecSetFromOptions(bA));
    //PetscCall(VecGetOwnershipRange(b,&Istart,&Iend));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&AA));
    PetscCall(MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,N,N));
    PetscCall(MatSetFromOptions(AA));
	PetscCall(MatSetUp(AA)); //PetscCall(MatSeqAIJSetPreallocation(A, 30, NULL)); 
	
	vector<double> G(Th2.nv); for (int i=0;i<Th2.nv;i++){G[i] = 2.*alpha*alpha*diffTExt[i];}
	adjointDiscontinuousVentcell3D(AA, bA, Th1, G, kappa1, u, Gamma1, Th2, DBC1, DBC2);
	adjointDiscontinuousVentcell3D(AA, bA, Th2, G, kappa2, u, Gamma2, Th1, DBC1, DBC2);

	PetscCall(VecAssemblyBegin(bA));
    PetscCall(VecAssemblyEnd(bA));
    PetscCall(MatAssemblyBegin(AA,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(AA,MAT_FINAL_ASSEMBLY));

	PetscCall(KSPCreate(PETSC_COMM_WORLD,&kspA));
    PetscCall(KSPSetOperators(kspA,AA,AA));
    PetscCall(KSPSetFromOptions(kspA));

	PetscCall(VecDuplicate(bA,&xA));
    PetscCall(VecSet(xA,0.0));
	PetscCall(KSPSolve(kspA,bA,xA));
    PetscCall(VecViewFromOptions(xA,NULL,"-vec_view"));

	PetscScalar *TA;
    PetscCall(VecGetArray(xA,&TA));
	//I use T00 to define vh1 = 0.;
	vector<double> vh1 = P1Function(Th1, &T00); vector<double> vh2 = P1Function(Th2, &T00);
	if (rank==0){
		for (int i=0;i<Th1.nv;i++){
			if (DBC1.vertices[i][0] >= 0){
				vh1[i] = TA[DBC1.vertices[i][0]];
			}
		}
		for (int i=0;i<Th2.nv;i++){
			if (DBC2.vertices[i][0] >= 0){
				vh2[i] = TA[Nv1+DBC2.vertices[i][0]]; 
			}
		}
		saveSolution(output+"/R1.sol",vh1);saveSolution(output+"/R2.sol",vh2);
	}

	PetscCall(KSPDestroy(&kspA));
    PetscCall(MatDestroy(&AA));
    PetscCall(VecDestroy(&xA));
    PetscCall(VecDestroy(&bA));

    PetscCall(PetscFinalize());
	
}
