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

const int region1 = 20; //fluid
const int region2 = 3; //isolant
const int label_interface = 100;
const int label_robin = 10;
const double k1=0.15; //0.15*1e-6; //thermal diffusivity D or alpha of water;
const double k2=0.1; //0.1*1e-6; //thermal diffusivity D or alpha of rubber;
const double lambda=1e-2;
const double rhocp = 1e6; //1e6instead of 1e-6 in diffusivity

const double epsilon = 1e-3;
const double ks = 111.;//111. 111.*1e-6; //thermal diffusivity D or alpha of cupper;
const double alpha = 2e1;

double T1(double x, double y, double z) { return 40.;} //40.
double T2(double x, double y, double z) { return 40.;} //40.
double F1(double x, double y, double z) { return 0.;}
double F2(double x, double y, double z) { return 0.;}
double Text(double x, double y, double z) { return 0.;}
double T00(double x, double y, double z) { return 0.;}

#include "linalg.h"
#include "mesh3D.h"
#include "trunc.h"
#include "medit.h"
#include "fespace.h"
#include "data.h"
#include "newIndicesWithoutBoundaryConditions.h"
//#include "laplaceDiscontinuousVentcell3D.h"
//#include "adjointDiscontinuousVentcell3DExt.h"
#include "nitscheDiscontinuousVentcell3D.h"
#include "nitscheAdjointDiscontinuousVentcell3DExt.h"

int main(int argc, char ** argv) {
	output = argv[1];
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

	//Bijection between triangles at the interface: t <-> t_tilde //
	//index triangle in Th1 || index triangle in Th2

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
	
	nitscheDiscontinuousVentcell3D(A, b, Th1, kappa1, u, Gamma1, Th2, &Text, &T1, DBC1, DBC2);
	nitscheDiscontinuousVentcell3D(A, b, Th2, kappa2, u, Gamma2, Th1, &Text, &T1, DBC1, DBC2);

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
	double J=0.; double Q=0.;
    for (int t=0; t<Th2.nt; t++){
        int label_triangle = Th2.triangles[t][Th2.d];
		double area = Th2.area[t];
		if (label_triangle == label_robin){	
			for (int p=0; p<Th2.d; p++){
				int pp = Th2.triangles[t][p];
				Q += (area/3.)*alpha*diffTExt[pp];
				for (int q=0; q<p; q++){
					int qq=Th2.triangles[t][q];
					double val = 0.5*diffTExt[pp] + 0.5*diffTExt[qq];
					J += (area/3.)*alpha*alpha*val*val;
				}
			}
		}
	}

	ofstream file(output+"/J.gp");
	file.precision(8);
	file << J << endl;
	file.close();

	ofstream file2(output+"/Q.gp");
	file2.precision(8);
	file2 << Q << endl;
	file2.close();

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
	
	vector<double> G(Th2.nv); for (int i=0;i<Th2.nv;i++){G[i] = alpha;}
	nitscheAdjointDiscontinuousVentcell3D(AA, bA, Th1, G, kappa1, u, Gamma1, Th2, DBC1, DBC2);
	nitscheAdjointDiscontinuousVentcell3D(AA, bA, Th2, G, kappa2, u, Gamma2, Th1, DBC1, DBC2);

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
