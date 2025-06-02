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
const double rhocp = 1e6; //1e6instead of 1e-6 in diffusivity

const double epsilon = 1e-3;
const double ks = 111.;// 111.*1e-6; //thermal diffusivity D or alpha of cupper;
const double alpha = 2e1;
const double sigma1 = 0.3; //we write sigma but actually sigma^2, the variance
const double sigma2 = 0.7;

double T1(double x, double y, double z) { return 40.;} //40.
double T2(double x, double y, double z) { return 40.;} //40.

double Text1(double x, double y, double z) { return 20.*x;} //change this, Text^a = 0, Text^b = 20;
double Text2(double x, double y, double z) { return 10.*z;} //change this, Text^a = 0, Text^b = 20;
double T00(double x, double y, double z) { return 0.;} //changing variables T = T-Td -> Text = T-Td

#include "linalg.h"
#include "mesh3D.h"
#include "trunc.h"
#include "medit.h"
#include "fespace.h"
#include "data.h"
#include "newIndicesWithoutBoundaryConditions.h"
#include "laplaceDiscontinuousVentcell3D.h"
#include "adjointDiscontinuousVentcell3DExt.h"

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

	//only for Th2
	/*int count=0; //all the points that are in Th2, that belongs to a DBC in Th1
	for (int i=0; i<Th2.nv; i++){
		int i_global = localToGlobalVertices2[i]; 
		if ( (DBC2.vertices[i][0] >=0 ) && (mapGlobalVertices1.count(i_global) > 0)){ //Th1 and Th2 share the vertex i
			int i1 = mapGlobalVertices1.at(i_global);
			if (DBC1.vertices[i1][0] == -1){ //if is Dirichlet in Th1
				DBC2.vertices[i][0] = -1; DBC2.vertices[i][1] = DBC1.vertices[i1][1];
				count++;
			}
			else {
				DBC2.vertices[i][0] -= count;
			}
		}
		else{
			if (DBC2.vertices[i][0] >=0 ){ 
				DBC2.vertices[i][0] -= count;
			}
		}
	}
	DBC2.nvDBC += count; 
	DBC2.nvNotDBC -= count;*/

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

	/*Deterministic part*/
	Mat A;
    Vec x,b;
    KSP ksp;
	PetscCall(PetscInitialize(&argc,&argv,NULL,"Solve laplacian\n"));

    PetscCall(VecCreate(PETSC_COMM_WORLD,&b));
    PetscCall(VecSetSizes(b,PETSC_DECIDE,N));
    PetscCall(VecSetFromOptions(b));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
    PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
    PetscCall(MatSetFromOptions(A));
	PetscCall(MatSetUp(A)); //PetscCall(MatSeqAIJSetPreallocation(A, 30, NULL)); 

	laplaceDiscontinuousVentcell3D(A, b, Th1, kappa1, u, Gamma1, Th2, &T00, &T1, DBC1, DBC2);
	laplaceDiscontinuousVentcell3D(A, b, Th2, kappa2, u, Gamma2, Th1, &T00, &T1, DBC1, DBC2);

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

	/*First random variable*/
	Mat A1;
    Vec x1,b1;
    KSP ksp1;
	PetscCall(PetscInitialize(&argc,&argv,NULL,"Solve laplacian\n"));

    PetscCall(VecCreate(PETSC_COMM_WORLD,&b1));
    PetscCall(VecSetSizes(b1,PETSC_DECIDE,N));
    PetscCall(VecSetFromOptions(b1));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&A1));
    PetscCall(MatSetSizes(A1,PETSC_DECIDE,PETSC_DECIDE,N,N));
    PetscCall(MatSetFromOptions(A1));
	PetscCall(MatSetUp(A1)); //PetscCall(MatSeqAIJSetPreallocation(A, 30, NULL)); 

	laplaceDiscontinuousVentcell3D(A1, b1, Th1, kappa1, u, Gamma1, Th2, &Text1, &T00, DBC1, DBC2);
	laplaceDiscontinuousVentcell3D(A1, b1, Th2, kappa2, u, Gamma2, Th1, &Text1, &T00, DBC1, DBC2);

	PetscCall(VecAssemblyBegin(b1));
    PetscCall(VecAssemblyEnd(b1));
    PetscCall(MatAssemblyBegin(A1,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A1,MAT_FINAL_ASSEMBLY));

	PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp1));
    PetscCall(KSPSetOperators(ksp1,A1,A1));
    PetscCall(KSPSetFromOptions(ksp1));

	PetscCall(VecDuplicate(b1,&x1));
    PetscCall(VecSet(x1,0.0));
	PetscCall(KSPSolve(ksp1,b1,x1));

	/*Second random variable*/
	Mat A2;
    Vec x2,b2;
    KSP ksp2;
    PetscCall(VecCreate(PETSC_COMM_WORLD,&b2));
    PetscCall(VecSetSizes(b2,PETSC_DECIDE,N));
    PetscCall(VecSetFromOptions(b2));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&A2));
    PetscCall(MatSetSizes(A2,PETSC_DECIDE,PETSC_DECIDE,N,N));
    PetscCall(MatSetFromOptions(A2));
	PetscCall(MatSetUp(A2)); //PetscCall(MatSeqAIJSetPreallocation(A, 30, NULL)); 

	laplaceDiscontinuousVentcell3D(A2, b2, Th1, kappa1, u, Gamma1, Th2, &Text2, &T00, DBC1, DBC2);
	laplaceDiscontinuousVentcell3D(A2, b2, Th2, kappa2, u, Gamma2, Th1, &Text2, &T00, DBC1, DBC2);

	PetscCall(VecAssemblyBegin(b2));
    PetscCall(VecAssemblyEnd(b2));
    PetscCall(MatAssemblyBegin(A2,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A2,MAT_FINAL_ASSEMBLY));

	PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp2));
    PetscCall(KSPSetOperators(ksp2,A2,A2));
    PetscCall(KSPSetFromOptions(ksp2));

	PetscCall(VecDuplicate(b2,&x2));
    PetscCall(VecSet(x2,0.0));
	PetscCall(KSPSolve(ksp2,b2,x2));

	/*Save solution*/
	PetscScalar *T; PetscScalar *T11; PetscScalar *T22;
    PetscCall(VecGetArray(x,&T)); PetscCall(VecGetArray(x1,&T11)); PetscCall(VecGetArray(x2,&T22));
	vector<double> uh1 = P1Function(Th1, &T1); vector<double> uh2 = P1Function(Th2, &T2);
	vector<double> uh11 = P1Function(Th1, &T00); vector<double> uh21 = P1Function(Th2, &T00);
	vector<double> uh12 = P1Function(Th1, &T00); vector<double> uh22 = P1Function(Th2, &T00);
	PetscMPIInt    rank;
	PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
	if (rank==0){
		for (int i=0;i<Th1.nv;i++){
			if (DBC1.vertices[i][0] >= 0){
				uh1[i] = T[DBC1.vertices[i][0]]; uh11[i] = T11[DBC1.vertices[i][0]]; uh12[i] = T22[DBC1.vertices[i][0]];
			}
		}
		for (int i=0;i<Th2.nv;i++){
			if (DBC2.vertices[i][0] >= 0){
				uh2[i] = T[Nv1+DBC2.vertices[i][0]]; uh21[i] = T11[Nv1+DBC2.vertices[i][0]]; uh22[i] = T22[Nv1+DBC2.vertices[i][0]];
			}
		}
		savemesh(Th2,output+"/Th2.mesh");
		saveSolution(output+"/T1.sol",uh1);saveSolution(output+"/T2.sol",uh2);
		saveSolution(output+"/T11.sol",uh11);saveSolution(output+"/T21.sol",uh21);
		saveSolution(output+"/T12.sol",uh12);saveSolution(output+"/T22.sol",uh22);
	}

	PetscCall(KSPDestroy(&ksp)); PetscCall(KSPDestroy(&ksp1)); PetscCall(KSPDestroy(&ksp2));
    PetscCall(MatDestroy(&A)); PetscCall(MatDestroy(&A1)); PetscCall(MatDestroy(&A2));
    PetscCall(VecDestroy(&x)); PetscCall(VecDestroy(&x1)); PetscCall(VecDestroy(&x2));
    PetscCall(VecDestroy(&b)); PetscCall(VecDestroy(&b1)); PetscCall(VecDestroy(&b2));

	/*Compute J*/
	/*With both experiments, we can compute the expectation*/
	vector<double> T01 = P1Function(Th2,&Text1); vector<double> T02 = P1Function(Th2,&Text2);
	vector<double> diffTExt(Th2.nv); for (int i=0; i<Th2.nv;i++){diffTExt[i] = uh2[i];} //-0
	vector<double> diffTExt1(Th2.nv); for (int i=0; i<Th2.nv;i++){diffTExt1[i] = uh21[i]-T01[i];}
	vector<double> diffTExt2(Th2.nv); for (int i=0; i<Th2.nv;i++){diffTExt2[i] = uh22[i]-T02[i];}
	double J=0.;
    for (int t=0; t<Th2.nt; t++){
        int label_triangle = Th2.triangles[t][Th2.d];
		double area = Th2.area[t];
		if (label_triangle == label_robin){	
			for (int p=0; p<Th2.d; p++){
				int pp = Th2.triangles[t][p];
				for (int q=0; q<p; q++){
					int qq=Th2.triangles[t][q];
					double val = 0.5*diffTExt[pp] + 0.5*diffTExt[qq];  double val1 = 0.5*diffTExt1[pp] + 0.5*diffTExt1[qq]; double val2 = 0.5*diffTExt2[pp] + 0.5*diffTExt2[qq];
					J += (area/3.)*alpha*alpha*(val*val + val1*val1*sigma1 + val2*val2*sigma2);
				}
			}
		}
	}
	if (rank==0){
		ofstream file(output+"/J.gp");
		file.precision(8);
		file << J << endl;
		file.close();
	}

	/*Adjoint*/

	/*Text = 0 deterministic*/
	Mat AA;
    Vec xA,bA;
    KSP kspA;

    PetscCall(VecCreate(PETSC_COMM_WORLD,&bA));
    PetscCall(VecSetSizes(bA,PETSC_DECIDE,N));
    PetscCall(VecSetFromOptions(bA));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&AA));
    PetscCall(MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,N,N));
    PetscCall(MatSetFromOptions(AA));
	PetscCall(MatSetUp(AA)); 
	
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

	/*Text1*/
	Mat AA1;
    Vec xA1,bA1;
    KSP kspA1;

    PetscCall(VecCreate(PETSC_COMM_WORLD,&bA1));
    PetscCall(VecSetSizes(bA1,PETSC_DECIDE,N));
    PetscCall(VecSetFromOptions(bA1));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&AA1));
    PetscCall(MatSetSizes(AA1,PETSC_DECIDE,PETSC_DECIDE,N,N));
    PetscCall(MatSetFromOptions(AA1));
	PetscCall(MatSetUp(AA1)); 
	
	vector<double> G1(Th2.nv); for (int i=0;i<Th2.nv;i++){G1[i] = 2.*alpha*alpha*diffTExt1[i];}
	adjointDiscontinuousVentcell3D(AA1, bA1, Th1, G1, kappa1, u, Gamma1, Th2, DBC1, DBC2);
	adjointDiscontinuousVentcell3D(AA1, bA1, Th2, G1, kappa2, u, Gamma2, Th1, DBC1, DBC2);

	PetscCall(VecAssemblyBegin(bA1));
    PetscCall(VecAssemblyEnd(bA1));
    PetscCall(MatAssemblyBegin(AA1,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(AA1,MAT_FINAL_ASSEMBLY));

	PetscCall(KSPCreate(PETSC_COMM_WORLD,&kspA1));
    PetscCall(KSPSetOperators(kspA1,AA1,AA1));
    PetscCall(KSPSetFromOptions(kspA1));

	PetscCall(VecDuplicate(bA1,&xA1));
    PetscCall(VecSet(xA1,0.0));
	PetscCall(KSPSolve(kspA1,bA1,xA1));

	/*Text2*/
	Mat AA2;
    Vec xA2,bA2;
    KSP kspA2;

    PetscCall(VecCreate(PETSC_COMM_WORLD,&bA2));
    PetscCall(VecSetSizes(bA2,PETSC_DECIDE,N));
    PetscCall(VecSetFromOptions(bA2));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&AA2));
    PetscCall(MatSetSizes(AA2,PETSC_DECIDE,PETSC_DECIDE,N,N));
    PetscCall(MatSetFromOptions(AA2));
	PetscCall(MatSetUp(AA2)); 
	
	vector<double> G2(Th2.nv); for (int i=0;i<Th2.nv;i++){G2[i] = 2.*alpha*alpha*diffTExt2[i];}
	adjointDiscontinuousVentcell3D(AA2, bA2, Th1, G2, kappa1, u, Gamma1, Th2, DBC1, DBC2);
	adjointDiscontinuousVentcell3D(AA2, bA2, Th2, G2, kappa2, u, Gamma2, Th1, DBC1, DBC2);

	PetscCall(VecAssemblyBegin(bA2));
    PetscCall(VecAssemblyEnd(bA2));
    PetscCall(MatAssemblyBegin(AA2,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(AA2,MAT_FINAL_ASSEMBLY));

	PetscCall(KSPCreate(PETSC_COMM_WORLD,&kspA2));
    PetscCall(KSPSetOperators(kspA2,AA2,AA2));
    PetscCall(KSPSetFromOptions(kspA2));

	PetscCall(VecDuplicate(bA2,&xA2));
    PetscCall(VecSet(xA2,0.0));
	PetscCall(KSPSolve(kspA2,bA2,xA2));

	PetscScalar *TA; PetscScalar *TA1; PetscScalar *TA2;
    PetscCall(VecGetArray(xA,&TA)); PetscCall(VecGetArray(xA1,&TA1)); PetscCall(VecGetArray(xA2,&TA2));
	//I use T00 to define vh1 = 0.;
	vector<double> vh1 = P1Function(Th1, &T00); vector<double> vh2 = P1Function(Th2, &T00);
	vector<double> vh11 = P1Function(Th1, &T00); vector<double> vh21 = P1Function(Th2, &T00);
	vector<double> vh12 = P1Function(Th1, &T00); vector<double> vh22 = P1Function(Th2, &T00);
	if (rank==0){
		for (int i=0;i<Th1.nv;i++){
			if (DBC1.vertices[i][0] >= 0){
				vh1[i] = TA[DBC1.vertices[i][0]]; vh11[i] = TA1[DBC1.vertices[i][0]]; vh12[i] = TA2[DBC1.vertices[i][0]];
			}
		}
		for (int i=0;i<Th2.nv;i++){
			if (DBC2.vertices[i][0] >= 0){
				vh2[i] = TA[Nv1+DBC2.vertices[i][0]]; vh21[i] = TA1[Nv1+DBC2.vertices[i][0]]; vh22[i] = TA2[Nv1+DBC2.vertices[i][0]];
			}
		}
		saveSolution(output+"/R1.sol",vh1);saveSolution(output+"/R2.sol",vh2);
		saveSolution(output+"/R11.sol",vh11);saveSolution(output+"/R21.sol",vh21);
		saveSolution(output+"/R12.sol",vh12);saveSolution(output+"/R22.sol",vh22);
	}

	PetscCall(KSPDestroy(&kspA)); PetscCall(KSPDestroy(&kspA1)); PetscCall(KSPDestroy(&kspA2));
    PetscCall(MatDestroy(&AA)); PetscCall(MatDestroy(&AA1)); PetscCall(MatDestroy(&AA2));
    PetscCall(VecDestroy(&xA)); PetscCall(VecDestroy(&xA1)); PetscCall(VecDestroy(&xA2));
    PetscCall(VecDestroy(&bA)); PetscCall(VecDestroy(&bA1)); PetscCall(VecDestroy(&bA2));

    PetscCall(PetscFinalize());
	
}
