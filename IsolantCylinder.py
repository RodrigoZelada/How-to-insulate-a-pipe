from pyfreefem import FreeFemRunner, readFFArray
from nullspace_optimizer import Optimizable, nlspace_solve
from pymedit import Mesh, P1Function, trunc, Mesh3D, cube, mmg3d, generate3DMesh, P1Function3D, trunc3DMesh
from pymedit import saveToVtk, advect, P1Vector3D, mshdist
import numpy as np
import shutil
import pandas as pd

output='output/caseCylinder' #1
import os
os.makedirs(output,exist_ok=True)

config={'OUTPUT':output}

N=8 #number of cores
i = format(0,'04d')
meshes='meshes'
M0 = Mesh3D(meshes+"/CylinderOrdered.mesh")
Th0 = mmg3d(M0,extra_args='-nomove -noswap -nosurf -noinsert',nr=False)
#Th0=M0

epsB = 0.05 #NO CAMBIAR, ES EL QUE USO PARA NAVIER STOKES Y GEOMETRÍA!
eps = 0.02 #SE PUEDE CAMBIAR, ES EL EPSILON PARA LA BARRERA
r=0.1
R=0.2

barrier = lambda x,y,z: min(max((y-0.5)**2 + (z-0.5)**2 - R**2, x-eps), max((y-0.5)**2 + (z-0.5)**2 - R**2,1-eps-x) ) #(z-0.5)**2 + (y-0.5)**2 - 0.11**2
phi1 = lambda x,y,z: (z-0.5)**2 + (y-0.5)**2 - 0.1**2

cutIn = lambda x,y,z: min((y-0.5)**2 + (z-0.5)**2 - R**2, eps-x)
cutOut = lambda x,y,z: min((y-0.5)**2 + (z-0.5)**2 - R**2, x-1+eps)
cut = lambda x,y,z: max(cutIn(x,y,z),cutOut(x,y,z))

phi2 = lambda x,y,z: max((z-0.5)**2 + (y-0.5)**2 - 0.45**2,cut(x,y,z)) #(z-0.5)**2 + (y-0.5)**2 - 0.45**2 #phi02(x,y,z) #(z-0.5)**2 + (y-0.5)**2 - 0.45**2

# Meshing parameters
hmin = 1e-2
hmax = 5e-2
hgrad = 1.3 # Gradation factor
hausd = 5e-3 # Approximation quality factor, leave as it is

# Local parameter to prescribe local mesh size on the interface boundary (corresponding to the level set)
paramsmmg = f"""Parameters
5

10 Triangles 5e-3 1e-2 1e-3
1 Triangles 5e-3 1e-2 1e-3
2 Triangles 5e-3 1e-2 1e-3
3 Triangles 1e-2 2.5e-2 1e-3
4 Triangles 1e-2 2.5e-2 1e-3
"""

phiinit1=P1Function3D(Th0, phi1)

# Mesh with mmg2d
M1 = mmg3d(Th0,hmin,hmax,hgrad,hausd,ls=True,sol=phiinit1,params=paramsmmg, debug=1)
M1.Boundaries[100] = M1.Boundaries.pop(10)
M1.triangles[M1.triangles[:,3] == 10,3] = 100
M1.tetrahedra[M1.tetrahedra[:,4] == 3, 4] = 20

refermmg = f"""
                LSReferences
                2
                2 3 2
                20 nosplit

                Parameters
                6

                100 Triangles 1e-2 2e-2 1e-3
                10 Triangles 1e-2 2e-2 1e-3
                1 Triangles 5e-3 1e-2 1e-3
                2 Triangles 5e-3 1e-2 1e-3
                12 Triangles 5e-3 1e-2 1e-3
                13 Triangles 5e-3 1e-2 1e-3
"""

phiinit2=P1Function3D(M1,phi2)
Thinit = mmg3d(M1,hmin,hmax,hgrad,hausd,ls=True,sol=phiinit2,params=refermmg,extra_args="-rmc 1e-8",debug=1)
Thinit.save(output+'/Th00.mesh')

preamble="""
func int readSolFile(mesh & Th, string fileName, real[int] & phi){
    ifstream f(fileName);
    string dummy="";
    while(dummy(0:2)!="Sol"){
        f>>dummy;
    }
    int n;
    f >> n;
    if(n!=Th.nv){
        cout << "Error : the number of vertices in the file "+fileName+" does not correspond to the mesh in memory"<< endl;
        exit(1);
    }
    f >> dummy;
    f >> dummy;
    for(int i=0;i<Th.nv;i++){
        f>>phi[i];
    }
}

func int saveArray(string fileName, real[int] &value){
    ofstream file(fileName);
    file.precision(16);
    file << value;
}

func int readData(string fileName, real[int] &data){
    {
        ifstream f(fileName);
        f>>data;
    }
}

include "macros.edp"
load "medit"
load "msh3"
load "iovtk"

int region1 = 20; //interior circle
int region2 = 3; //exterior donuts
int empty = 2;
real k2=0.1; //0.1*1e-6; //thermal diffusivity D or alpha of rubber;
real vtarget=0.0942477796; 
real nu = 1e-2; //1e-3; //1e-6
real beta=2e1; //1e2
real corr=1.;
"""

mesh_code='''
load "PETSc"
macro dimension()3 //EOM
include "macro_ddm.idp"

mesh3 ThBox=readmesh3("$MESH");
savevtk("$OUTPUT/ThBox.vtu",ThBox);
mesh3 Th=trunc(ThBox, ((region == region1) || (region == region2)));
savemesh(Th,"$OUTPUT/Th.mesh");
{
    ofstream f("$OUTPUT/V.gp");
    f.precision(16);
    f << int3d(Th,region2,qforder=1)(1.) - vtarget  << endl;
}
'''

solve_kappa = """
mesh3 ThGlobal = Th;

int[int] n2o;
macro ThN2O() n2o // this tells buildDmesh to keep the local to global correspondence
buildDmesh(Th)

fespace Gh(Th, [P1,P1,P1]);
fespace GhGlobal(ThGlobal, [P1,P1,P1]);
fespace Ph(Th, P1);
fespace PhGlobal(ThGlobal, P1);
fespace PhBox(ThBox, P1);

Gh [gx, gy, gz];
GhGlobal [gxGlobal, gyGlobal, gzGlobal], [gxsum, gysum, gzsum];

Ph kappa, dOmega;
PhGlobal kappaGlobal, kappasum, dOmegaGlobal, nx, ny, nz;

PhBox dOmegaBox;
dOmegaBox[] = readsol("$OUTPUT/d1.sol");
dOmegaGlobal = dOmegaBox;

int[int] subIdxG  = restrict(Gh, GhGlobal, n2o);
int[int] subIdx  = restrict(Ph, PhGlobal, n2o);

dOmega[] = dOmegaGlobal[](subIdx);

varf probgsd([gx,gy,gz],[thetax, thetay,thetaz]) = int3d(Th,qforder=3)( tr([gx,gy,gz])*[thetax, thetay, thetaz]      )
                                                   +int3d(Th,qforder=2)( tr(grad(dOmega))*[thetax, thetay, thetaz]  );

Mat AG; //Global distributed matrix
{
macro def(i)[i, i#B, i#C] //EOM
macro init(i)[i, i, i] // EOM
createMat(Th, AG, [P1,P1,P1])
}

AG=probgsd(Gh,Gh,tgv=-1);
set(AG, sparams = "-ksp_max_it 99 -ksp_rtol 1e-6 -pc_type hypre"); //"-pc_type lu"
real[int] bG=probgsd(0,Gh,tgv=-1);
gx[] = AG^-1 * bG;

Mat A;
createMat(Th,A,P1)
varf probkappa(kappa, psi) = int3d(Th,qforder=3)(9*Th.hmax*Th.hmax*tr(grad(kappa))*grad(psi) + kappa*psi)
                             +int3d(Th,qforder=2)(div(gx,gy,gz)*psi);


A = probkappa(Ph, Ph, tgv=-1);
set(A, sparams = "-pc_type lu " );
real[int] b = probkappa(0, Ph, tgv=-1);
kappa[] = A^-1 * b;

gx[] .*= AG.D;
gxGlobal[](subIdxG) = gx[]; //local function from #0 interpolated by 0, on the global domain
mpiAllReduce(gxGlobal[], gxsum[], mpiCommWorld, mpiSUM); //problem because of global elements => need to use partition of unity

kappa[] .*= A.D;
kappaGlobal[](subIdx) = kappa[]; //local function from #0 interpolated by 0, on the global domain
mpiAllReduce(kappaGlobal[], kappasum[], mpiCommWorld, mpiSUM); //problem because of global elements => need to use partition of unity


if (mpirank==0){
    //save the solution
    savesol("$OUTPUT/kappa.sol",ThGlobal,kappasum);
    savesol("$OUTPUT/signeddistance.sol",ThGlobal,dOmegaGlobal);
}


"""

solve_ns = """
load "PETSc"
macro dimension()3 //EOM
include "macro_ddm.idp"

mesh3 Th=readmesh3("$OUTPUT/Th00.mesh");
mesh3 Th1 = trunc(Th, (region==region1));
mesh3 Th1Global = Th1;
savemesh(Th1Global,"$OUTPUT/ThFluid00.mesh");

int[int] n2o1;
macro Th1N2O() n2o1 // this tells buildDmesh to keep the local to global correspondence
buildDmesh(Th1)

fespace Vh1(Th1, [P2,P2,P2,P1]);
Vh1 [ux, uy, uz, p], [dux, duy, duz, dp];

varf stokes([ux, uy, uz, p], [wx, wy, wz, qq]) =
    int3d(Th1, qforder=3 )(2*nu*tr(EPS(ux,uy,uz))*EPS(wx,wy,wz) - p*div(wx,wy,wz) - qq*div(ux,uy,uz))
   -int3d(Th1, qforder=3)(p*qq*1e-8)
    +on(100,ux=0.,uy=0.,uz=0.)
   +on(1,ux=(0.1^2 - (y-0.5)^2 - (z-0.5)^2)/0.01, uy=0., uz=0.);

Mat Ans;
{
  macro def(i)[i, i#B, i#C, i#D] //
  macro init(i)[i, i, i, i] //
  createMat(Th1, Ans, [P2,P2,P2,P1])
}

Ans = stokes(Vh1, Vh1, tgv=-1);
set(Ans, sparams="-pc_type lu");

real[int] bns = stokes(0, Vh1, tgv=-1);
ux[] = Ans^-1 * bns;

Mat AOseen;
{
macro def(i)[i, i#B, i#C, i#D] //
macro init(i)[i, i, i, i] //
createMat(Th1, AOseen, [P2,P2,P2,P1]);
}
real[int] bOseen(Vh1.ndof);

real err=1.;
real tol = 1e-4;
int N=3;

varf Oseen([dux,duy,duz,dp],[vx,vy,vz,qq]) =
int3d(Th1)(2*nu*tr(EPS(dux,duy,duz))*EPS(vx,vy,vz)
+ tr(UgradV(dux,duy,duz,ux,uy,uz))*[vx,vy,vz]
+ tr(UgradV(ux,uy,uz,dux,duy,duz))*[vx,vy,vz]
- div(dux,duy,duz)*qq - div(vx,vy,vz)*dp
- 1e-8*dp*qq)
+int3d(Th1)(2*nu*tr(EPS(ux,uy,uz))*EPS(vx,vy,vz)
    + tr(UgradV(ux,uy,uz,ux,uy,uz))*[vx,vy,vz]
    - div(ux,uy,uz)*qq - div(vx,vy,vz)*p
    - 1e-8*p*qq)
+on(1,100, dux=0.,duy=0.,duz=0.);

int n=0;
while ( (err > tol) && (n<N)) {
  AOseen = Oseen(Vh1, Vh1, tgv=-1);
  set(AOseen, sparams = "-pc_type lu");
  bOseen = Oseen(0, Vh1, tgv=-1);

  dux[] = AOseen^-1 * bOseen;
  ux[] -= dux[];
  real errLoc = sqrt(int3d(Th1)(tr(GRAD(dux,duy,duz))*GRAD(dux,duy,duz) + tr([dux,duy,duz])*[dux,duy,duz]) / int3d(Th1)(tr(GRAD(ux,uy,uz))*GRAD(ux,uy,uz) + tr([ux,uy,uz])*[ux,uy,uz]));

  mpiReduce(errLoc, err, processor(0), mpiSUM);

  if (mpirank == 0){cout << "err = " << err << endl;}
  n++;
}
if (mpirank==0){
    cout << "Navier-Stokes is finished, now create global solution" << endl;
}
fespace VhGlobal1(Th1Global, [P2,P2,P2,P1]);
int[int] subIdxV  = restrict(Vh1, VhGlobal1, n2o1);

VhGlobal1 [uxGlobal,uyGlobal,uzGlobal,pGlobal], [uxsum,uysum,uzsum,psum];
ux[] .*= Ans.D;
if (mpirank==0){
    cout << "Ans.D" << endl;
}
uxGlobal[](subIdxV) = ux[]; //local function from #0 interpolated by 0, on the global domain
if (mpirank==0){
    cout << "uxGlobal[](subIdxV)" << endl;
}
mpiAllReduce(uxGlobal[], uxsum[], mpiCommWorld, mpiSUM); //problem because of global elements => need to use partition of unity

if (mpirank==0){
    cout << "Before save u" << endl;
    savevtk("$OUTPUT/u.vtu",Th1Global,uxsum,uysum,uzsum);
    saveArray("$OUTPUT/ux00.gp",uxsum[]);
    //real unmax = int2d(Th1Global,2)(1.*(uxsum<0.));
    real unmin=uxsum[].max; //1e5;
    for (int k=0; k<Th1Global.nbe; k++){
        if (Th1Global.be(k).label==2){
            int i0 = Th1Global.be(k)[0]; int i1 = Th1Global.be(k)[1]; int i2 = Th1Global.be(k)[2];
            real x0 = Th1Global(i0).x;  real y0 = Th1Global(i0).y; real z0 = Th1Global(i0).z; 
            real x1 = Th1Global(i1).x;  real y1 = Th1Global(i1).y; real z1 = Th1Global(i1).z; 
            real x2 = Th1Global(i2).x;  real y2 = Th1Global(i2).y; real z2 = Th1Global(i2).z;
            real un0 = uxsum(x0,y0,z0); real un1 = uxsum(x1,y1,z1); real un2 = uxsum(x2,y2,z2);
            real unminloc = min(un0, min(un1,un2)); 
            unmin = min(unminloc,unmin);
        }
    }
}

"""


interpolate_ns="""
                mesh3 Thold=readmesh3("$OUTPUT/ThFluid00.mesh");
                mesh3 ThBox=readmesh3("$MESH");
                mesh3 Th=trunc(ThBox, ((region == region1) || (region == region2)));
                savemesh(Th,"$OUTPUT/Th.mesh");
                mesh3 Thnew=trunc(Th, (region==region1));

                fespace Phold(Thold, [P2,P2,P2,P1]);
                fespace Phnew(Thnew, [P2,P2,P2,P1]);

                Phold [uxold,uyold,uzold,pold];
                Phnew [uxnew,uynew,uznew,pnew];

                readData("$OUTPUT/ux00.gp",uxold[]);
                [uxnew,uynew,uznew,pnew] = [uxold,uyold,uzold,pold];
                //[uxnew,uynew,uznew,pnew] = [(0.1^2 - (y-0.5)^2 - (z-0.5)^2)/0.01,0.,0.,0.];
                saveArray("$OUTPUT/ux.gp",uxnew[]);
                savesol("$OUTPUT/ux.sol",Thnew,uxnew);
                savesol("$OUTPUT/uy.sol",Thnew,uynew);
                savesol("$OUTPUT/uz.sol",Thnew,uznew);
            """

d1 = P1Function3D(Thinit, phi=phi1)
d2 = P1Function3D(Thinit, phi=phi2)
dOmega1 = mshdist(Thinit, phi=d1, ncpu=N)
dOmega2 = mshdist(Thinit, phi=d2, ncpu=N)
dOmega1.save(output+'/d1.sol')
dOmega2.save(output+'/d2.sol')
dOmega2.save(output+'/d2.o.sol')

FreeFemRunner([preamble, solve_ns],config,run_dir=output,run_file='solve_ns2.edp',debug=1).execute({'MESH':output+'/Th00.mesh'},ncpu=N)

Qlist = []

def computeU(mesh):
    FreeFemRunner([preamble, mesh_code, solve_kappa],config,run_dir=output,run_file='solve.edp',debug=1).execute({'MESH':mesh},ncpu=N)
    FreeFemRunner([preamble, interpolate_ns],config,run_dir=output,run_file='interpolate_ns.edp',debug=1).execute({'MESH':mesh})
    os.system(f"time mpiexec -n 1 ./MatrixBinaryPetscQ {output} -pc_type lu -pc_factor_mat_solver_type superlu_dist -ksp_monitor ")
    
    with open(output+'/J.gp','r') as f:
        J = float(f.readlines()[0])
    with open(output+'/V.gp','r') as f:
        V = float(f.readlines()[0])

    with open(output+'/Q.gp','r') as f:
        Q = float(f.readlines()[0])
    open(output+'/J.txt', 'w').close()
    Qlist.append(J)
    df = pd.DataFrame(data=np.array(Qlist),columns=['Q'])
    df.to_csv(output+'/J.txt', header=None, index=None, sep=' ', mode='a')
    return (Q,V)

code_sensitivity=r"""
load "PETSc"
macro dimension()3 //EOM
include "macro_ddm.idp"

mesh3 ThBox=readmesh3("$MESH");
mesh3 ThBoxGlobal = ThBox;
mesh3 Th2 = readmesh3("$OUTPUT/Th2.mesh");

int[int] n2oBox;
macro ThBoxN2O() n2oBox // this tells buildDmesh to keep the local to global correspondence
buildDmesh(ThBox)

// Variational space
fespace Ph2(Th2, P1);
fespace PhBox(ThBoxGlobal, P1);

PhBox dOmegaBox;
Ph2 T2,R2,T0;

T2[] = readsol("$OUTPUT/T2.sol");
R2[] = readsol("$OUTPUT/R2.sol");
T0 = 0.;
dOmegaBox[] = readsol("$OUTPUT/d2.sol");

if (mpirank==0){
    int[int] Order = [1];
    string DataName = "T";
    savevtk("$OUTPUT/T.vtu",Th2,T2,dataname=DataName, order=Order);
    //savevtk("$OUTPUT/R.vtu",Th2,R2);
}

fespace Gh1(ThBox,[P1,P1,P1]); // space for deformation field theta
Gh1 [thetaxJ,thetayJ,thetazJ], [thetaxH,thetayH,thetazH];
real gammaR = 10. * $hmin;


Gh1 [gx,gy,gz];
fespace PhParallel(ThBox,P1);
PhParallel kappa, dOmega;
int[int] subIdxK  = restrict(PhParallel, PhBox, n2oBox);

dOmega[] = dOmegaBox[](subIdxK);
fespace PhBox0(ThBox, P0);
PhBox0 nx2,ny2,nz2,norsqrt;
nx2 = dx(dOmega);
ny2 = dy(dOmega);
nz2 = dz(dOmega);
//norsqrt = max(sqrt(nx2^2 + ny2^2 + nz2^2),1e-6);
norsqrt = sqrt(nx2^2 + ny2^2 + nz2^2) + 1e-6;
nx2 = nx2 / norsqrt ;
ny2 = ny2 / norsqrt ;
nz2 = nz2 / norsqrt ;

varf probgsd([gx,gy,gz],[thetax, thetay,thetaz]) = int3d(ThBox,qforder=3)( tr([gx,gy,gz])*[thetax, thetay, thetaz]      )
                                                   +int3d(ThBox,qforder=2)( tr(grad(dOmega))*[thetax, thetay, thetaz]  );

Mat AG; //Global distributed matrix
{
macro def(i)[i, i#B, i#C] //EOM
macro init(i)[i, i, i] // EOM
createMat(ThBox, AG, [P1,P1,P1])
}

AG=probgsd(Gh1,Gh1,tgv=-1);
set(AG, sparams = "-ksp_max_it 99 -ksp_rtol 1e-6 -pc_type hypre"); //"-pc_type lu"
real[int] bG=probgsd(0,Gh1,tgv=-1);
gx[] = AG^-1 * bG;

Mat AK;
createMat(ThBox,AK,P1)
varf probkappa(kappa, psi) = int3d(ThBox,qforder=3)(9*ThBox.hmax*ThBox.hmax*tr(grad(kappa))*grad(psi) + kappa*psi)
                             +int3d(ThBox,qforder=2)(div(gx,gy,gz)*psi);


AK = probkappa(PhParallel, PhParallel, tgv=-1);
set(AK, sparams = "-pc_type lu " );
real[int] bK = probkappa(0, PhParallel, tgv=-1);
kappa[] = AK^-1 * bK;

func cut=min(max((y-0.5)^2 + (z-0.5)^2 - 0.2^2, x-0.02), max((y-0.5)^2 + (z-0.5)^2 - 0.2^2,0.98-x) );
PhParallel cutp1=cut<0;

/* //Using the functional J = int2d(Th,10)(beta^2*(T-T0)^2)

varf rieszJ([dpxJ,dpyJ,dpzJ],[thetaxp,thetayp,thetazp])= int3d(ThBox,qforder=3,optimize=0)(gammaR^2*GRAD(dpxJ,dpyJ,dpzJ)'*GRAD(thetaxp,thetayp,thetazp)+[dpxJ,dpyJ,dpzJ]'*[thetaxp,thetayp,thetazp])
                                            +int3d(ThBox)(tgv*cutp1*[dpxJ,dpyJ,dpzJ]'*[thetaxp,thetayp,thetazp])     
                                            +int2d(ThBox,10,optimize=0)(dot(thetaxp,thetayp,thetazp,nx2,ny2,nz2)*(beta^2*T2^2*(kappa - 4.*beta/k2) + beta*T2*R2*(2.*beta/k2 - kappa) - k2*tr(grad(T2))*grad(R2) ) )
                                           +on(1, 2, 5, 6, 7, 8, 9, 12, 13, 100, dpxJ=0.,dpyJ=0.,dpzJ=0.);

varf rieszH([dpxVol,dpyVol,dpzVol],[thetaxp,thetayp,thetazp])=int3d(ThBox,qforder=3,optimize=0)(gammaR^2*GRAD(dpxVol,dpyVol,dpzVol)'*GRAD(thetaxp,thetayp,thetazp)+[dpxVol,dpyVol,dpzVol]'*[thetaxp,thetayp,thetazp])
                                                            +int3d(ThBox)(tgv*cutp1*[dpxVol,dpyVol,dpzVol]'*[thetaxp,thetayp,thetazp])     
                                                            +int2d(ThBox,10,qforder=2,optimize=0)(dot(thetaxp,thetayp,thetazp,nx2,ny2,nz2))
                                                             +on(1, 2, 5, 6, 7, 8, 9, 12, 13, 100, dpxVol=0.,dpyVol=0.,dpzVol=0.);                                                             


varf rieszJ([dpxJ,dpyJ,dpzJ],[thetaxp,thetayp,thetazp])= int3d(ThBox,qforder=3,optimize=0)(gammaR^2*GRAD(dpxJ,dpyJ,dpzJ)'*GRAD(thetaxp,thetayp,thetazp)+[dpxJ,dpyJ,dpzJ]'*[thetaxp,thetayp,thetazp])
                                            +int2d(ThBox,10,qforder=3,optimize=0)(corr*divtauN(thetaxp,thetayp,thetazp,nx2,ny2,nz2)*beta^2*(T2-T0)^2 - corr*2.*beta^2*(T2-T0)*tr(grad(T0))*[thetaxp,thetayp,thetazp])
                                            +int3d(ThBox,region2,qforder=1,optimize=0)(-corr*div(thetaxp,thetayp,thetazp)*k2*tr(grad(T2))*grad(R2)+ corr*2*k2*tr(grad(R2))*MatrixByVector(EPS(thetaxp,thetayp,thetazp), grad(T2)))
                                            +int2d(ThBox,10,qforder=3,optimize=0)(-corr*divtauN(thetaxp,thetayp,thetazp,nx2,ny2,nz2)*beta*(T2-T0)*R2 + corr*beta*R2*tr(grad(T0))*[thetaxp,thetayp,thetazp])
                                           +on(1, 2, 5, 6, 7, 8, 9, 12, 13, 100, dpxJ=0.,dpyJ=0.,dpzJ=0.);

varf rieszH([dpxVol,dpyVol,dpzVol],[thetaxp,thetayp,thetazp])=int3d(ThBox,qforder=3,optimize=0)(gammaR^2*GRAD(dpxVol,dpyVol,dpzVol)'*GRAD(thetaxp,thetayp,thetazp)+[dpxVol,dpyVol,dpzVol]'*[thetaxp,thetayp,thetazp])
                                                            
                                                            +int3d(ThBox,region2,qforder=1,optimize=0)(div(thetaxp,thetayp,thetazp))
                                                             +on(1, 2, 5, 6, 7, 8, 9, 12, 13, 100, dpxVol=0.,dpyVol=0.,dpzVol=0.);                                                             
*/

//Surface shape derivative
varf rieszJ([dpxJ,dpyJ,dpzJ],[thetaxp,thetayp,thetazp])= int3d(ThBox,qforder=3,optimize=0)(gammaR^2*GRAD(dpxJ,dpyJ,dpzJ)'*GRAD(thetaxp,thetayp,thetazp)+[dpxJ,dpyJ,dpzJ]'*[thetaxp,thetayp,thetazp])
                                            
                                            +int2d(ThBox,10,optimize=0)(dot(thetaxp,thetayp,thetazp,nx2,ny2,nz2)*(beta*T2*(kappa - 2.*beta/k2) + beta*T2*R2*(2.*beta/k2 - kappa) - k2*tr(grad(T2))*grad(R2) ) )
                                           +on(1, 2, 5, 6, 7, 8, 9, 12, 13, 100, dpxJ=0.,dpyJ=0.,dpzJ=0.);

varf rieszH([dpxVol,dpyVol,dpzVol],[thetaxp,thetayp,thetazp])=int3d(ThBox,qforder=3,optimize=0)(gammaR^2*GRAD(dpxVol,dpyVol,dpzVol)'*GRAD(thetaxp,thetayp,thetazp)+[dpxVol,dpyVol,dpzVol]'*[thetaxp,thetayp,thetazp])
                                                            
                                                            +int2d(ThBox,10,qforder=2,optimize=0)(dot(thetaxp,thetayp,thetazp,nx2,ny2,nz2))
                                                             +on(1, 2, 5, 6, 7, 8, 9, 12, 13, 100, dpxVol=0.,dpyVol=0.,dpzVol=0.);                                                             


//Volume shape derivative
/*varf rieszJ([dpxJ,dpyJ,dpzJ],[thetaxp,thetayp,thetazp])= int3d(ThBox,qforder=3,optimize=0)(gammaR^2*GRAD(dpxJ,dpyJ,dpzJ)'*GRAD(thetaxp,thetayp,thetazp)+[dpxJ,dpyJ,dpzJ]'*[thetaxp,thetayp,thetazp])
                                            +int2d(ThBox,10,qforder=3,optimize=0)(corr*divtauN(thetaxp,thetayp,thetazp,nx2,ny2,nz2)*beta*(T2-T0) - corr*beta*tr(grad(T0))*[thetaxp,thetayp,thetazp])
                                            +int3d(ThBox,region2,qforder=1,optimize=0)(-corr*div(thetaxp,thetayp,thetazp)*k2*tr(grad(T2))*grad(R2)+ corr*2*k2*tr(grad(R2))*MatrixByVector(EPS(thetaxp,thetayp,thetazp), grad(T2)))
                                            +int2d(ThBox,10,qforder=3,optimize=0)(-corr*divtauN(thetaxp,thetayp,thetazp,nx2,ny2,nz2)*beta*(T2-T0)*R2 + corr*beta*R2*tr(grad(T0))*[thetaxp,thetayp,thetazp])
                                           +on(1, 2, 5, 6, 7, 8, 9, 12, 13, 100, dpxJ=0.,dpyJ=0.,dpzJ=0.);

varf rieszH([dpxVol,dpyVol,dpzVol],[thetaxp,thetayp,thetazp])=int3d(ThBox,qforder=3,optimize=0)(gammaR^2*GRAD(dpxVol,dpyVol,dpzVol)'*GRAD(thetaxp,thetayp,thetazp)+[dpxVol,dpyVol,dpzVol]'*[thetaxp,thetayp,thetazp])
                                                            +int3d(ThBox,region2,qforder=1,optimize=0)(div(thetaxp,thetayp,thetazp))
                                                             +on(1, 2, 5, 6, 7, 8, 9, 12, 13, 100, dpxVol=0.,dpyVol=0.,dpzVol=0.);                                                             
*/

 //Using the real heat loss functional Q = int2d(Th,10)(beta*(T-T0)) 


Mat A; //Global distributed matrix
{
macro def(i)[i, i#B, i#C] //EOM
macro init(i)[i, i, i] // EOM
createMat(ThBox, A, [P1,P1,P1])
}

A=rieszJ(Gh1,Gh1,tgv=-1);
set(A, sparams = "-ksp_max_it 99 -ksp_rtol 1e-10 -pc_type hypre"); //"-pc_type lu"
real[int] diffH=rieszH(0,Gh1);
thetaxH[] = A^-1 * diffH;
real[int] diffJ=rieszJ(0,Gh1);
thetaxJ[] = A^-1*diffJ;

fespace Gh1Global(ThBoxGlobal, [P1,P1,P1]);
Gh1Global [thetaxJGlob,thetayJGlob,thetazJGlob],  [thetaxHGlob,thetayHGlob,thetazHGlob], [thetaxJGlobSum,thetayJGlobSum,thetazJGlobSum], [thetaxHGlobSum,thetayHGlobSum,thetazHGlobSum];

int[int] subIdx = restrict(Gh1, Gh1Global, n2oBox);

//I want to convert gradJ from local to global!

thetaxJ[] .*= A.D;
thetaxJGlob[](subIdx) = thetaxJ[];
mpiAllReduce(thetaxJGlob[], thetaxJGlobSum[], mpiCommWorld, mpiSUM);

thetaxH[] .*= A.D;
thetaxHGlob[](subIdx) = thetaxH[];
mpiAllReduce(thetaxHGlob[], thetaxHGlobSum[], mpiCommWorld, mpiSUM);

//1st option
real[int] diffJGlob(Gh1Global.ndof), diffGGlob(Gh1Global.ndof), diffHGlob(Gh1Global.ndof), diffJGlobSum(Gh1Global.ndof), diffGGlobSum(Gh1Global.ndof), diffHGlobSum(Gh1Global.ndof);
diffJGlob(subIdx) = diffJ;
diffHGlob(subIdx) = diffH;

mpiAllReduce(diffJGlob, diffJGlobSum, mpiCommWorld, mpiSUM);
mpiAllReduce(diffHGlob, diffHGlobSum, mpiCommWorld, mpiSUM);

if (mpirank==0){
    saveArray("$OUTPUT/diffJ.gp",diffJGlobSum);
    saveArray("$OUTPUT/diffV.gp",diffHGlobSum);
    saveArray("$OUTPUT/gradJ.gp",thetaxJGlobSum[]);
    saveArray("$OUTPUT/gradV.gp",thetaxHGlobSum[]);

    //savevtk("$OUTPUT/gradJ.vtu",ThBoxGlobal, [-thetaxJGlobSum,-thetayJGlobSum,-thetazJGlobSum]);
}

"""
def sensitivity(mesh):
    FreeFemRunner([preamble,code_sensitivity],config,run_dir=output,run_file='sensitivities.edp',debug=1).execute(config={'MESH':mesh,'hmin':hmin},ncpu=N)
    dJ=readFFArray(output+'/diffJ.gp')
    dV=readFFArray(output+'/diffV.gp')
    gradJ=readFFArray(output+'/gradJ.gp')
    gradV=readFFArray(output+'/gradV.gp')
    return (dJ,dV,gradJ,gradV)


def shapeUpdate(mesh, xi):
    M = Mesh3D(mesh)
    lvlset = P1Function3D(M, output+'/d2.o.sol')
    phi = mshdist(M,phi=lvlset,ncpu=N) # Initial level set

    # Advection of the level set
    vel = P1Vector3D(M,[xi[::3],xi[1::3],xi[2::3]])
    phiNew = advect(M,phi,vel,T=1.) #T=1.

    # Enforce non optimizable regions
    #cutP1 = P1Function3D(M,cut)
    #phiNew.sol = np.maximum(phiNew.sol,cutP1.sol)
    #barrierP1 = P1Function3D(M,barrier)
    #phiNew.sol = np.minimum(phiNew.sol, barrierP1.sol)
    phiNew.save(output+"/d2.sol")

    codeRegularize="""
    load "PETSc"
    macro dimension()3 //EOM
    include "macro_ddm.idp"
    
    mesh3 Th=readmesh3("$MESH");
    mesh3 ThGlobal = Th;

    int[int] n2o;
    macro ThN2O() n2o // this tells buildDmesh to keep the local to global correspondence
    buildDmesh(Th)

    fespace Fh(Th,P1); fespace FhGlobal(ThGlobal,P1);
    FhGlobal phi0Global, phiGlobal, phisum; Fh phi,phi0;
    phi0Global[]=readsol("$OUTPUT/d2.sol");

    real gamma=0.01*$hmin;

    varf regularize(phi,psi)=int3d(Th)(gamma^2*(dx(phi)*dx(psi)+dy(phi)*dy(psi)+dz(phi)*dz(psi) )+phi*psi)
                            +int3d(Th)(phi0*psi);

    int[int] subIdx  = restrict(Fh, FhGlobal, n2o);
    phi0[] = phi0Global[](subIdx);

    Mat A;
    createMat(Th,A,P1)

    A = regularize(Fh, Fh, tgv=-1);
    set(A, sparams = "-pc_type lu " );
    real[int] b = regularize(0, Fh, tgv=-1);
    phi[] = A^-1 * b;

    phi[] .*= A.D;
    phiGlobal[](subIdx) = phi[]; //local function from #0 interpolated by 0, on the global domain
    mpiAllReduce(phiGlobal[], phisum[], mpiCommWorld, mpiSUM); //problem because of global elements => need to use partition of unity

    if (mpirank==0){
        //save the solution
        savesol("$OUTPUT/d2.sol",ThGlobal,phisum);
    }
    """
    FreeFemRunner([preamble,codeRegularize],config,run_dir=output,run_file='regularize.edp',debug=1).execute({'MESH':mesh,'hmin':hmin},ncpu=N)

    # Remesh
    #10 Triangles 1e-2 2e-2 1e-3
    refermmg = f"""
                LSReferences
                2
                2 3 200
                20 nosplit

                Parameters
                6

                100 Triangles 5e-3 1e-2 1e-3
                10 Triangles 1e-2 2e-2 1e-3
                1 Triangles 5e-3 1e-2 1e-3
                2 Triangles 5e-3 1e-2 1e-3
                12 Triangles 5e-3 1e-2 1e-3
                13 Triangles 5e-3 1e-2 1e-3
    """
    newM=mmg3d(M,hmin,hmax,hgrad,hausd,ls=True,sol=phiNew,params=refermmg,extra_args="-rmc 1e-8 ",debug=1)
    newM.tetrahedra[newM.tetrahedra[:,4] == 200, 4] = 2
    newM.save(output+'/ThBox.o.mesh')
    
    interpolate="""
                    load "medit"
                    load "msh3"
                    mesh3 Thold=readmesh3("$MESH");
                    mesh3 Thnew=readmesh3("$OUTPUT/ThBox.o.mesh");

                    fespace Phold(Thold, P1);
                    fespace Phnew(Thnew, P1);

                    Phold d2old;
                    Phnew d2new;

                    d2old[] = readsol("$OUTPUT/d2.sol");
                    d2new = d2old;
                    savesol("$OUTPUT/d2.sol",Thnew,d2new);
                """
    #FreeFemRunner(interpolate,config,run_dir=output,run_file='interpolate.edp',debug=1).execute({'MESH':mesh})

    #phiNew = P1Function3D(newM, output+'/d2.sol')
    #d2 = mshdist(newM, phi=phiNew, ncpu=N)

    #modif
    newM.tetrahedra[newM.tetrahedra[:,4] == 20, 4] = 3
    d2 = mshdist(newM, ncpu=N)

    d2.save(output+'/d2.sol')
    phiNew = P1Function3D(newM, phi1)
    d1 = mshdist(newM, phi=phiNew, ncpu=N)
    d1.save(output+'/d1.sol')

    return output+'/ThBox.o.mesh'


class PureThermalHeatConduction(Optimizable):
    def __init__(self):
        super().__init__()
        self.ucomputed = False
        self.sensitivity_computed = False
        self.obj = None
        self.nconstraints = 1
        self.nineqconstraints = 0

    # Initialization
    def x0(self):
        return output+"/Th00.mesh"

    # Objective function and constraints
    def evalObjective(self, x):
        if not self.ucomputed:
            (J,H)=computeU(x)
            self.ucomputed = True
            self.obj = (J,H)
        return self.obj

    def J(self,x):
        (J,H)=self.evalObjective(x)
        return J

    def G(self,x):
        (J,H)=self.evalObjective(x)
        return [H]

    # Shape derivatives, sensitivity of objective and constraint
    def evalSensitivities(self,x):
        if not self.sensitivity_computed:
            (dJ,dH,gradJ,gradH)=sensitivity(x)
            self.sensitivity_computed = True
            self.sensitivities = (dJ,dH,gradJ,gradH)
        return self.sensitivities

    def dJ(self,x):
        (dJ,dH,gradJ,gradH) = self.evalSensitivities(x)
        return dJ

    def dG(self,x):
        (dJ,dH,gradJ,gradH) = self.evalSensitivities(x)
        return [dH]

    #Gradient and transpose
    def dJT(self,x):
        (dJ,dH,gradJ,gradH) = self.evalSensitivities(x)
        return gradJ

    def dGT(self,x):
        (dJ,dH,gradJ,gradH) = self.evalSensitivities(x)
        return np.asarray([gradH]).T

    #Retraction : shape update
    def retract(self, x, dx):
        # Assume that the next computations will be performed on a new mesh
        self.sensitivity_computed = False
        self.ucomputed = False
        return shapeUpdate(x,dx)

    #Accept step : save the current result
    def accept(self,results):
        mesh = results['x'][-1]
        niter = len(results['J'])-1
        itf = format(niter,'03d')
        shutil.copyfile(mesh,output+'/Th_'+itf+'.mesh')
        shutil.copyfile(output+'/d2.sol',output+'/d2.o.sol')
        M = Mesh3D(mesh)
        if niter >= 1:
            shutil.copyfile(output+'/T.vtu',output+'/T_'+format(niter-1,'03d')+'.vtu')
            #shutil.copyfile(output+'/gradJ.vtu',output+'/gradJ_'+itf+'.vtu')
        shutil.copyfile(output+'/ThBox.vtu',output+'/ThBox_'+itf+'.vtu')
        results['x'][-1] = output+'/Th_'+itf+'.mesh'

        open(output+'/Q.txt', 'w').close()
        df = pd.DataFrame(data=np.array(results['J']),columns=['J'])
        df.to_csv(output+'/Q.txt', header=None, index=None, sep=' ', mode='a')

        open(output+'/G.txt', 'w').close()
        df = pd.DataFrame(data=np.array(results['G']),columns=['G'])
        df.to_csv(output+'/G.txt', header=None, index=None, sep=',', mode='a')


optSettings = {'dt':hmin,'alphaJ':2.,'alphaC':1., 'maxit':150, 'provide_gradient': True, 'maxtrials':1,
              'itnormalisation': 3}
results=nlspace_solve(PureThermalHeatConduction(), optSettings)

