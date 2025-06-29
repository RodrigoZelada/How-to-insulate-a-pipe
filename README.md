# How-to-insulate-a-pipe

Requirements:
- petsc version >= 3.11.3 configured with SupoerLU_DIST : https://petsc.org/release/manualpages/Mat/MATSOLVERSUPERLU_DIST/
- FreeFem++ (configured with PETSc): https://doc.freefem.org/introduction/installation.html#compilation-on-ubuntu
- python >= 3.8, with pymedit, pyfreefem, nullspace_optimizer: https://gitlab.com/florian.feppon/pyfreefem, https://gitlab.com/florian.feppon/pymedit, https://gitlab.com/florian.feppon/null-space-optimizer
- ISCDtoolbox (medit, advection, mshdist): https://github.com/iscdtoolbox](https://github.com/ISCDtoolbox/Mshdist, https://github.com/ISCDtoolbox/Advection, https://github.com/ISCDtoolbox/Medit

  Once you have installed everything, you have to compile the C++ code by doing:

  make MatrixBinaryPetsc
  make MatrixBinaryPetscQ

  and then, you can run the test cases by doing:
  
  python IsolantCylinder.py
  (or the others .py files)

The base codes used are cited in our article:
F. Caubet, C. Conca, M. Dambrine, and R. Zelada. How to insulate a pipe?.
