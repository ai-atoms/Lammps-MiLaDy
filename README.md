# Lammps-MiLaDy
ML interatomic potentials developed with MiLaDy package and Lammps module to perform calculations with them 


Installation:
-----------------


1. Download the last version of Lammps

```
git clone --recursive  https://github.com/lammps/lammps.git  lammps.git 
```

2. Download this package:

```
git clone --recursive  git@github.com:ai-atoms/Lammps-MiLaDy.git Lammps-MiLaDy.git
```

3. copy USER-MILADY interface and milady library in main Lammps src directory and lib directory, respectively:

```
cp -rp Lammps-MiLaDy.git/USER-MILADY  lammps.git/src/
cp -rp Lammps-MiLaDy.git/milady       lammps.git/lib/
```

4. install MiLaDy in Lammps

```
cd lammps.git/src 
make yes-user-milady
```

5. Is time to choose the comppilator intel or gfortran for MiLaDy. Note that there is no any restriction on the choice that you have for the compilation 
of Lammps (we use default parameters and g++ compiler). Let's take the case Intel Fortran. There are two substeps.

5a. Make the appropiate links for Makefile and milady library for intel fortran. 
```
cp lammps.git/lib/milady
ln -s Makefile.lammps.intel  Makefile.lammps
ln -s libmilady.a_intel_mkl_oneAPI2021  libmilady.a
```

5b. Edit the Makefile.lammps for your architecture. MiLaDy library need MKL and some Intel libraries from Intel compiler. You need to localize the 
root directory for MKL and intel64 libraries. We use oneAPI Intel distribution. And here is our choice: 

```
 MKLROOT=/opt/intel/oneapi/mkl/latest/
 LIBCOMP=/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64/
```
There are similar path for any older distribution of MKL and Intel Fortran such as Intel Composer, Intel Parallel Studio  etc

If you have doubts. Write us. We are happy to help you !!!!

