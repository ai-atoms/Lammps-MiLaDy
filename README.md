# Lammps-MiLaDy
This directory provides ML interatomic potentials developed with MiLaDy package and Lammps module to perform calculations with them 


Using MiLaDy potentials in Lammps:
--------------------------------------
Running Lammps with MiLaDy potentials is easy: set `pair_style  milady` and indicate the potential file in your Lammps input file:

```
#Lammps input

pair_style milady
pair_coeff * * LML.pot Fe 
```
More deailed input can be found in `Examples` of this repositiry


Installation:
-----------------


1. Download the last version of Lammps

```
git clone --recursive  https://github.com/lammps/lammps.git  lammps.git 
```

2. Download Lammps-MiLaDy package:

```
git clone --recursive  git@github.com:ai-atoms/Lammps-MiLaDy.git Lammps-MiLaDy.git
```

3. Copy USER-MILADY interface and milady library in main src directory and lib directory, respectively:

```
cp -rp Lammps-MiLaDy.git/USER-MILADY  lammps.git/src/
cp -rp Lammps-MiLaDy.git/milady       lammps.git/lib/
```

4. Install `MiLaDy` in `Lammps`

```
cd lammps.git/src 
make yes-user-milady
```

5. Is time to choose the comppilator `intel` or `gfortran` for `MiLaDy`. Note that there is no restriction in the choice that you have for the compilation 
of `Lammps` (we use default parameters and `g++` compiler). Now, let's take the case of `MiLaDy`integrated in `Lammps`  using Intel Fortran `ifort`. There are two steps to do that:

5a. Make the appropiate links for Makefile and milady library for `ifort`. 
```
cd lammps.git/lib/milady
ln -s Makefile.lammps.intel  Makefile.lammps
ln -s libmilady.a_intel_mkl_oneAPI2021  libmilady.a
```

5b. Edit the `Makefile.lammps` for your architecture. `MiLaDy` library need MKL and some Intel libraries from `Intel Fortran` compiler. You need to localize the 
root directory for MKL and intel64 libraries. We use `oneAPI Intel` distribution. Here are our choices: 

```
 MKLROOT=/opt/intel/oneapi/mkl/latest/
 LIBCOMP=/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64/
```
There are similar path for any older distribution of MKL and Intel Fortran such as Intel Composer, Intel Parallel Studio  etc

If you have doubts: write us. We are happy to help you !!!!

6. Compile again Lammps (e.g. using traditional `make mpi`). 
7. To test your Lammps-MiLady, run the examples provided in provided in `Examples` of this repositiry


