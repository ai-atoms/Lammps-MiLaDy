# Lammps-MiLaDy
This directory provides ML interatomic potentials developed with MiLaDy package and Lammps module to perform calculations with them. 

```
,---.    ,---.-./`)   .---.       ____    ______        ____     __  
|    \  /    \ .-.')  | ,_|     .'  __ `.|    _ `''.    \   \   /  / 
|  ,  \/  ,  / `-' \,-./  )    /   '  \  \ _ | ) _  \    \  _. /  '  
|  |\_   /|  |`-'`"`\  '_ '`)  |___|  /  |( ''_'  ) |     _( )_ .'   
|  _( )_/ |  |.---.  > (_)  )     _.-`   | . (_) `. | ___(_ o _)'    
| (_ o _) |  ||   | (  .  .-'  .'   _    |(_    ._) '|   |(_,_)'     
|  (_,_)  |  ||   |  `-'`-'|___|  _( )_  |  (_.\.' / |   `-'  /      
|  |      |  ||   |   |        \ (_ o _) /       .'   \      /       
'--'      '--''---'   `--------`'.(_,_).''-----'`      `-..-'        
```
MiLaDy library has many contributors, such as:

Alexandra M. Goryaeva, Université Paris-Saclay, CEA, Service de Recherches de Métallurgie Physique, 91191, Gif-sur-Yvette, France (`alexandra . goryaeva at cea dot fr`) 

Thomas  D. Swinburne, Aix-Marseille Université, CNRS, CINaM UMR 7325, Campus de Luminy, 13288 Marseille, France (`swinburne  at cinam dot univ - mrs dot fr`)

Clovis Lapointe, Université Paris-Saclay, CEA, Service de Recherches de Métallurgie Physique, 91191, Gif-sur-Yvette, France (`clovis . lapointe at cea dot fr`) 

Mihai-Cosmin Marinica, Université Paris-Saclay, CEA, Service de Recherches de Métallurgie Physique, 91191, Gif-sur-Yvette, France(`mihai - cosmin . marinica at cea dot fr`)

The interface with Lammps was written by TDS and MCM. 


Using MiLaDy potentials in Lammps:
--------------------------------------
Running Lammps with MiLaDy potentials is easy: set `pair_style  milady` and indicate the potential file in your Lammps input file:

```
#Lammps input

pair_style milady
pair_coeff * * Fe_LML.pot Fe 
```
The ready-to-use input files can be found in `Examples` of this repository


Installation:
-----------------


1. Download the last version of `Lammps`

```
git clone --recursive  https://github.com/lammps/lammps.git  lammps.git 
```

2. Download our `Lammps-MiLaDy` package:

```
git clone --recursive  git@github.com:ai-atoms/Lammps-MiLaDy.git Lammps-MiLaDy.git
```

3. Copy `USER-MILADY` interface and `milady` library in main `src` and `lib` directories of `Lammps`, respectively:

```
cp -rp Lammps-MiLaDy.git/USER-MILADY  lammps.git/src/
cp -rp Lammps-MiLaDy.git/milady       lammps.git/lib/
```

4. Install `MiLaDy` in `Lammps`

```
cd lammps.git/src 
make yes-user-milady
```

5. Choose the comppilator for `MiLaDy`: `intel` or `gfortran`. Note that there is no restriction in the choice that you have for the compilation 
of `Lammps` (we use default parameters and `g++` compiler). 
    Now, let's take the case of `MiLaDy`integrated in `Lammps`  using Intel Fortran `ifort`. There are two steps to do that:
    
    5a. Make the appropiate links for `Makefile` and `milady` library for `ifort`. 
    
    ```
    cd lammps.git/lib/milady
    ln -s Makefile.lammps.intel  Makefile.lammps
    ln -s libmilady.a_intel_mkl_oneAPI2021  libmilady.a
    ```
    
    5b. Edit the `Makefile.lammps` for your architecture. `MiLaDy` library uses `MKL` and some other Intel libraries from `Intel Fortran` compiler. You need to   localize the root directory for `MKL` and `intel64` libraries. We use `oneAPI Intel` free distribution. Here are our choices:

   ```
    MKLROOT=/opt/intel/oneapi/mkl/latest/
    LIBCOMP=/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64/
   ```
   The paths will be similar for any older distribution of `MKL` and `Intel Fortran`,  such as `Intel Composer`, `Intel Parallel Studio`,  etc.
   If you have doubts: write us. We are happy to help you !!!!

6. Compile your `Lammps` (e.g. using traditional `make mpi`). 
7. To test your `Lammps-MiLady`, run the examples provided in `Examples` of this repositiry


