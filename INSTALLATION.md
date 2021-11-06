# Lammps-MiLaDy

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

## [About Lammps-MiLaDy](README.md)
## [Using MiLaDy potentials in Lammps](USAGE.md)

## Installation

**Requirements:**

- Fortran compiler: `gfortran` OR `ifort`
- C++ compiler: `g++` 
- MKL library from Intel 

We have tested with following versions:

```
  g++        (version >= 9.0.0)
  gfortran   (version >= 9.0.0)
  ifort      (version >= 2018.0.0)
  MKL        (version >= 2018.0.0) with gfortran support
```


**Instructions**

1. Clone the public repository of `Lammps`:

```
git clone --recursive  https://github.com/lammps/lammps.git  lammps.git 
```

2. Clone our `Lammps-MiLaDy` repository:

```
git clone --recursive  https://github.com/ai-atoms/Lammps-MiLaDy.git Lammps-MiLaDy.git

```

3. Copy `USER-MILADY` interface and `milady` library to main `src` and `lib` directories of `Lammps`, respectively:

```
cp -rp Lammps-MiLaDy.git/USER-MILADY  lammps.git/src/
cp -rp Lammps-MiLaDy.git/milady       lammps.git/lib/
```

4. Install `MiLaDy` in `Lammps`

```
cd lammps.git/src 
make yes-user-milady
```

5. Choose the compilator for `MiLaDy`: `intel` or `gfortran`. Note that there is no restriction in the choice that you have for the compilation 
of `Lammps` (we use default parameters and `g++` compiler). 
Below, we provide the 2-step example using Intel Fortran `ifort`. The case of `gfortran` is similar.
    
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


6. Compile `Lammps`. For more details on`Lammps` compilation, see [Lammps documentation](https://docs.lammps.org/Install.html)
 ```
 make mpi
```
 
7. Test your `Lammps-MiLady`, run the examples provided in `Examples` of this repository

IMPORTANT:  Many thanks to users that have reported some typos in this installation page: Marie Landeiro Dos Reis and Antoine Kraych ! 

