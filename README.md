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

3. copy USER-MILADY interface in main src Lammps directory:

```
cp -rp Lammps-MiLaDy.git/USER-MILADY  lammps.git/src
```

4. install MiLaDy in Lammps

```
cd lammps.git/src 
make yes-user-milady 
```
