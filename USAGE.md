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

## Using MiLaDy potentials in Lammps

Running Lammps with MiLaDy potentials is easy: set `pair_style  milady` and indicate the potential file in your Lammps input file:

```
#Lammps input

pair_style milady
pair_coeff * * Fe_LML.pot Fe 
```
The ready-to-use input files can be found in `Examples` of this repository

## [Installation](INSTALLATION.md)
