In this folder, we provide some examples which can be directly used or slightly modified for practical applications. 

The name rule of each example is 'lattice_-model'
+ lattice: square, honeycomb, triangle, kagome, CuO2, bi-square, ...
+ model: hubbard, spin-orbit coupling (soc), general Coulomb (UVJ), phonon(Holstein, Peierls, Frohlich...)

+ The program prepare_input can help one prepare the input file. However, it's not very smart now, e.g. checking the validity of parameters.

### models without sign problem which can be simulated in this code

arbitrary-filled artitrary lattice: 
+ SU(2N) Holstein/Peierls/negative-U model

half-filled bipartite lattice: 
+ SU(2N) Hubbard/Jperp-model
+ Kane-Mele_-Hubbard model
+ SU(N) Jdimer-model (when N=1, it becomes spinless V-model)
+ SU(N) Holstein/Peierls/negative-U model
? UV-model

Their suitable combinations, e.g.:
+ SU(2N) Kane-Mele-Hubbard-Jperp-Peierls model
