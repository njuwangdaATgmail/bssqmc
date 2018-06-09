I am designing a BSS auxiliary field DQMC program.

### How to design a user-friendly interface?

+ The fermionic lattice can be input by giving unit cell structure and its spatial arrangement. For periodic boundary condition, (La,Lb,Lc) are enough to pin down the lattice. For a general open boundary condition, one can realize it as a large super unit cell. 
+ The quatic fermion-fermion interactions are always be decoupled using descrete HS fields. The interaction values can be used to obtain gamma and lambda, which of couse should also be given by hand, e.g. a file or an external subroutine. In addition, we also need provide the basis $[i_1,i_2,\cdots]$ of fmat, which can also be given by a file or external subroutine. *For better interface, I prefer to use files as input, which will not break the completeness of the source code.*
+ The continuous boson fields can be used to realize phonon or other types of boson fields. *Can we avoid to write an external subroutine and then recompile the source code?*

### One solution: *dqmc.x* and *dqmc.a*
+ If the boson action has already been realized within the package (e.g. Ising or phonon or phi4), we can run the executable program *dqmc.x* directly. 
+ Otherwise, we should use *dqmc.a* as a library.

### How to input the boson field?
+ nb_field can be filled by finding site+shiftb. 
+ mask_field_site should be correctly set to realize the 'checkerboard' structure. How to realize it? For example, one can input (orb,supa,moda,supb,modb), then only when mod(a,supa)==moda and mod(b,supb)==modb, the mask_fielda_site is set true.
+ fmat of cousre can be input directly.

### List of benchmarks and targets

+ SU(N) Hubbard model: U is decoupled in charge channel
+ SU(2) Hubbard-Holstein model: U is decoupled in magnetic channel and Vph is in charge channel.
+ bilayer UVJ-ph-model: all interactions are on the rungs.
+ ...
+ a general multi-orbital model with/without spin orbit couplings
+ other boson fields should be able to be easily added, *e.g.* spin/charge/nematic-fluctuations, gauge fields. 

