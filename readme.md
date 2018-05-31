I am designing a BSS auxiliary field DQMC program.

### How to design a user-friendly interface?

+ The fermionic lattice can be input by giving unit cell structure and its spatial arrangement. For periodic boundary condition, (La,Lb,Lc) are enough to pin down the lattice. For a general open boundary condition, one can realize it as a large super unit cell. 
+ The quatic fermion-fermion interactions are always be decoupled using descrete HS fields. The interaction values can be used to obtain gamma and lambda, which of couse should also be given by hand, e.g. a file or an external subroutine. In addition, we also need provide the basis $[i_1,i_2,\cdots]$ of fmat, which can also be given by a file or external subroutine. *For better interface, I prefer to use files as input, which will not break the completeness of the source code.*
+ The continuous boson fields can be used to realize phonon or other types of boson fields. *Can we avoid to write an external subroutine and then recompile the source code?*

### List of benchmarks

+ SU(N) Hubbard model: 
+ SU(2) Hubbard-Holstein model: 
+ bilayer UVJ-model
+ 