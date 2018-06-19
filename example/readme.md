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
  + small V will greatly reduce the sign problem of Hubbard-Holstein model?

Their suitable combinations, e.g.:
+ SU(2N) Kane-Mele-Hubbard-Jperp-Peierls model

### questions
+ Is Majorana representation necessary to eliminate sign prblem in some specific models? Can we realize all the Majorana-QMC models?

### Multi-orbital Coulomb interactions
+ Kanamori parameters: (see e.g. Hotta, 2006; XY Deng, thesis) 
  + J=J' only without spin-orbit coupling
  + U=V+J+J' only for rotational symmetrical system (e.g. among Eg or t2g orbitals, or isolated atom)

+ Take J=J', these two terms can be combin as a square form:
$$
\frac{1}{2}\left(\sum_\sigma a^\dagger_\sigma b_\sigma+h.c.\right)^2=\left(a_\uparrow^\dagger a_\downarrow^\dagger b_\downarrow b_\uparrow +h.c.\right)-\left( \sum_{\sigma \sigma'} a_\sigma^\dagger a_{\sigma '}b_{\sigma'}^\dagger b_\sigma \right)
$$
based on which complex HS decoupling can be performed. Spin is still a good quantum number. 
+ In special case of $J'=0$, the Hund's coupling term can be written as $-J\vec{S}\cdot\vec{S}$, which can be decoupled via real HS field. This may weaken the sign problem. 
+ The density-density interactions can be decoupled in density channel directly. 

$$
-\frac{V}{2}\sum_{a<b} \left(n_a-n_b\right)^2-\frac{U}{2}\sum_a\left(n_{a\uparrow}-n_{a\downarrow}\right)^2
$$

+ In fact, the V-term can also be written as $V(\sum_a n_a-N)^2$, which however yields a complex HS field and may bring severe sign problem. Nevertheless, it deserves a try.