All UVJ-terms are decoupled using complex HS fields.
$$
\frac{V}{2} \left(\sum_{a\sigma} n_{a\sigma} -N\right)^2+\frac{U-V}{2}\sum_a\left(\sum_\sigma n_{a\sigma}-1\right)^2+\frac{J}{2}\sum_{a<b}\left(\sum_{\sigma} a^\dagger_\sigma b_\sigma+h.c.\right)^2.
$$


## benchmarked by exact diagonalization

We consider a two-site 3-orbital model with U=4, V=2, J=1 and compare the QMC results with exact diagonalization. 
QMC parameters are: $T=0$, $\beta=16$, $\Delta\tau=0.05$, nwarmup=1000, nmeasure=10000, nbin=10, ncore=2 (see input file *dqmc.in_benchmark*)

|                                                  |       QMC        |    ED     |
| :----------------------------------------------: | :--------------: | :-------: |
| $\langle a_{1\sigma}^\dagger a_{2\sigma}\rangle$ | $0.233\pm0.001$  | $0.2336$  |
| $\langle n_{1a\uparrow}n_{1a\downarrow}\rangle$  | $0.039\pm0.004$  | $0.0405$  |
| $\langle n_{1a\uparrow}n_{2a\downarrow}\rangle$  | $0.385\pm0.003$  | $0.3850$  |
| $\langle n_{1a\uparrow}n_{1b\downarrow}\rangle$  | $0.180\pm0.002$  | $0.1799$  |
| $\langle n_{1a\uparrow}n_{2b\downarrow}\rangle$  | $0.356\pm0.002$  | $0.3573$  |
|        $\langle S_{1a}^+S_{1a}^-\rangle$         | $0.460\pm0.004$  | $0.4595$  |
|        $\langle S_{1a}^+S_{2a}^-\rangle$         | $-0.308\pm0.003$ | $-0.3099$ |
|        $\langle S_{1a}^+S_{1b}^-\rangle$         | $0.137\pm0.002$  | $0.1377$  |
|        $\langle S_{1a}^+S_{2b}^-\rangle$         | $-0.211\pm0.002$ | $-0.2125$ |
|        $\langle d_{12a} d_{12a} \rangle $        |  $3.08\pm0.02$   | $3.0940$  |
|        $\langle d_{12a} d_{12b} \rangle $        | $1.032\pm0.007$  | $1.0349$  |

