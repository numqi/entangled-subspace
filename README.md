# entangled subspace

This repository is to reproduce results in the paper [arxiv-link](https://arxiv.org/abs/2311.10353) "Quantifying Subspace Entanglement with Geometric Measures"

PS: this folder previously located in `project/entangled-subspace` of the `github/numqi/numqi` repository.

1. link
   * [github/tensorspace](https://github.com/thetensor-space/TensorSpace) magma based
2. completely entangled subspace (CES)
   * maximal dimension in space $\bigotimes_{i=1}^{k}\mathcal{H}_k$: $d_1d_2\cdots d_k-(d_1+\cdots+d_k)+k-1$ [doi-link](https://doi.org/10.1007/BF02829441)
3. perfectly entangled subspace

`demo_schmidt_rank.py`, `GenTiles1`

| dim | `loss(1)` | `loss(2)` | time (second) |
| :-: | :-: | :-: | :-: |
| `4` | $0.030$ | $9.2\times 10^{-15}$ | `0.82` |
| `8` | $0.016$ | $2.2\times 10^{-12}$ | `2.5` |
| `16` | $0.0079$ | $2.3\times 10^{-12}$ | `3.4` |
| `32` | $0.0038$ | $4.3\times 10^{-12}$ | `20.0` |
| `64` | $0.0018$ | $5.8\times 10^{-12}$ | `352` |

`draft_seesaw.py/demo_compare_with_seesaw_dicke()`

| (n,k) | analytical | err(GD) | err(seesaw) | time(GD) | time(seesaw) |
| :-: | :-: | :-: | :-: | :-: | :-: |
| (5,1) | $0.5904$ | $5.6\times 10^{-13}$ | $2.3\times 10^{-13}$ | `0.098` | `0.0092` |
| (5,2) | $0.6544$ | $1.1\times 10^{-13}$ | $5.1\times 10^{-13}$ | `0.088` | `0.0095` |
| (5,3) | $0.6544$ | $9.8\times 10^{-13}$ | $3.4\times 10^{-13}$ | `0.123` | `0.0099` |
| (6,1) | $0.5981$ | $1.8\times 10^{-13}$ | $1.8\times 10^{-13}$ | `0.120` | `0.017` |
| (6,2) | $0.6708$ | $2.2\times 10^{-13}$ | $3.4\times 10^{-13}$ | `0.120` | `0.016` |
| (6,3) | $0.6875$ | $6.8\times 10^{-14}$ | $3.8\times 10^{-13}$ | `0.116` | `0.017` |
| (7,1) | $0.6034$ | $8.8\times 10^{-14}$ | $2.3\times 10^{-13}$ | `0.128` | `0.018` |
| (7,2) | $0.6813$ | $3.6\times 10^{-12}$ | $1.7\times 10^{-13}$ | `0.138` | `0.018` |
| (7,3) | $0.7062$ | $7.9\times 10^{-12}$ | $2.7\times 10^{-13}$ | `0.167` | `0.018` |
