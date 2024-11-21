# entangled subspace

[![DOI](https://zenodo.org/badge/639765941.svg)](https://zenodo.org/doi/10.5281/zenodo.7995115)

This repository is to reproduce results in the paper [arxiv-link](https://arxiv.org/abs/2311.10353) "Quantifying Subspace Entanglement with Geometric Measures"

🚀 Exciting News! We've launched the `numqi` package [github/numqi](https://github.com/numqi/numqi), combining all the functionalities of this repository and even more! 🌟 To dive into these features, just install `numqi` using `pip install numqi`, and explore the relevant functions within the package. 🛠️

## quick start

This project is based on `numqi` package [github-link](https://github.com/numqi/numqi), which is a package for numerical optimization for quantum information. To install `numqi`, run the following command (see [numqi-doc](https://numqi.github.io/numqi/installation/) for detailed installation instruction):

```bash
pip install numqi
```

Then, clone this repository and cd into the directory

```bash
git clone git@github.com:numqi/entangled-subspace.git
cd entangled-subspace
```

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

| (n,k) | analytical | err(GD) | err(seesaw) | err(PPT) | time(GD) | time(seesaw) | time(PPT) |
| :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| (5,1)| 0.5904 | $2.1\times 10^{-13}$ | $2.3\times 10^{-13}$ | $1.4\times 10^{-10}$ | 0.103 | 0.0101 | 1.73 |
| (5,2)| 0.6544 | $2.7\times 10^{-14}$ | $7.8\times 10^{-13}$ | $5.1\times 10^{-11}$ | 0.087 | 0.0091 | 1.64 |
| (5,3)| 0.6544 | $1.3\times 10^{-13}$ | $7.0\times 10^{-13}$ | $5.0\times 10^{-11}$ | 0.084 | 0.0091 | 1.59 |
| (6,1)| 0.5981 | $2.1\times 10^{-13}$ | $2.7\times 10^{-13}$ | $8.9\times 10^{-10}$ | 0.129 | 0.017 | 25.1 |
| (6,2)| 0.6708 | $1.7\times 10^{-14}$ | $4.4\times 10^{-13}$ | $3.4\times 10^{-11}$ | 0.122 | 0.018 | 26.5 |
| (6,3)| 0.6875 | $1.2\times 10^{-13}$ | $6.5\times 10^{-13}$ | $2.6\times 10^{-13}$ | 0.166 | 0.018 | 25.7 |
| (7,1)| 0.6034 | $1.2\times 10^{-14}$ | $1.1\times 10^{-12}$ | None | 0.173 | 0.019 | None |
| (7,2)| 0.6813 | $2.2\times 10^{-13}$ | $2.4\times 10^{-13}$ | None | 0.163 | 0.019 | None |
| (7,3)| 0.7062 | $4.2\times 10^{-14}$ | $1.6\times 10^{-13}$ | None | 0.150 | 0.020 | None |
