import os
import pickle
import time
import numpy as np
import opt_einsum
from tqdm import tqdm

import numqi

np_rng = np.random.default_rng()
hf_randc = lambda *sz: np_rng.normal(size=sz) + 1j * np_rng.normal(size=sz)
hf_norm = lambda x: x / np.linalg.norm(x.reshape(-1), ord=2)


def get_GME_pure_seesaw(np0:np.ndarray, converge_eps:float=1e-7, num_repeat:int=1, maxiter=1000, psi_list=None, seed=None):
    # seesaw algorithm doi.org/10.1103/PhysRevA.84.022323
    N0 = np0.ndim
    assert (N0>=2) and all(x>=2 for x in np0.shape)
    np0 = np0 / np.linalg.norm(np0.reshape(-1), ord=2)
    if num_repeat>1:
        assert psi_list is None, "num_repeat>1 is not supported with psi_list"
    dim_list = np0.shape
    contract_expr_list = []
    for ind0 in range(N0):
        tmp0 = [((dim_list[x],),(x,)) for x in range(N0) if x!=ind0]
        tmp1 = [y for x in tmp0 for y in x]
        contract_expr_list.append(opt_einsum.contract_expression(np0.conj(), list(range(N0)), *tmp1, [ind0], constants=[0]))
    ret_repeat_list = []
    for _ in range(num_repeat):
        if (num_repeat>1) or (psi_list is None):
            np_rng = np.random.default_rng(seed)
            psi_list = []
            for ind0 in range(N0):
                tmp0 = np_rng.normal(size=dim_list[ind0]) + 1j * np_rng.normal(size=dim_list[ind0])
                psi_list.append(tmp0/np.linalg.norm(tmp0, ord=2))
        else:
            assert all(x.shape==(y,) for x,y in zip(psi_list, dim_list)), "Inconsistent shape"
        for _ in range(maxiter):
            fidelity_list = np.zeros(N0, dtype=np.float64)
            for ind0 in range(N0):
                tmp0 = [psi_list[x] for x in range(N0) if x!=ind0]
                tmp1 = contract_expr_list[ind0](*tmp0).conj()
                ret_fidelity = np.vdot(tmp1, tmp1).real
                tmp2 = tmp1 / np.sqrt(ret_fidelity)
                fidelity_list[ind0] = abs(np.vdot(tmp2, psi_list[ind0]))**2
                psi_list[ind0] = tmp2
            if (1-fidelity_list.min())<converge_eps:
                break
        ret_repeat_list.append((max(0,1-ret_fidelity), psi_list))
    ret = min(ret_repeat_list, key=lambda x: x[0])
    return ret


def get_GME_subspace_seesaw(np0:np.ndarray, converge_eps:float=1e-7, num_repeat:int=1, maxiter=1000, psi_list=None, seed=None):
    # seesaw algorithm doi.org/10.1103/PhysRevA.84.022323
    N0 = np0.ndim-1
    assert (N0>=2) and all(x>=2 for x in np0.shape[1:])
    dim_list = np0.shape[1:]
    np0 = numqi.matrix_space.reduce_vector_space(np0.reshape(np0.shape[0],-1), zero_eps=1e-10).reshape(-1, *dim_list)
    N1 = np0.shape[0]
    if N1==1:
        assert False, "N1==1 is not supported" #TODO
    if num_repeat>1:
        assert psi_list is None, "num_repeat>1 is not supported with psi_list"
    dim_list = np0.shape[1:]
    contract_expr_list = []
    for ind0 in range(N0):
        tmp0 = [((dim_list[x],),(x+1,)) for x in range(N0) if x!=ind0]
        tmp1 = [y for x in tmp0 for y in x]
        contract_expr_list.append(opt_einsum.contract_expression(np0.conj(), list(range(N0+1)), *tmp1, [0,ind0+1], constants=[0]))
    ret_repeat_list = []
    for _ in range(num_repeat):
        if (num_repeat>1) or (psi_list is None):
            np_rng = np.random.default_rng(seed)
            psi_list = []
            for ind0 in range(N0):
                tmp0 = np_rng.normal(size=dim_list[ind0]) + 1j * np_rng.normal(size=dim_list[ind0])
                psi_list.append(tmp0/np.linalg.norm(tmp0, ord=2))
        else:
            assert all(x.shape==(y,) for x,y in zip(psi_list, dim_list)), "Inconsistent shape"
        for _ in range(maxiter):
            fidelity_list = np.zeros(N0, dtype=np.float64)
            for ind0 in range(N0):
                tmp0 = [psi_list[x] for x in range(N0) if x!=ind0]
                tmp1 = contract_expr_list[ind0](*tmp0)
                EVL,EVC = np.linalg.eigh(tmp1.T.conj() @ tmp1)
                EVCi = EVC[:,-1]
                fidelity_list[ind0] = abs(np.vdot(EVCi, psi_list[ind0]))**2
                gme = max(0, 1-EVL[-1])
                psi_list[ind0] = EVCi
            if (1-fidelity_list.min())<converge_eps:
                break
        ret_repeat_list.append((gme, psi_list))
    ret = min(ret_repeat_list, key=lambda x: x[0])
    return ret


def test_get_GME_pure_seesaw_dicke():
    num_qubit = 4
    klist = list(range(num_qubit+1))

    ret_ = np.array([numqi.state.get_qubit_dicke_state_GME(num_qubit,x) for x in klist])
    np0_list = [numqi.dicke.Dicke(num_qubit-x,x) for x in klist]
    ret0 = np.array([get_GME_pure_seesaw(x.reshape([2]*num_qubit), converge_eps=1e-10, maxiter=1000, num_repeat=1)[0] for x in np0_list])
    assert np.abs(ret_-ret0).max() < 1e-7

def test_get_GME_pure_seesaw_Wtype():
    tmp0 = np_rng.normal(size=(23,3))
    coeff_list = tmp0 / np.linalg.norm(tmp0, ord=2, axis=1)[:,None]

    ret_ = np.array([numqi.state.get_Wtype_state_GME(*x) for x in coeff_list])
    kwargs = dict(converge_eps=1e-10, maxiter=1000, num_repeat=1)
    ret0 = np.array([get_GME_pure_seesaw(numqi.state.Wtype(x).reshape(2,2,2), **kwargs)[0] for x in coeff_list])
    assert np.abs(ret_-ret0).max() < 1e-7


def demo_compare_gradient_descent_with_seesaw():
    case_list = [(2,2,2), (2,2,4),(2,2,6),(2,3,4),(2,3,6),(2,3,8),(3,3,6),(3,3,8),(3,4,7),(4,4,7),(4,5,10)]
    kwargs_gd = dict(theta0='uniform', num_repeat=3, tol=1e-12, print_every_round=0)
    kwargs_seesaw = dict(converge_eps=1e-10, num_repeat=3, maxiter=2000)
    for dimA,dimB,dimC in case_list:
        np_list = numqi.matrix_space.get_completed_entangled_subspace((dimA, dimB, dimC), kind='quant-ph/0409032')[0]
        model = numqi.matrix_space.DetectCanonicalPolyadicRankModel((dimA, dimB, dimC), rank=1)
        model.set_target(np_list)
        tmp0 = time.time()
        theta_optim1 = numqi.optimize.minimize(model, **kwargs_gd)
        t0 = time.time() - tmp0
        tmp0 = time.time()
        gme,psi_list = get_GME_subspace_seesaw(np_list, **kwargs_seesaw)
        t1 = time.time() - tmp0
        print(f'[{dimA}x{dimB}x{dimC}] loss(GD)={theta_optim1.fun}, delta={theta_optim1.fun-gme}, time(GD)= {t0:.3f}, time(Seesaw)={t1:.3f}')


def demo_compare_with_seesaw_dicke():
    datapath = 'data/compare_with_seesaw_dicke.pkl'
    if os.path.exists(datapath):
        with open(datapath, 'rb') as fid:
            import pickle
            all_data = pickle.load(fid)
    else:
        n_k_list = [(5,1),(5,2),(5,3),(6,1),(6,2),(6,3),(7,1),(7,2),(7,3)]
        all_data = dict()
        kwargs_gd = dict(theta0='uniform', tol=1e-10, num_repeat=3, print_every_round=0)
        kwargs_seesaw = dict(converge_eps=1e-10, num_repeat=3, maxiter=2000)
        for num_qubit,dicke_k in n_k_list:
            np0 = numqi.dicke.Dicke(num_qubit-dicke_k, dicke_k).reshape([2]*num_qubit)
            ret_ = numqi.state.get_qubit_dicke_state_GME(num_qubit, dicke_k)
            model = numqi.matrix_space.DetectCanonicalPolyadicRankModel([2]*num_qubit, rank=1)
            model.set_target(np0)
            model() # precompute
            tmp0 = time.time()
            ret_gd = numqi.optimize.minimize(model, **kwargs_gd).fun
            t0 = time.time() - tmp0
            tmp0 = time.time()
            ret_seesaw = get_GME_pure_seesaw(np0, **kwargs_seesaw)[0]
            t1 = time.time() - tmp0
            all_data[(num_qubit,dicke_k)] = (ret_gd, ret_seesaw, ret_, t0, t1)
        with open(datapath, 'wb') as fid:
            import pickle
            pickle.dump(all_data, fid)

    for (num_qubit,dicke_k), (ret_gd, ret_seesaw, ret_, t0, t1) in all_data.items():
        print(f'[Dicke({num_qubit},{dicke_k})] analytical={ret_} error(GD)={ret_gd-ret_}, error(seesaw)={ret_seesaw-ret_}, time(GD)={t0:.5f}, time(Seesaw)={t1:.5f}')


if __name__=='__main__':
    # demo_compare_gradient_descent_with_seesaw()
    demo_compare_with_seesaw_dicke()
