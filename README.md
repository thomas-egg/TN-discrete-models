# Tensor Network Methods for Statistical Physics
Tensor Networks provide a convenient way to factorize high dimensional tensors into a contraction over a chain or **network** of lower dimensional tensors. For example,
given some high dimensional tensor $`T`$ with six indices, we can express this tensor as a sum over lower dimensional tensors by:

$$
T^{s_1 s_2 s_3 s_4 s_5 s_6} = \sum_{{\alpha}} A_{\alpha_1}^{s_1} A_{\alpha_1 \alpha_2}^{s_2} A_{\alpha_2 \alpha_3}^{s_3} A_{\alpha_3 \alpha_4}^{s_4} A_{\alpha_4 \alpha_5}^{s_5} A_{\alpha_5 \alpha_6}^{s_6}
$$

Indices $`s`$ are physical indices, and indices $`\alpha`$ are nonphysical **bond indices** that we sum over. This particular factorization is a matrix product state/tensor train (MPSS/TT) factorization. Tensor network methods are generally used in Quantum Many-Body physics, but I explore realizations of tensor networks for statistical physics problems in this repository.

## Ising Model and TRG
### The 2D Ising Model
The 2D Ising model is one of the most well-studied models in statistical physics. This model describes a 2D lattice of interacting spins fully described by the following Hamiltonian:

$$
  H = - J \sum_{i=1}^N \sum_{j=1}^N \left[ \sigma_{i,j} \sigma_{i+1,j} + \sigma_{i,j} \sigma_{i,j+1} \right] - h \sum_{i=1}^N \sum_{j=1}^N \sigma_{i,j}
$$

where $`\sigma_{i,j} \in {-1,1}`$ is the spin at site i, j in the 2D lattice. The partition function of this model was solved by Lars Onsager in 1944 but we will study how to solve for the partition function numerically using the **Tensor Renormalization Group (TRG)** algorithm. 
### TRG
TRG is a coarse-graining process by which a numerical estimate for the partition function of some lattice model in the thermodynamic limit can be solved for numerically. At a high-level, the steps of the algorithm can be outlined as follows:
- Express the lattice as a lattice of tensors, each with 4 physical indices corresponding to the 4 nearest neighbor spins
- In each tensor, store the Boltzmann weight of every possible spin configuration (for Ising: $`2^4`$)
- Factorize each of these tensors (SVD) and contract over new indices.
- Repeat this decomposition -> contraction process until only **one** tensor remains
- The double trace of this tensor is the partition function!

## References
Original solution for 2D Ising: https://journals.aps.org/pr/abstract/10.1103/PhysRev.65.117

Original TRG: https://arxiv.org/pdf/cond-mat/0611687

TRG for Ising Model: https://tensornetwork.org/trg/
