# Tensor Network Methods for Statistical Physics
Tensor Networks provide a convenient way to factorize high dimensional tensors into a contraction over a chain or **network** of lower dimensional tensors. For example,
given some high dimensional tensor $`T`$ with six indices, we can express this tensor as a sum over lower dimensional tensors by:
$$
  T^{s_1 s_2 s_3 s_4 s_5 s_6 = \Sum_{{\alpha}} A_{\alpha1}^{s_1} A_{\alpha1 \alpha2}^{s_2} A_{\alpha2 \alpha3}^{s_3} A_{\alpha3 \alpha4}^{s_4} A_{\alpha4 \alpha5}^{s_5} A_{\alpha5 \alpha6}^{s_6}
$$
