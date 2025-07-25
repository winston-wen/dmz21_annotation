注释加解密流程, 考证group update.

DONE:
* 验证 $$a, b, \Delta_k, \Delta_p$$ 的关系.
* 考证出group和group update分别对应数学上的 $$\left<g\right>$$ 和 $$\left<f\right>$$.
* 考证出 `discrete_log_f` 对应[CL15, Proposition 1]中的 $$L(m)$$.
* 考证出选取参数 $$p, q$$ 的原理, 以及 $$g$$ 的阶数.
* 考证出 $$\varphi_p$$ 位于[CL15, Appendix B.1].
* 考证出 $$\varphi_p$$ 及其逆运算是怎么算的.
* 考证出 $$\mathfrak{r}$$ 是怎么构造的[HJPT98, Sec 3.1].
* 复述为什么 $$\left<f\right>$$ 中的元素具有简单的离散代数, 在什么条件下简单?

简述: 因为 $$f^m=[(p^2, m^{-1}p)]$$. 在任何条件下都简单, 因为明文是b参数的模p逆元.
但是, 乘了 $$h^r$$ 之后, 就只对私钥持有者简单.

TODO:

* 找到类群零元内元素都是主理想, 零元外元素都不是主理想的理论依据.
* 理解理想相乘的高效算法, 如NUCOMP. 在Rust版本中可以先不做, 用现成的. 在Go版本中要做.