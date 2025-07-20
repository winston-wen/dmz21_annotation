注释加解密流程, 考证group update.

DONE:
* 验证 $$a, b, \Delta_k, \Delta_p$$ 的关系.
* 考证出group和group update分别对应数学上的 $$\left<g\right>$$ 和 $$\left<f\right>$$.
* 考证出 `discrete_log_f` 对应[CL15, Proposition 1]中的 $$L(m)$$.
* 考证出选取参数 $$p, q$$ 的原理, 以及 $$g$$ 的阶数.

TODO:
* 考证 $$\varphi_p^{-1}$$ 是怎么算的.
* 考证 $$\mathfrak{r}$$ 是怎么构造的[HJPT98, Sec 3.1].
* 清晰地复述出为什么 `discrete_log_f` 是那么做的.