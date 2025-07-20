注释加解密流程, 考证group update.

DONE:
* 验证 $$a, b, \Delta_k, \Delta_p$$ 的关系.
* 考证出group和group update分别对应数学上的 $$\left<g\right>$$ 和 $$\left<f\right>$$.
* 考证出 `discrete_log_f` 对应[CL15, Proposition 1]中的 $$L(m)$$.
* 考证出选取参数 $$p, q$$ 的原理, 以及 $$g$$ 的阶数.
* 考证出 $$\varphi_p$$ 位于[CL15, Appendix B.1].

TODO:
* 考证 $$\varphi_p$$ 及其逆运算是怎么算的. 特别留意在[CL15, Appendix B.1]的最后一段中,
    * 理想的数乘 & 相加 & 相乘, 怎么表示为二次型?
    * 理想与整数互质满足什么性质, 有什么可计算的判定方法?
* 考证 $$\mathfrak{r}$$ 是怎么构造的[HJPT98, Sec 3.1].
* 复述为什么 $$\left<f\right>$$ 中的元素具有简单的离散代数, 在什么条件下简单?
