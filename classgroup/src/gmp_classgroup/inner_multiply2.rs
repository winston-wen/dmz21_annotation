impl GmpClassGroup {
pub fn into_raw(self) -> (Mpz, Mpz) {
    (self.a, self.b)
}

// 出处: [Cohen1993, Algorithm 5.4.9] NUCOMP算法, 计算两个二次型的复合.
// 原理: [Cohen1993, Definition 5.4.6, Section 5.2] 二次型的复合就是理想的乘.
// TODO: 看起来不像NUCOMP, 因为里面一个if-else都没有.
// TODO: 看起来更像[Cohen1993, Algorithm 5.4.7]
fn inner_multiply2(&mut self, rhs: &Self, ctx: &mut Ctx) {
    self.assert_valid();
    rhs.assert_valid();

    // $$g = (b_1 + b_2) / 2$$
    ffi::mpz_add(&mut ctx.congruence_context.g, &self.b, &rhs.b);
    ffi::mpz_fdiv_q_ui_self(&mut ctx.congruence_context.g, 2);

    // $$h = (b_2 - b_1) / 2$$
    ffi::mpz_sub(&mut ctx.h, &rhs.b, &self.b);
    ffi::mpz_fdiv_q_ui_self(&mut ctx.h, 2);

    debug_assert!(&ctx.h + &ctx.congruence_context.g == rhs.b);
    debug_assert!(&ctx.congruence_context.g - &ctx.h == self.b);

    // w := gcd(a1, a2, g)
    ffi::three_gcd(&mut ctx.w, &self.a, &rhs.a, &ctx.congruence_context.g);

    // s = a1/w
    ffi::mpz_fdiv_q(&mut ctx.s, &self.a, &ctx.w);

    // t = a2/w
    ffi::mpz_fdiv_q(&mut ctx.t, &rhs.a, &ctx.w);

    // 至此已经可以计算 $$A = st = a^1a^2/w^2$$.
    // 对照一下[Cohen1993, Lemma 5.4.5], 发现 $$A$$ 少了系数 $$d_0$$.
    // 据此推测, 该函数假设输入中的至少一个二次型是 primitive (系数互质) 的.

    // u = g/w
    ffi::mpz_fdiv_q(&mut ctx.u, &ctx.congruence_context.g, &ctx.w);

    // a = t*u = a2*g/w^2
    ffi::mpz_mul(&mut ctx.a, &ctx.t, &ctx.u);

    // b = h*u + s*c1 = (g*h+a1*c1)/w
    ffi::mpz_mul(&mut ctx.b, &ctx.h, &ctx.u);
    ffi::mpz_mul(&mut ctx.m, &ctx.s, &self.c);
    ctx.b += &ctx.m; // &mut ctx.b

    // m = s*t = a1*a2 / w^2
    ffi::mpz_mul(&mut ctx.m, &ctx.s, &ctx.t); 

    // 求解 mu 使得 t*u*mu = h*u + s*c1 (mod s)
    ctx.congruence_context.solve_linear_congruence(
        &mut ctx.mu,
        Some(&mut ctx.v), // v = m / gcd(a, m) = a1 / gcd(a1, g)
        &ctx.a, // a = t*u
        &ctx.b, // b = h*u + s*c1
        &ctx.m, // m = s*t
    );

    // a = t*v = a1*a2 / (gcd(a1,a2,g) * gcd(a1, g))
    ffi::mpz_mul(&mut ctx.a, &ctx.t, &ctx.v);

    // b = h - t * mu
    ffi::mpz_mul(&mut ctx.m, &ctx.t, &ctx.mu);
    ffi::mpz_sub(&mut ctx.b, &ctx.h, &ctx.m);

    // m = s = a1 / w
    ctx.m.set(&ctx.s); // &mut ctx.m

    // 求解 lambda 使得 t*v*lambda = h - t * mu (mod s)
    ctx.congruence_context.solve_linear_congruence(
        &mut ctx.lambda,
        Some(&mut ctx.sigma),
        &ctx.a, // a = t*v
        &ctx.b, // b = h - t*mu
        &ctx.m, // m = s
    );

    // k = mu + v*lambda 
    ffi::mpz_mul(&mut ctx.a, &ctx.v, &ctx.lambda);
    ffi::mpz_add(&mut ctx.k, &ctx.mu, &ctx.a);

    // l = (k*t - h)/s
    ffi::mpz_mul(&mut ctx.l, &ctx.k, &ctx.t);
    ffi::mpz_sub(&mut ctx.v, &ctx.l, &ctx.h); 
    ffi::mpz_fdiv_q(&mut ctx.l, &ctx.v, &ctx.s); 
    // 此后再没用过 `ctx.v`. 验算时可以放心地使用第一次赋值的结果.

    // m = (t*u*k - h*u - c*s) / s*t
    ffi::mpz_mul(&mut ctx.m, &ctx.t, &ctx.u);
    ctx.m *= &ctx.k; // &mut ctx.m1
    ffi::mpz_mul(&mut ctx.a, &ctx.h, &ctx.u);
    ctx.m -= &ctx.a; // &mut ctx.m
    ffi::mpz_mul(&mut ctx.a, &self.c, &ctx.s);
    ctx.m -= &ctx.a; // &mut ctx.m
    ffi::mpz_mul(&mut ctx.a, &ctx.s, &ctx.t);
    ffi::mpz_fdiv_q(&mut ctx.lambda, &ctx.m, &ctx.a);

    // A = s*t
    ffi::mpz_mul(&mut self.a, &ctx.s, &ctx.t);

    // B = w*u - k*t - l*s = g+h-2k*t = b2 - 2k*t
    ffi::mpz_mul(&mut self.b, &ctx.w, &ctx.u);
    ffi::mpz_mul(&mut ctx.a, &ctx.k, &ctx.t);
    self.b -= &ctx.a; // &mut self.b
    ffi::mpz_mul(&mut ctx.a, &ctx.l, &ctx.s);
    self.b -= &ctx.a; // &mut self.b

    // C = k*l - w*m
    ffi::mpz_mul(&mut self.c, &ctx.k, &ctx.l);
    ffi::mpz_mul(&mut ctx.a, &ctx.w, &ctx.lambda);
    self.c -= &ctx.a; // &mut self.c

    self.inner_reduce(ctx);
}