impl GmpClassGroup {
pub fn into_raw(self) -> (Mpz, Mpz) {
    (self.a, self.b)
}

fn inner_multiply2(&mut self, rhs: &Self, ctx: &mut Ctx) {
    self.assert_valid();
    rhs.assert_valid();

    ffi::mpz_add(&mut ctx.congruence_context.g, &self.b, &rhs.b);
    ffi::mpz_fdiv_q_ui_self(&mut ctx.congruence_context.g, 2);

    ffi::mpz_sub(&mut ctx.h, &rhs.b, &self.b);
    ffi::mpz_fdiv_q_ui_self(&mut ctx.h, 2);

    debug_assert!(&ctx.h + &ctx.congruence_context.g == rhs.b);
    debug_assert!(&ctx.congruence_context.g - &ctx.h == self.b);

    ffi::three_gcd(&mut ctx.w, &self.a, &rhs.a, &ctx.congruence_context.g);
    ffi::mpz_fdiv_q(&mut ctx.s, &self.a, &ctx.w);
    ffi::mpz_fdiv_q(&mut ctx.t, &rhs.a, &ctx.w);
    ffi::mpz_fdiv_q(&mut ctx.u, &ctx.congruence_context.g, &ctx.w);
    ffi::mpz_mul(&mut ctx.a, &ctx.t, &ctx.u);
    ffi::mpz_mul(&mut ctx.b, &ctx.h, &ctx.u);
    ffi::mpz_mul(&mut ctx.m, &ctx.s, &self.c);
    ctx.b += &ctx.m;

    ffi::mpz_mul(&mut ctx.m, &ctx.s, &ctx.t); 

    ctx.congruence_context.solve_linear_congruence(
        &mut ctx.mu,
        Some(&mut ctx.v),
        &ctx.a,
        &ctx.b,
        &ctx.m,
    );

    ffi::mpz_mul(&mut ctx.a, &ctx.t, &ctx.v);

    ffi::mpz_mul(&mut ctx.m, &ctx.t, &ctx.mu);
    ffi::mpz_sub(&mut ctx.b, &ctx.h, &ctx.m);

    ctx.m.set(&ctx.s);

    ctx.congruence_context.solve_linear_congruence(
        &mut ctx.lambda,
        Some(&mut ctx.sigma),
        &ctx.a,
        &ctx.b,
        &ctx.m,
    );

    ffi::mpz_mul(&mut ctx.a, &ctx.v, &ctx.lambda);
    ffi::mpz_add(&mut ctx.k, &ctx.mu, &ctx.a);

    ffi::mpz_mul(&mut ctx.l, &ctx.k, &ctx.t);
    ffi::mpz_sub(&mut ctx.v, &ctx.l, &ctx.h);
    ffi::mpz_fdiv_q(&mut ctx.l, &ctx.v, &ctx.s);

    ffi::mpz_mul(&mut ctx.m, &ctx.t, &ctx.u);
    ctx.m *= &ctx.k;
    ffi::mpz_mul(&mut ctx.a, &ctx.h, &ctx.u);
    ctx.m -= &ctx.a;
    ffi::mpz_mul(&mut ctx.a, &self.c, &ctx.s);
    ctx.m -= &ctx.a;
    ffi::mpz_mul(&mut ctx.a, &ctx.s, &ctx.t);
    ffi::mpz_fdiv_q(&mut ctx.lambda, &ctx.m, &ctx.a);

    ffi::mpz_mul(&mut self.a, &ctx.s, &ctx.t);

    ffi::mpz_mul(&mut self.b, &ctx.w, &ctx.u);
    ffi::mpz_mul(&mut ctx.a, &ctx.k, &ctx.t);
    self.b -= &ctx.a;
    ffi::mpz_mul(&mut ctx.a, &ctx.l, &ctx.s);
    self.b -= &ctx.a;

    ffi::mpz_mul(&mut self.c, &ctx.k, &ctx.l);
    ffi::mpz_mul(&mut ctx.a, &ctx.w, &ctx.lambda);
    self.c -= &ctx.a;

    self.inner_reduce(ctx);
}