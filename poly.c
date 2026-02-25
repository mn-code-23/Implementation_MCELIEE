#include <stdlib.h>
#include <string.h>
#include "poly.h"

/* ==================== HELPERS INTERNES ==================== */
void poly_copy(poly_t *dst, const poly_t *src) {
  if (!dst || !src || !dst->coeff || !src->coeff) return;
  dst->deg = src->deg;
  size_t cp = (size_t)(src->deg + 1);
  if (cp > (size_t)dst->size) cp = (size_t)dst->size;
  memcpy(dst->coeff, src->coeff, cp * sizeof(gf_t));
  for (size_t i = cp; i < (size_t)dst->size; i++) dst->coeff[i] = 0;
}

static void poly_zero(poly_t *p) {
  if (!p || !p->coeff) return;
  memset(p->coeff, 0, (size_t)p->size * sizeof(gf_t));
  p->deg = -1;
}

/* ===== Opérations polynômiales de base ===== */

poly_t poly_alloc(int max_deg) {
  poly_t p;
  if (max_deg < 0) max_deg = 0;
  p.size = (size_t)max_deg + 1;
  p.deg = -1;
  p.coeff = (gf_t *)calloc(p.size, sizeof(gf_t));
  return p;
}

void poly_free(poly_t *p) {
  if (!p) return;
  free(p->coeff);
  p->coeff = NULL;
  p->size = 0;
  p->deg = -1;
}

void poly_normalize(poly_t *p) {
  if (!p || !p->coeff) return;
  while (p->deg >= 0 && p->coeff[p->deg] == 0) {
    p->deg--;
  }
}

void poly_add(poly_t *r, const poly_t *a, const poly_t *b) {
  if (!r || !a || !b) return;
  int maxd = (a->deg > b->deg) ? a->deg : b->deg;
  if ((size_t)maxd + 1 > r->size) return;
  for (int i = 0; i <= maxd; i++) {
    gf_t av = (i <= a->deg) ? a->coeff[i] : 0;
    gf_t bv = (i <= b->deg) ? b->coeff[i] : 0;
    r->coeff[i] = (gf_t)(av ^ bv);
  }
  r->deg = maxd;
  poly_normalize(r);
}

void poly_mul(poly_t *r, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx) {
  if (!r || !a || !b || !ctx) return;
  int deg = a->deg + b->deg;
  if (deg < 0) {
    poly_zero(r);
    return;
  }
  if ((size_t)deg + 1 > r->size) deg = (int)r->size - 1;
  for (int i = 0; i <= deg; i++) r->coeff[i] = 0;
  for (int i = 0; i <= a->deg; i++) {
    for (int j = 0; j <= b->deg; j++) {
      if (i + j > deg) break;
      gf_t prod = gf_mul(ctx, a->coeff[i], b->coeff[j]);
      r->coeff[i + j] = gf_add(r->coeff[i + j], prod);
    }
  }
  r->deg = deg;
  poly_normalize(r);
}

void poly_divrem(poly_t *q, poly_t *r, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx) {
  if (!q || !r || !a || !b || !ctx) return;
  poly_t rem = poly_alloc(a->deg);
  poly_copy(&rem, a);
  poly_zero(q);
  while (rem.deg >= b->deg && rem.deg >= 0 && b->deg >= 0) {
    int shift = rem.deg - b->deg;
    gf_t coef = gf_mul(ctx, rem.coeff[rem.deg], gf_inv(ctx, b->coeff[b->deg]));
    if (shift < (int)q->size) q->coeff[shift] = coef;
    poly_t temp = poly_alloc(rem.deg);
    poly_zero(&temp);
    for (int i = 0; i <= b->deg; i++) {
      if (i + shift < (int)temp.size)
        temp.coeff[i + shift] = gf_mul(ctx, coef, b->coeff[i]);
    }
    poly_add(&rem, &rem, &temp); /* soustraction = addition */
    poly_normalize(&rem);
    poly_free(&temp);
  }
  poly_copy(r, &rem);
  poly_normalize(q);
  poly_free(&rem);
}

void poly_mod(poly_t *r, const poly_t *a, const poly_t *g, const gf_ctx_t *ctx) {
  if (!r || !a || !g || !ctx) return;
  poly_t qtmp = poly_alloc(a->deg);
  poly_divrem(&qtmp, r, a, g, ctx);
  poly_free(&qtmp);
}

void poly_gcd(poly_t *d, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx) {
  if (!d || !a || !b || !ctx) return;
  poly_t aa = poly_alloc(a->deg);
  poly_t bb = poly_alloc(b->deg);
  poly_copy(&aa, a);
  poly_copy(&bb, b);
  while (!(bb.deg < 0 || (bb.deg == 0 && bb.coeff[0] == 0))) {
    poly_t rtmp = poly_alloc(bb.deg);
    poly_divrem(NULL, &rtmp, &aa, &bb, ctx);
    poly_copy(&aa, &bb);
    poly_copy(&bb, &rtmp);
    poly_free(&rtmp);
  }
  poly_copy(d, &aa);
  poly_free(&aa);
  poly_free(&bb);
}

gf_t poly_eval(const poly_t *p, gf_t alpha, const gf_ctx_t *ctx) {
  if (!p || !ctx) return 0;
  gf_t res = 0;
  gf_t pow = 1;
  for (int i = 0; i <= p->deg; i++) {
    gf_t term = gf_mul(ctx, p->coeff[i], pow);
    res = gf_add(res, term);
    pow = gf_mul(ctx, pow, alpha);
  }
  return res;
}

static void poly_set_one(poly_t *p) {
  poly_zero(p);
  if (p->size > 0) {
    p->coeff[0] = 1;
    p->deg = 0;
  }
}

/* ==================== XGCD CORRIGÉ ==================== */

/* Forward declaration required by decrypt.c */

void poly_xgcd(poly_t *d, poly_t *u, poly_t *v, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx) {
  if (!d || !u || !v || !a || !b || !ctx) return;

  /* allocations temporaires de degré maximal a->deg + b->deg */
  int max_deg = a->deg + b->deg;
  poly_t old_r = poly_alloc(max_deg);
  poly_t r     = poly_alloc(max_deg);
  poly_t old_s = poly_alloc(max_deg);
  poly_t s     = poly_alloc(max_deg);
  poly_t old_t = poly_alloc(max_deg);
  poly_t t     = poly_alloc(max_deg);
  poly_t quotient = poly_alloc(max_deg);
  poly_t remainder = poly_alloc(max_deg);
  poly_t temp = poly_alloc(max_deg);
  poly_t product = poly_alloc(max_deg);

  /* initialisation */
  poly_copy(&old_r, a);
  poly_copy(&r, b);
  poly_set_one(&old_s);
  poly_zero(&s);
  poly_zero(&old_t);
  poly_set_one(&t);

  while (r.deg >= 0 && (r.deg > 0 || r.coeff[0] != 0)) {
    poly_divrem(&quotient, &remainder, &old_r, &r, ctx);

    poly_copy(&old_r, &r);
    poly_copy(&r, &remainder);

    /* === CORRECTION ICI === */
    /* s_new = old_s + quotient * s     (caractéristique 2) */
    poly_mul(&product, &quotient, &s, ctx);
    poly_add(&temp, &old_s, &product);
    poly_copy(&old_s, &s);
    poly_copy(&s, &temp);

    /* t_new = old_t + quotient * t */
    poly_mul(&product, &quotient, &t, ctx);
    poly_add(&temp, &old_t, &product);
    poly_copy(&old_t, &t);
    poly_copy(&t, &temp);
  }

  poly_copy(d, &old_r);
  poly_copy(u, &old_s);
  poly_copy(v, &old_t);

  /* libération */
  poly_free(&old_r);
  poly_free(&r);
  poly_free(&old_s);
  poly_free(&s);
  poly_free(&old_t);
  poly_free(&t);
  poly_free(&quotient);
  poly_free(&remainder);
  poly_free(&temp);
  poly_free(&product);
}

/* forward prototypes */
static void poly_xpow_mod(poly_t *res, uint64_t exp, const poly_t *modp, const gf_ctx_t *ctx);

/* ==================== RAND IRREDUCIBLE (avec vrai test) ==================== */
static int poly_has_root(const poly_t *p, const gf_ctx_t *ctx) {
  gf_t q = (gf_t)1 << ctx->m;
  for (gf_t a = 0; a < q; a++) {
    if (poly_eval(p, a, ctx) == 0) return 1;
  }
  return 0;
}

static int poly_is_irreducible(const poly_t *g, const gf_ctx_t *ctx) {
  if (g->deg <= 0 || g->coeff[g->deg] != 1) return 0;
  if (poly_has_root(g, ctx)) return 0;

  int t = g->deg;
  for (int d = 1; d <= t/2; d++) {
    uint64_t qd = (uint64_t)1 << (ctx->m * d);
    // X^qd mod g
    poly_t xqd = poly_alloc(g->deg - 1);
    // implémentation simple de X^qd (binary exp) → je te donne la fonction ci-dessous
    poly_xpow_mod(&xqd, qd, g, ctx);               // fonction ajoutée plus bas

    poly_t f = poly_alloc(g->deg - 1);
    poly_copy(&f, &xqd);
    if (f.deg >= 1)
      f.coeff[1] = gf_add(f.coeff[1], 1);
    else {
      f.coeff[1] = 1; f.deg = 1;
    }

    poly_t gcd_res = poly_alloc(g->deg);
    poly_gcd(&gcd_res, g, &f, ctx);
    int bad = (gcd_res.deg > 0);
    poly_free(&xqd); poly_free(&f); poly_free(&gcd_res);
    if (bad) return 0;
  }
  return 1;
}

int poly_rand_irreducible(poly_t *g, int t, const gf_ctx_t *ctx, int (*u8rnd)()) {
  if (!g || g->size < t+1) return -1;
  int tries = 0;
  do {
    for (int i = 0; i < t; i++) {
      gf_t c = 0;
      for (int j = 0; j < ctx->m; j++) {
        c |= ((gf_t)(u8rnd() & 1) << j);
      }
      g->coeff[i] = c;
    }
    g->coeff[t] = 1;
    g->deg = t;
    tries++;
  } while (!poly_is_irreducible(g, ctx) && tries < 1000);
  return (tries < 1000) ? 0 : -1;
}

/* ==================== X^exp mod g (utilisé pour le test irréductibilité) ==================== */
static void poly_xpow_mod(poly_t *res, uint64_t exp, const poly_t *modp, const gf_ctx_t *ctx) {
  poly_set_one(res);
  if (exp == 0) return;

  poly_t curr = poly_alloc(modp->deg - 1);
  curr.deg = 1; curr.coeff[1] = 1;   // X

  uint64_t e = exp;
  while (e) {
    if (e & 1) {
      poly_t prod = poly_alloc(modp->deg - 1);
      poly_mul(&prod, res, &curr, ctx);
      poly_mod(res, &prod, modp, ctx);
      poly_free(&prod);
    }
    poly_t sq = poly_alloc(modp->deg - 1);
    poly_mul(&sq, &curr, &curr, ctx);
    poly_mod(&curr, &sq, modp, ctx);
    poly_free(&sq);
    e >>= 1;
  }
  poly_free(&curr);   // curr libérée à la fin
}

/* ==================== RACINE CARRÉE MOD g (Alg. 5) ==================== */
void poly_sqrt_mod(poly_t *theta, const poly_t *Q, const poly_t *g, const gf_ctx_t *ctx) {
  if (!theta || !Q || !g || !ctx) return;
  poly_copy(theta, Q);
  if (theta->deg < 0) return;

  int mt = ctx->m * g->deg;
  for (int i = 0; i < mt - 1; i++) {          // mt-1 carrés
    poly_t sq = poly_alloc(g->deg - 1);
    poly_mul(&sq, theta, theta, ctx);
    poly_mod(theta, &sq, g, ctx);
    poly_free(&sq);
  }
}

/* ==================== ÉQUATION CLÉ (Alg. 6) ==================== */
void poly_solve_key_equation(poly_t *sigma, const poly_t *theta, const poly_t *g, const gf_ctx_t *ctx) {
  if (!sigma || !theta || !g || !ctx) return;

  int t_half = g->deg / 2;

  poly_t r0 = poly_alloc(g->deg); poly_copy(&r0, g);
  poly_t r1 = poly_alloc(g->deg); poly_copy(&r1, theta);
  poly_t q   = poly_alloc(g->deg);
  poly_t rem = poly_alloc(g->deg);

  while (r1.deg > t_half) {
    poly_divrem(&q, &rem, &r0, &r1, ctx);
    poly_copy(&r0, &r1);
    poly_copy(&r1, &rem);
  }
  poly_copy(sigma, &r1);

  poly_free(&r0); poly_free(&r1); poly_free(&q); poly_free(&rem);
}

/* ==================== RECHERCHE DE RACINES (Alg. 7) ==================== */
int poly_root_finding(const poly_t *f, gf_t *roots, size_t max_roots, const gf_ctx_t *ctx, int (*u8rnd)()) {
  (void)u8rnd;
  if (!f || !roots) return 0;

  size_t count = 0;
  gf_t q = (gf_t)1 << ctx->m;
  for (gf_t a = 0; a < q && count < max_roots; a++) {
    if (poly_eval(f, a, ctx) == 0) {
      roots[count++] = a;
    }
  }
  return (int)count;
}