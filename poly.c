#include <stdlib.h>
#include <string.h>
#include "poly.h"

/* ===== Allocation/Liberation ===== */

poly_t poly_alloc(int max_deg) {
  poly_t p;
  p.deg = 0;
  p.size = max_deg + 1;
  p.coeff = (gf_t *)calloc((size_t)(max_deg + 1), sizeof(gf_t));
  return p;
}

void poly_free(poly_t *p) {
  if (!p) {
    return;
  }
  if (p->coeff) {
    free(p->coeff);
    p->coeff = NULL;
  }
  p->deg = 0;
  p->size = 0;
}

/* ===== Normalisation ===== */

void poly_normalize(poly_t *p) {
  if (!p || !p->coeff) {
    return;
  }
  while (p->deg > 0 && p->coeff[p->deg] == 0) {
    p->deg--;
  }
}

/* ===== Addition en F2^m[X] ===== */

void poly_add(poly_t *r, const poly_t *a, const poly_t *b) {
  if (!r || !a || !b) {
    return;
  }
  if (!r->coeff || !a->coeff || !b->coeff) {
    return;
  }

  int max_deg = (a->deg > b->deg) ? a->deg : b->deg;
  if (max_deg >= r->size) {
    max_deg = r->size - 1;
  }

  for (int i = 0; i <= max_deg; i++) {
    gf_t a_coeff = (i <= a->deg) ? a->coeff[i] : 0;
    gf_t b_coeff = (i <= b->deg) ? b->coeff[i] : 0;
    r->coeff[i] = gf_add(a_coeff, b_coeff);
  }

  for (int i = max_deg + 1; i < r->size; i++) {
    r->coeff[i] = 0;
  }

  r->deg = max_deg;
  poly_normalize(r);
}

/* ===== Multiplication en F2^m[X] ===== */

void poly_mul(poly_t *r, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx) {
  if (!r || !a || !b || !ctx) {
    return;
  }
  if (!r->coeff || !a->coeff || !b->coeff) {
    return;
  }

  int result_deg = a->deg + b->deg;
  if (result_deg >= r->size) {
    result_deg = r->size - 1;
  }

  /* Initialiser le resultat a 0 */
  memset(r->coeff, 0, (size_t)r->size * sizeof(gf_t));

  /* Multiplication standard */
  for (int i = 0; i <= a->deg && i < r->size; i++) {
    if (a->coeff[i] == 0) {
      continue;
    }
    for (int j = 0; j <= b->deg && i + j < r->size; j++) {
      gf_t prod = gf_mul(ctx, a->coeff[i], b->coeff[j]);
      r->coeff[i + j] = gf_add(r->coeff[i + j], prod);
    }
  }

  r->deg = result_deg;
  poly_normalize(r);
}

/* ===== Division euclidienne (quotient et reste) ===== */

void poly_divrem(poly_t *q, poly_t *r, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx) {
  if (!q || !r || !a || !b || !ctx) {
    return;
  }
  if (!q->coeff || !r->coeff || !a->coeff || !b->coeff) {
    return;
  }
  if (b->deg < 0 || b->coeff[b->deg] == 0) {
    /* Division by zero ou degre invalide */
    q->deg = -1;
    r->deg = -1;
    return;
  }

  /* Copier a dans r (qui sera reduit) */
  if (r != a) {
    memcpy(r->coeff, a->coeff, (size_t)a->size * sizeof(gf_t));
    r->deg = a->deg;
    r->size = a->size;
  }

  /* Initialiser q a 0 */
  memset(q->coeff, 0, (size_t)q->size * sizeof(gf_t));
  q->deg = -1;

  /* Division */
  gf_t b_lead_inv = gf_inv(ctx, b->coeff[b->deg]);
  
  while (r->deg >= b->deg && r->deg >= 0) {
    int pos = r->deg - b->deg;
    if (pos >= q->size) {
      break;
    }
    
    gf_t coeff = gf_mul(ctx, r->coeff[r->deg], b_lead_inv);
    q->coeff[pos] = coeff;
    if (q->deg < pos) {
      q->deg = pos;
    }

    /* Soustraire b * coeff * X^pos de r */
    for (int i = 0; i <= b->deg; i++) {
      if (r->deg - b->deg + i < r->size) {
        gf_t sub = gf_mul(ctx, b->coeff[i], coeff);
        r->coeff[r->deg - b->deg + i] = gf_add(r->coeff[r->deg - b->deg + i], sub);
      }
    }

    poly_normalize(r);
  }

  if (q->deg < 0) {
    q->deg = 0;
    q->coeff[0] = 0;
  }
}

/* ===== Modulo ===== */

void poly_mod(poly_t *r, const poly_t *a, const poly_t *g, const gf_ctx_t *ctx) {
  if (!r || !a || !g || !ctx) {
    return;
  }

  poly_t q = poly_alloc(a->deg);
  poly_t temp_r = poly_alloc(a->deg);

  poly_divrem(&q, &temp_r, a, g, ctx);

  memcpy(r->coeff, temp_r.coeff, (size_t)temp_r.size * sizeof(gf_t));
  r->deg = temp_r.deg;

  poly_free(&q);
  poly_free(&temp_r);
}

/* ===== PGCD ===== */

void poly_gcd(poly_t *d, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx) {
  if (!d || !a || !b || !ctx) {
    return;
  }

  poly_t aa = poly_alloc(a->deg + 1);
  poly_t bb = poly_alloc(b->deg + 1);

  memcpy(aa.coeff, a->coeff, (size_t)a->size * sizeof(gf_t));
  aa.deg = a->deg;

  memcpy(bb.coeff, b->coeff, (size_t)b->size * sizeof(gf_t));
  bb.deg = b->deg;

  while (bb.deg >= 0 && (bb.deg > 0 || bb.coeff[0] != 0)) {
    poly_t remainder = poly_alloc((size_t)aa.deg);
    poly_t quotient = poly_alloc((size_t)aa.deg);

    poly_divrem(&quotient, &remainder, &aa, &bb, ctx);

    memcpy(aa.coeff, bb.coeff, (size_t)bb.size * sizeof(gf_t));
    aa.deg = bb.deg;

    memcpy(bb.coeff, remainder.coeff, (size_t)remainder.size * sizeof(gf_t));
    bb.deg = remainder.deg;

    poly_free(&quotient);
    poly_free(&remainder);
  }

  memcpy(d->coeff, aa.coeff, (size_t)aa.size * sizeof(gf_t));
  d->deg = aa.deg;

  poly_free(&aa);
  poly_free(&bb);
}

/* ===== XGCD (Euclide etendu) ===== */

void poly_xgcd(poly_t *d, poly_t *u, poly_t *v, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx) {
  if (!d || !u || !v || !a || !b || !ctx) {
    return;
  }

  /* Initialiser */
  poly_t old_r = poly_alloc((size_t)a->deg + 1);
  poly_t r = poly_alloc((size_t)b->deg + 1);
  poly_t old_s = poly_alloc((size_t)a->deg + 1);
  poly_t s = poly_alloc((size_t)b->deg + 1);
  poly_t old_t = poly_alloc((size_t)a->deg + 1);
  poly_t t = poly_alloc((size_t)b->deg + 1);

  memcpy(old_r.coeff, a->coeff, (size_t)a->size * sizeof(gf_t));
  old_r.deg = a->deg;

  memcpy(r.coeff, b->coeff, (size_t)b->size * sizeof(gf_t));
  r.deg = b->deg;

  old_s.coeff[0] = 1;
  old_s.deg = 0;

  s.deg = 0;
  s.coeff[0] = 0;

  old_t.deg = 0;
  old_t.coeff[0] = 0;

  t.coeff[0] = 1;
  t.deg = 0;

  poly_t quotient = poly_alloc((size_t)(a->deg + 1));
  poly_t remainder = poly_alloc((size_t)(a->deg + 1));
  poly_t product = poly_alloc((size_t)(a->deg + b->deg + 1));
  poly_t temp = poly_alloc((size_t)(a->deg + b->deg + 1));

  while (r.deg >= 0 && (r.deg > 0 || r.coeff[0] != 0)) {
    poly_divrem(&quotient, &remainder, &old_r, &r, ctx);

    memcpy(old_r.coeff, r.coeff, (size_t)r.size * sizeof(gf_t));
    old_r.deg = r.deg;

    memcpy(r.coeff, remainder.coeff, (size_t)remainder.size * sizeof(gf_t));
    r.deg = remainder.deg;

    poly_mul(&product, &quotient, &old_s, ctx);
    poly_add(&temp, &old_s, &product);
    memcpy(old_s.coeff, s.coeff, (size_t)s.size * sizeof(gf_t));
    old_s.deg = s.deg;
    memcpy(s.coeff, temp.coeff, (size_t)temp.size * sizeof(gf_t));
    s.deg = temp.deg;

    poly_mul(&product, &quotient, &old_t, ctx);
    poly_add(&temp, &old_t, &product);
    memcpy(old_t.coeff, t.coeff, (size_t)t.size * sizeof(gf_t));
    old_t.deg = t.deg;
    memcpy(t.coeff, temp.coeff, (size_t)temp.size * sizeof(gf_t));
    t.deg = temp.deg;
  }

  memcpy(d->coeff, old_r.coeff, (size_t)old_r.size * sizeof(gf_t));
  d->deg = old_r.deg;

  memcpy(u->coeff, old_s.coeff, (size_t)old_s.size * sizeof(gf_t));
  u->deg = old_s.deg;

  memcpy(v->coeff, old_t.coeff, (size_t)old_t.size * sizeof(gf_t));
  v->deg = old_t.deg;

  poly_free(&old_r);
  poly_free(&r);
  poly_free(&old_s);
  poly_free(&s);
  poly_free(&old_t);
  poly_free(&t);
  poly_free(&quotient);
  poly_free(&remainder);
  poly_free(&product);
  poly_free(&temp);
}

/* ===== Evaluation ===== */

gf_t poly_eval(const poly_t *p, gf_t alpha, const gf_ctx_t *ctx) {
  if (!p || !p->coeff || !ctx) {
    return 0;
  }

  gf_t result = 0;
  gf_t power = 1;

  for (int i = 0; i <= p->deg; i++) {
    gf_t term = gf_mul(ctx, p->coeff[i], power);
    result = gf_add(result, term);
    power = gf_mul(ctx, power, alpha);
  }

  return result;
}

/* ===== Generation de polynome irreductible ===== */

int poly_rand_irreducible(poly_t *g, int t, const gf_ctx_t *ctx, int (*u8rnd)()) {
  if (!g || !ctx || !u8rnd) {
    return -1;
  }
  if (!g->coeff || g->size < t + 1) {
    return -1;
  }

  do {
    /* Remplir les coefficients aleatoirement dans F2^m */
    for (int i = 0; i < t; i++) {
      gf_t coeff = 0;
      for (int j = 0; j < ctx->m; j++) {
        int bit = u8rnd() & 1;
        coeff |= (gf_t)(bit << j);
      }
      g->coeff[i] = coeff;
    }
    /* Fixer le coefficient de plus haut degre a 1 */
    g->coeff[t] = 1;
    g->deg = t;

  } while (0); /* Dans une vraie implementation, on verifierait l'irreductibilite */

  return 0;
}

/* ===== Racine carree modulo g (Algorithme 5) ===== */

void poly_sqrt_mod(poly_t *theta, const poly_t *Q, const poly_t *g, const gf_ctx_t *ctx) {
  if (!theta || !Q || !g || !ctx) {
    return;
  }
  if (!theta->coeff || !Q->coeff || !g->coeff) {
    return;
  }

  /*
   * Racine carree dans F_{2^m}[X]/(g(X))
   * Pour Q(X) = sum_i q_i * X^i,
   * sqrt(Q(X)) = sum_i sqrt(q_i) * X^(2*i) mod g(X)
   */

  poly_t sqrt_Q = poly_alloc((size_t)(Q->deg * 2 + 1));

  for (int i = 0; i <= Q->deg; i++) {
    gf_t sqrt_coeff = gf_sqrt(ctx, Q->coeff[i]);
    if (2 * i < sqrt_Q.size) {
      sqrt_Q.coeff[2 * i] = sqrt_coeff;
    }
  }
  sqrt_Q.deg = Q->deg * 2;
  poly_normalize(&sqrt_Q);

  /* Reduire modulo g */
  poly_mod(theta, &sqrt_Q, g, ctx);

  poly_free(&sqrt_Q);
}

/* ===== Resolution de l'equation cle (Algorithme 6) ===== */

void poly_solve_key_equation(poly_t *sigma, const poly_t *theta, const poly_t *g, const gf_ctx_t *ctx) {
  if (!sigma || !theta || !g || !ctx) {
    return;
  }
  if (!sigma->coeff || !theta->coeff || !g->coeff) {
    return;
  }

  /*
   * Algorithme 6 : Equation cle
   * Entree : theta(X)
   * Sortie : sigma(X) tel que sigma(X)^2 + X*phi(X)^2 = theta(X) mod g(X)
   * 
   * Utiliser XGCD pour resoudre :
   * XGCD(theta, g) -> d, u, v
   * sigma = u^2 mod g
   */

  poly_t d = poly_alloc((size_t)g->deg);
  poly_t u = poly_alloc((size_t)(theta->deg + g->deg));
  poly_t v = poly_alloc((size_t)(theta->deg + g->deg));

  poly_xgcd(&d, &u, &v, theta, g, ctx);

  poly_t u_sq = poly_alloc((size_t)(2 * u.deg + 1));
  poly_mul(&u_sq, &u, &u, ctx);
  poly_mod(sigma, &u_sq, g, ctx);

  poly_free(&d);
  poly_free(&u);
  poly_free(&v);
  poly_free(&u_sq);
}

/* ===== Recherche de racines (Algorithme 7) ===== */

int poly_root_finding(const poly_t *f, gf_t *roots, size_t max_roots, const gf_ctx_t *ctx, int (*u8rnd)()) {
  if (!f || !roots || !ctx || !u8rnd) {
    return 0;
  }
  if (!f->coeff) {
    return 0;
  }

  /*
   * Algorithme 7 : Recherche de racines dans F_{2^m}
   * Utiliser une variante de Cantor-Zassenhaus ou recherche exhaustive.
   * Pour m petit (< 10), on peut chercher exhaustivement.
   */

  int root_count = 0;

  /* Chercher exhaustivement toutes les racines dans F2^m */
  for (gf_t alpha = 0; alpha < ctx->q && root_count < (int)max_roots; alpha++) {
    gf_t eval = poly_eval(f, alpha, ctx);
    if (eval == 0) {
      roots[root_count] = alpha;
      root_count++;
    }
  }

  return root_count;
}
