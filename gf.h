#ifndef GF_H
#define GF_H

#include <stdint.h>
#include "param.h"

/*
 * Representation d'un element de F2^m dans un entier.
 * Chaque bit represente un coefficient b_i de Y^i.
 */
typedef uint16_t gf_t;

/*
 * Contexte du corps F2^m : tables log/exp facultatives.
 */
typedef struct {
  int m;              /* degre du corps */
  int q;              /* q = 2^m */
  int q_minus_1;      /* q-1 */
  uint32_t prim_poly; /* polynome irreductible */
  gf_t *exp;          /* exp[i] = alpha^i */
  gf_t *log;          /* log[x] = i tel que x = alpha^i */
} gf_ctx_t;

/*
 * Initialise le contexte GF(2^m) (tables, reduction, etc.).
 * Operation mathematique : construction du corps F2[Y]/(f(Y)).
 */
int gf_init(gf_ctx_t *ctx, const mc_params_t *p);

/* Libere les tables allouees par gf_init. */
void gf_clear(gf_ctx_t *ctx);

/*
 * Addition dans F2^m : a + b = a XOR b.
 */
static inline gf_t gf_add(gf_t a, gf_t b) { return (gf_t)(a ^ b); }

/*
 * Multiplication dans F2^m.
 * Operation mathematique : produit modulo f(Y).
 */
gf_t gf_mul(const gf_ctx_t *ctx, gf_t a, gf_t b);

/*
 * Carre dans F2^m.
 * Operation mathematique : a^2 modulo f(Y).
 */
gf_t gf_square(const gf_ctx_t *ctx, gf_t a);

/*
 * Inversion dans F2^m.
 * Operation mathematique : a^{-1} tel que a * a^{-1} = 1.
 */
gf_t gf_inv(const gf_ctx_t *ctx, gf_t a);

/*
 * Racine carree dans F2^m.
 * Operation mathematique : sqrt(a) tel que (sqrt(a))^2 = a.
 */
gf_t gf_sqrt(const gf_ctx_t *ctx, gf_t a);

/*
 * Puissance a^e dans F2^m.
 */
gf_t gf_pow(const gf_ctx_t *ctx, gf_t a, int e);

#endif /* GF_H */
