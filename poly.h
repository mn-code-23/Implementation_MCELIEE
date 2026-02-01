#ifndef POLY_H
#define POLY_H

#include <stddef.h>
#include "gf.h"

/*
 * Polynome sur F2^m : p(X) = sum_{i=0..deg} coeff[i] * X^i.
 */
typedef struct {
  int deg;      /* degre logique */
  int size;     /* taille du tableau coeff */
  gf_t *coeff;  /* coefficients dans F2^m */
} poly_t;

/* Allocation / liberation */
poly_t poly_alloc(int max_deg);
void poly_free(poly_t *p);

/*
 * Normalisation : ajuste le degre en ignorant les coefficients nuls.
 * Operation mathematique : deg(p) = max{i | coeff[i] != 0}.
 */
void poly_normalize(poly_t *p);

/*
 * Addition : r = a + b (dans F2^m[X]).
 */
void poly_add(poly_t *r, const poly_t *a, const poly_t *b);

/*
 * Multiplication : r = a * b (dans F2^m[X]).
 */
void poly_mul(poly_t *r, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx);

/*
 * Division euclidienne : a = q * b + r, deg(r) < deg(b).
 */
void poly_divrem(poly_t *q, poly_t *r, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx);

/*
 * Reste modulo g : r = a mod g.
 */
void poly_mod(poly_t *r, const poly_t *a, const poly_t *g, const gf_ctx_t *ctx);

/*
 * PGCD : d = gcd(a, b) dans F2^m[X].
 */
void poly_gcd(poly_t *d, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx);

/*
 * XGCD : u*a + v*b = d = gcd(a,b).
 * Operation mathematique : algorithme d'Euclide etendu.
 */
void poly_xgcd(poly_t *d, poly_t *u, poly_t *v, const poly_t *a, const poly_t *b, const gf_ctx_t *ctx);

/*
 * Evaluation : p(alpha) dans F2^m.
 */
gf_t poly_eval(const poly_t *p, gf_t alpha, const gf_ctx_t *ctx);

/*
 * Generation d'un polynome irreductible de degre t (Goppa).
 */
int poly_rand_irreducible(poly_t *g, int t, const gf_ctx_t *ctx, int (*u8rnd)());

/*
 * Racine carree modulo g (Alg. 5 du guide).
 * Operation mathematique : sqrt(Q(X)) mod g(X).
 */
void poly_sqrt_mod(poly_t *theta, const poly_t *Q, const poly_t *g, const gf_ctx_t *ctx);

/*
 * Resolution de l'equation cle (Alg. 6 du guide).
 * Entree : theta(X) ; sortie : sigma(X).
 * Operation mathematique : sigma = gamma^2 + X * phi^2.
 */
void poly_solve_key_equation(poly_t *sigma, const poly_t *theta, const poly_t *g, const gf_ctx_t *ctx);

/*
 * Recherche des racines dans F2^m (Alg. 7 du guide).
 * Retourne le nombre de racines, et remplit roots[].
 */
int poly_root_finding(const poly_t *f, gf_t *roots, size_t max_roots, const gf_ctx_t *ctx, int (*u8rnd)());

#endif /* POLY_H */
