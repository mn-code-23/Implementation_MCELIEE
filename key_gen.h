#ifndef KEY_GEN_H
#define KEY_GEN_H

#include <stddef.h>
#include "param.h"
#include "gf.h"
#include "poly.h"
#include "matrix.h"

/*
 * Cle publique : G_pub et t.
 */
typedef struct {
  mc_params_t params;
  int t;
  binmat_t G_pub; /* k x n */
} mc_public_key_t;

/*
 * Cle secrete : support L, polynome g, matrices S et permutation P.
 */
typedef struct {
  mc_params_t params;
  gf_t *L;  /* support de taille n */
  poly_t g; /* polynome de Goppa degre t */
} mc_secret_key_t;

/*
 * Construction du support L (n elements distincts de F2^m).
 */
int build_support(gf_t *L, size_t n, const gf_ctx_t *ctx, int (*u8rnd)());

/*
 * Construire Y et Z puis H = Y * Z (Alg. 1, etapes 1-3).
 * - Y_{i,j} = alpha_j^i pour i=0..t-1
 * - Z = diag(g(alpha_j)^{-1})
 */
void build_parity_check(gf_t *H, const gf_t *L, size_t n, const poly_t *g, int t, const gf_ctx_t *ctx);

/*
 * Expansion binaire H' : remplace chaque coeff dans F2^m par un vecteur binaire de taille m.
 * Operation mathematique : H in F2^m -> H' in F2.
 */
void expand_parity_check(binmat_t *H_bin, const gf_t *H, size_t n, int t, int m);

/*
 * Gauss-Jordan mod 2 sur les lignes.
 * - Entree : M in F2^{r x n}.
 * - Sorties : M en RREF et pivots enregistres.
 */
size_t gauss_jordan_mod2(binmat_t *m, size_t *pivot_cols, size_t max_pivots);

/*
 * Systematisation : obtenir H'_r = [A | I_r].
 * - col_perm : permutation appliquee aux colonnes (taille n).
 *   Convention : col_perm[new_col] = old_col.
 * - Retour : 1 si succes (rang = r), 0 sinon.
 */
int systematize(binmat_t *h, size_t r, size_t n, size_t *col_perm);

/*
 * Construction de G = [I_k | A^T].
 */
void build_generator(const binmat_t *A, binmat_t *G);

/*
 * KeyGen complet (Alg. 1 du guide).
 * - Remplit pk et sk.
 */
int mc_keygen(const mc_params_t *params, mc_public_key_t *pk, mc_secret_key_t *sk, int (*u8rnd)());

#endif /* KEY_GEN_H */
