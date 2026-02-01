#ifndef DECRYPT_H
#define DECRYPT_H

#include <stddef.h>
#include <stdint.h>
#include "key_gen.h"
#include "poly.h"

/*
 * Patterson Decode (Alg. 4-7 du guide).
 * Entree : c in F2^n, sk = (L, g)
 * Sortie : y = c XOR e ou echec.
 */
int mc_decode(uint8_t *y, const uint8_t *c, const mc_secret_key_t *sk, const gf_ctx_t *ctx, int (*u8rnd)());

/*
 * Dechiffrement (Alg. 3 du guide).
 * - Decode puis recupere m directement (pas de S, pas de P).
 */
int mc_decrypt(uint8_t *m, const uint8_t *c, const mc_secret_key_t *sk, const gf_ctx_t *ctx, int (*u8rnd)());

/*
 * Syndrome S_c(X) = sum_i c_i * (X - L_i)^{-1} mod g(X).
 */
void compute_syndrome(poly_t *Sc, const uint8_t *c, const gf_t *L, size_t n, const poly_t *g, const gf_ctx_t *ctx);

/*
 * Pre-calcul des inverses (X - L_i)^{-1} mod g(X).
 */
void precompute_inverses(poly_t *inv_table, const gf_t *L, size_t n, const poly_t *g, const gf_ctx_t *ctx);

/*
 * Inverse d'un polynome modulo g(X).
 * Operation mathematique : a^{-1} mod g.
 */
void poly_inv_mod(poly_t *inv, const poly_t *a, const poly_t *g, const gf_ctx_t *ctx);

/*
 * Etape Alg. 5 : racine carree modulo g(X).
 */
void patterson_sqrt_step(poly_t *theta, const poly_t *T, const poly_t *g, const gf_ctx_t *ctx);

/*
 * Etape Alg. 6 : resolution de l'equation cle.
 */
void patterson_key_equation(poly_t *sigma, const poly_t *theta, const poly_t *g, const gf_ctx_t *ctx);

/*
 * Etape Alg. 7 : recherche de racines dans F2^m.
 * Retourne le nombre de racines.
 */
int patterson_root_finding(const poly_t *sigma, gf_t *roots, size_t max_roots, const gf_ctx_t *ctx, int (*u8rnd)());

/*
 * Construction du vecteur d'erreur a partir des racines.
 * e_i = 1 si L_i est racine.
 */
void build_error_vector(uint8_t *e, const gf_t *L, size_t n, const gf_t *roots, size_t root_count);

#endif /* DECRYPT_H */
