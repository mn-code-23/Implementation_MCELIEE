#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "param.h"
#include "gf.h"
#include "key_gen.h"
#include "encrypt.h"
#include "decrypt.h"
#include "rng.h"

static int u8rnd(void) {
  unsigned char b = 0;
  if (randombytes(&b, 1) != RNG_SUCCESS) {
    return -1;
  }
  return (int)b;
}

/* =================== tests unitaires =================== */

static int test_gf(const mc_params_t *params) {
  gf_ctx_t ctx = {0};
  if (gf_init(&ctx, params) != 0) return -1;
  for (gf_t a = 1; a < ctx.q; a++) {
    gf_t inv = gf_inv(&ctx, a);
    if (gf_mul(&ctx, a, inv) != 1) {
      gf_clear(&ctx);
      return -1;
    }
  }
  gf_clear(&ctx);
  return 0;
}

static int test_poly(const gf_ctx_t *ctx, int (*rnd)()) {
  /* construire deux polynomes p et d alatoires et verifier division */
  poly_t p = poly_alloc(5);
  poly_t d = poly_alloc(4);
  poly_t q = poly_alloc(6);
  poly_t r = poly_alloc(6);
  if (!p.coeff || !d.coeff) return -1;

  /* degre d >= 1 */
  d.deg = 1 + (rnd() % 3);
  for (int i = 0; i <= d.deg; i++) d.coeff[i] = (gf_t)(rnd() & ((1 << ctx->m) - 1));
  d.coeff[d.deg] = 1;

  p.deg = rnd() % 5;
  for (int i = 0; i <= p.deg; i++) p.coeff[i] = (gf_t)(rnd() & ((1 << ctx->m) - 1));

  poly_divrem(&q, &r, &p, &d, ctx);

  /* reconstruit p' = q*d + r */
  poly_t prod = poly_alloc(q.deg + d.deg);
  poly_mul(&prod, &q, &d, ctx);
  poly_add(&prod, &prod, &r);

  /* comparer coefficient par coefficient */
  int ok = 1;
  if (prod.deg != p.deg) ok = 0;
  else {
    for (int i = 0; i <= p.deg; i++) {
      if (prod.coeff[i] != p.coeff[i]) { ok = 0; break; }
    }
  }
  if (r.deg >= d.deg) ok = 0;

  poly_free(&p); poly_free(&d); poly_free(&q); poly_free(&r); poly_free(&prod);
  return ok ? 0 : -1;
}

static int test_system_and_orthog(const mc_params_t *params, const gf_ctx_t *ctx, int (*rnd)()) {
  /* recréer ce que fait mc_keygen pour H' et G puis tester */
  size_t n = params->n;
  int t = params->t;
  gf_t *L = malloc(n * sizeof(gf_t));
  if (!L) return -1;
  build_support(L, n, ctx, rnd);
  poly_t g = poly_alloc(t);
  poly_rand_irreducible(&g, t, ctx, rnd);

  gf_t *H = malloc((size_t)t * n * sizeof(gf_t));
  build_parity_check(H, L, n, &g, t, ctx);

  binmat_t H_bin = binmat_alloc((size_t)(ctx->m * t), n);
  expand_parity_check(&H_bin, H, n, t, ctx->m);

  size_t *col_perm = calloc(n, sizeof(size_t));
  if (!systematize(&H_bin, (size_t)(ctx->m * t), n, col_perm)) {
    free(L); free(H); free(col_perm); return -1;
  }

  /* test que le bloc droit est I */
  size_t r = (size_t)(ctx->m * t);
  size_t k = n - r;
  int ok = 1;
  for (size_t i = 0; i < r; i++) {
    for (size_t j = 0; j < r; j++) {
      int bit = binmat_get(&H_bin, i, k + j);
      if ((i == j && bit != 1) || (i != j && bit != 0)) {
        ok = 0;
        break;
      }
    }
    if (!ok) break;
  }
  if (!ok) {
    free(L); free(H); free(col_perm);
    binmat_free(&H_bin);
    return -1;
  }

  /* construire G_pub */
  binmat_t A = binmat_alloc(r, k);
  binmat_copy_block(&A, 0, 0, &H_bin, 0, 0, r, k);
  binmat_t G = binmat_alloc(k, n);
  build_generator(&A, &G);

  /* vérifier orthogonalité : H_bin * G^T == 0 (r x k matrix) */
  binmat_t Gt = binmat_alloc(n, k);
  for (size_t i = 0; i < k; i++) {
    for (size_t j = 0; j < n; j++) {
      binmat_set(&Gt, j, i, binmat_get(&G, i, j));
    }
  }
  binmat_t prod = binmat_alloc(r, k);
  binmat_mul(&prod, &H_bin, &Gt);
  for (size_t i = 0; i < r && ok; i++) {
    for (size_t j = 0; j < k; j++) {
      if (binmat_get(&prod, i, j) != 0) { ok = 0; break; }
    }
  }

  free(L); free(H); free(col_perm);
  binmat_free(&H_bin); binmat_free(&A); binmat_free(&G);
  binmat_free(&Gt); binmat_free(&prod);
  poly_free(&g);
  return ok ? 0 : -1;
}

int main(void) {
  setbuf(stdout, NULL);
  setbuf(stderr, NULL);
  
  printf("\n=== TEST MCELIECEMS 2025 ===\n\n");
  
  mc_params_t params = MC_PARAMS_TOY;
  printf("Parametres: m=%d, t=%d, n=%d, k=%d\n\n", params.m, params.t, params.n, params.k);
  
  gf_ctx_t gf = {0};
  mc_public_key_t pk = {0};
  mc_secret_key_t sk = {0};

  printf("[0/6] Tests GF, polynôme, systématisation, orthogonalité...\n");
  fflush(stdout);

  if (gf_init(&gf, &params) != 0) {
    fprintf(stderr, "ERREUR: gf_init failed\n");
    return 1;
  }
  if (test_gf(&params) != 0) {
    fprintf(stderr, "ECHEC: test_gf\n");
    return 1;
  }
  if (test_poly(&gf, u8rnd) != 0) {
    fprintf(stderr, "ECHEC: test_poly\n");
    return 1;
  }
  if (test_system_and_orthog(&params, &gf, u8rnd) != 0) {
    fprintf(stderr, "ECHEC: test_system_and_orthog\n");
    return 1;
  }
  printf("OK\n\n");
  fflush(stdout);

  printf("[1/4] Initialisation du corps de Galois F_{2^%d}...\n", params.m);
  fflush(stdout);

  printf("[2/4] Generation de la cle publique et privee...\n");
  fflush(stdout);
  clock_t t1 = clock();
  
  if (mc_keygen(&params, &pk, &sk, u8rnd) != 0) {
    fprintf(stderr, "ERREUR: mc_keygen failed\n");
    gf_clear(&gf);
    return 1;
  }
  
  clock_t t2 = clock();
  double keygen_time = (double)(t2 - t1) / CLOCKS_PER_SEC;
  printf("OK (%.2f sec)\n\n", keygen_time);
  fflush(stdout);

  uint8_t *m = (uint8_t *)calloc((size_t)params.k, sizeof(uint8_t));
  uint8_t *c = (uint8_t *)calloc((size_t)params.n, sizeof(uint8_t));
  uint8_t *m2 = (uint8_t *)calloc((size_t)params.k, sizeof(uint8_t));
  if (!m || !c || !m2) {
    fprintf(stderr, "ERREUR: allocation memory failed\n");
    free(m); free(c); free(m2);
    gf_clear(&gf);
    return 1;
  }

  for (size_t i = 0; i < (size_t)params.k; i++) {
    int r = u8rnd();
    m[i] = (r < 0) ? 0 : (uint8_t)(r & 1);
  }
  
  printf("[3/4] Chiffrement d'un message de %d bits...\n", params.k);
  fflush(stdout);
  t1 = clock();
  
  if (mc_encrypt(c, m, &pk, u8rnd) != 0) {
    fprintf(stderr, "ERREUR: mc_encrypt failed\n");
    free(m); free(c); free(m2);
    gf_clear(&gf);
    return 1;
  }
  
  t2 = clock();
  double enc_time = (double)(t2 - t1) / CLOCKS_PER_SEC;
  printf("OK (%.2f sec)\n\n", enc_time);
  fflush(stdout);

  printf("[4/4] Dechiffrement et verification...\n");
  fflush(stdout);
  t1 = clock();
  
  if (mc_decrypt(m2, c, &sk, &gf, u8rnd) != 0) {
    fprintf(stderr, "ERREUR: mc_decrypt failed\n");
    free(m); free(c); free(m2);
    gf_clear(&gf);
    return 1;
  }
  
  t2 = clock();
  double dec_time = (double)(t2 - t1) / CLOCKS_PER_SEC;
  
  if (memcmp(m, m2, (size_t)params.k) != 0) {
    fprintf(stderr, " ECHEC: decrypt(encrypt(m)) != m\n");
    free(m); free(c); free(m2);
    gf_clear(&gf);
    return 1;
  }
  
  printf(" OK (%.2f sec)\n\n", dec_time);
  fflush(stdout);

  printf("=== RESULTAT: SUCCESS ===\n");
  printf("Temps total: %.2f sec\n\n", keygen_time + enc_time + dec_time);
  fflush(stdout);

  free(m);
  free(c);
  free(m2);
  gf_clear(&gf);
  
  return 0;
}
