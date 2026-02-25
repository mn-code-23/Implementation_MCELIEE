#include <stdlib.h>
#include <string.h>
#include "decrypt.h"
#include "util.h"

/* ===== Syndrome (Algorithme 4, etape 1) ===== */

void compute_syndrome(poly_t *Sc, const uint8_t *c, const gf_t *L, size_t n, const poly_t *g, const gf_ctx_t *ctx) {
  if (!Sc || !c || !L || !g || !ctx) {
    return;
  }
  if (!Sc->coeff || !g->coeff || !ctx->exp || !ctx->log) {
    return;
  }

  /*
   * Syndrome : S_c(X) = sum_{i=0}^{n-1} c_i * (X - L_i)^{-1} mod g(X)
   * 
   * Calcul iteratif : on accumule les termes dans Sc.
   */

  memset(Sc->coeff, 0, (size_t)Sc->size * sizeof(gf_t));
  Sc->deg = -1;

  for (size_t i = 0; i < n; i++) {
    int c_bit = (c[i / 8] >> (i % 8)) & 1;
    if (!c_bit) {
      continue; /* Sauter si c_i = 0 */
    }

    /* Construire le polynome (X - L_i)^{-1} mod g(X) */
    /* (X - L_i)^{-1} = 1 / (X - L_i) */
    /* On utilise des calculs dans F2^m[X] / g(X) */

    gf_t L_i = L[i];

    /* Construire le polynome (X - L_i) = X - L_i (en binaire, c'est X + L_i) */
    poly_t X_minus_Li = poly_alloc(1);
    X_minus_Li.coeff[0] = L_i;
    X_minus_Li.coeff[1] = 1; /* Coefficient de X */
    X_minus_Li.deg = 1;

    /* Inverser (X - L_i) modulo g(X) */
    poly_t inv = poly_alloc((size_t)g->deg);
    poly_inv_mod(&inv, &X_minus_Li, g, ctx);

    /* Ajouter inv a Sc : Sc += inv */
    poly_add(Sc, Sc, &inv);
    poly_mod(Sc, Sc, g, ctx);

    poly_free(&X_minus_Li);
    poly_free(&inv);
  }

  poly_normalize(Sc);
}

/* ===== Pre-calcul des inverses ===== */

void precompute_inverses(poly_t *inv_table, const gf_t *L, size_t n, const poly_t *g, const gf_ctx_t *ctx) {
  if (!inv_table || !L || !g || !ctx) {
    return;
  }

  for (size_t i = 0; i < n; i++) {
    gf_t L_i = L[i];

    /* Construire (X - L_i) */
    poly_t X_minus_Li = poly_alloc(1);
    X_minus_Li.coeff[0] = L_i;
    X_minus_Li.coeff[1] = 1;
    X_minus_Li.deg = 1;

    /* Inverser modulo g */
    poly_inv_mod(&inv_table[i], &X_minus_Li, g, ctx);

    poly_free(&X_minus_Li);
  }
}

/* ===== Inverse d'un polynome modulo g ===== */

void poly_inv_mod(poly_t *inv, const poly_t *a, const poly_t *g, const gf_ctx_t *ctx) {
  if (!inv || !a || !g || !ctx) {
    return;
  }

  /*
   * Utiliser l'algorithme d'Euclide etendu :
   * XGCD(a, g) -> d, u, v
   * u * a + v * g = d = 1
   * u * a = 1 (mod g)
   * Donc inv = u
   */

  poly_t d = poly_alloc((size_t)g->deg);
  poly_t u = poly_alloc((size_t)(a->deg + g->deg));
  poly_t v = poly_alloc((size_t)(a->deg + g->deg));

  poly_xgcd(&d, &u, &v, a, g, ctx);

  /* Copier u dans inv */
  memcpy(inv->coeff, u.coeff, (size_t)u.size * sizeof(gf_t));
  inv->deg = u.deg;

  poly_free(&d);
  poly_free(&u);
  poly_free(&v);
}

/* ===== Racine carree modulo g (Algorithme 5) ===== */

void patterson_sqrt_step(poly_t *theta, const poly_t *T, const poly_t *g, const gf_ctx_t *ctx) {
  if (!theta || !T || !g || !ctx) {
    return;
  }

  /* Utiliser poly_sqrt_mod qui est deja implementee */
  poly_sqrt_mod(theta, T, g, ctx);
}

/* ===== Equation cle (Algorithme 6) ===== */

void patterson_key_equation(poly_t *sigma, const poly_t *theta, const poly_t *g, const gf_ctx_t *ctx) {
  if (!sigma || !theta || !g || !ctx) {
    return;
  }

  /* Utiliser poly_solve_key_equation qui est deja implementee */
  poly_solve_key_equation(sigma, theta, g, ctx);
}

/* ===== Recherche de racines (Algorithme 7) ===== */

int patterson_root_finding(const poly_t *sigma, gf_t *roots, size_t max_roots, const gf_ctx_t *ctx, int (*u8rnd)()) {
  if (!sigma || !roots || !ctx || !u8rnd) {
    return 0;
  }

  (void)u8rnd; /* Non utilise pour la recherche exhaustive */

  /* Utiliser poly_root_finding qui est deja implementee */
  return poly_root_finding(sigma, roots, max_roots, ctx, u8rnd);
}

/* ===== Construction du vecteur d'erreur ===== */

void build_error_vector(uint8_t *e, const gf_t *L, size_t n, const gf_t *roots, size_t root_count) {
  if (!e || !L || !roots) {
    return;
  }

  /* Initialiser e a 0 */
  memset(e, 0, (n + 7) / 8);

  /* Pour chaque racine, trouver sa position dans L et fixer le bit */
  for (size_t i = 0; i < root_count; i++) {
    gf_t root = roots[i];
    for (size_t j = 0; j < n; j++) {
      if (L[j] == root) {
        int bit_idx = j % 8;
        int byte_idx = j / 8;
        e[byte_idx] |= (1u << bit_idx);
        break;
      }
    }
  }
}

/* ===== Decodage Patterson (Algorithmes 4-7) ===== */

int mc_decode(uint8_t *y, const uint8_t *c, const mc_secret_key_t *sk, const gf_ctx_t *ctx, int (*u8rnd)()) {
  if (!y || !c || !sk || !ctx || !u8rnd) {
    return -1;
  }
  if (!sk->L || !sk->g.coeff || !ctx->exp || !ctx->log) {
    return -1;
  }

  /*
   * Algorithme 4 : Decodage d'un code de Goppa binaire via Patterson
   */

  size_t n = (size_t)sk->params.n;
  int t = sk->params.t;

  /* 1. Calculer le syndrome S_c(X) */
  poly_t Sc = poly_alloc(t);
  if (!Sc.coeff) {
    return -1;
  }
  compute_syndrome(&Sc, c, sk->L, n, &sk->g, ctx);

  /* 2. Test codeword : si S_c = 0, alors y = c */
  if (Sc.deg < 0 || (Sc.deg == 0 && Sc.coeff[0] == 0)) {
    memcpy(y, c, (n + 7) / 8);
    poly_free(&Sc);
    return 0;
  }

  /* 3. Inverser le syndrome : T = S_c^{-1} mod g */
  poly_t T = poly_alloc(t);
  if (!T.coeff) {
    poly_free(&Sc);
    return -1;
  }
  poly_inv_mod(&T, &Sc, &sk->g, ctx);

  /* 4. Racine carree : theta = sqrt(T + X) mod g */
  poly_t T_plus_X = poly_alloc(t);
  if (!T_plus_X.coeff) {
    poly_free(&T);
    poly_free(&Sc);
    return -1;
  }
  /* Dans mc_decode, construction de T + X */
  poly_copy(&T_plus_X, &T);                     // au lieu du poly_add erroné
  if (T_plus_X.deg < 1) {
    T_plus_X.deg = 1;
    T_plus_X.coeff[1] = 1;
  } else {
    T_plus_X.coeff[1] = gf_add(T_plus_X.coeff[1], 1);
  }

  poly_t theta = poly_alloc(t);
  if (!theta.coeff) {
    poly_free(&T_plus_X);
    poly_free(&T);
    poly_free(&Sc);
    return -1;
  }
  patterson_sqrt_step(&theta, &T_plus_X, &sk->g, ctx);

  /* 5. Resolution de l'equation cle : sigma = ... */
  poly_t sigma = poly_alloc(t);
  if (!sigma.coeff) {
    poly_free(&theta);
    poly_free(&T_plus_X);
    poly_free(&T);
    poly_free(&Sc);
    return -1;
  }
  patterson_key_equation(&sigma, &theta, &sk->g, ctx);

  /* 6. Recherche de racines */
  gf_t roots[256];
  int root_count = patterson_root_finding(&sigma, roots, 256, ctx, u8rnd);
  
  if (root_count < 0 || root_count > t) {
    poly_free(&sigma);
    poly_free(&theta);
    poly_free(&T_plus_X);
    poly_free(&T);
    poly_free(&Sc);
    return -1;
  }

  /* 7. Construction du vecteur d'erreur */
  uint8_t *e = (uint8_t *)calloc((n + 7) / 8, sizeof(uint8_t));
  if (!e) {
    poly_free(&sigma);
    poly_free(&theta);
    poly_free(&T_plus_X);
    poly_free(&T);
    poly_free(&Sc);
    return -1;
  }
  build_error_vector(e, sk->L, n, roots, (size_t)root_count);

  /* 8. Sortie : y = c XOR e */
  for (size_t i = 0; i < (n + 7) / 8; i++) {
    y[i] = c[i] ^ e[i];
  }

  free(e);
  poly_free(&sigma);
  poly_free(&theta);
  poly_free(&T_plus_X);
  poly_free(&T);
  poly_free(&Sc);

  return 0;
}

/* ===== Dechiffrement (Algorithme 3) ===== */

int mc_decrypt(uint8_t *m, const uint8_t *c, const mc_secret_key_t *sk, const gf_ctx_t *ctx, int (*u8rnd)()) {
  if (!m || !c || !sk || !ctx || !u8rnd) {
    return -1;
  }
  if (!sk->L || !sk->g.coeff) {
    return -1;
  }

  /*
   * Algorithme 3 : Dechiffrement
   * 1. Decoder c pour obtenir y = c XOR e via Patterson.
   * 2. Recuperer m directement de y (sans S, sans P).
   *    m = y[0..k-1]
   */

  size_t k = (size_t)sk->params.k;
  size_t n = (size_t)sk->params.n;

  /* 1. Decoder */
  uint8_t *y = (uint8_t *)malloc((n + 7) / 8);
  if (!y) {
    return -1;
  }

  if (mc_decode(y, c, sk, ctx, u8rnd) != 0) {
    free(y);
    return -1;
  }

  /* 2. Recuperer m = y[0..k-1] */
  memcpy(m, y, (k + 7) / 8);
  
  /* Masquer les bits au-dela de k */
  int trailing_bits = k % 8;
  if (trailing_bits > 0) {
    int last_byte_idx = (k - 1) / 8;
    m[last_byte_idx] &= (1u << trailing_bits) - 1;
  }

  free(y);
  return 0;
}
