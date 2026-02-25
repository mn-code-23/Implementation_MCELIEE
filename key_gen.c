#include <stdlib.h>
#include <string.h>
#include "key_gen.h"

/* ===== Construction du support L ===== */

int build_support(gf_t *L, size_t n, const gf_ctx_t *ctx, int (*u8rnd)()) {
  if (!L || !ctx || !u8rnd) {
    return -1;
  }

  /* Remplir L avec n elements aleatoires distincts de F2^m */
  size_t filled = 0;
  while (filled < n) {
    gf_t elem = 0;
    for (int i = 0; i < ctx->m; i++) {
      int bit = u8rnd() & 1;
      elem |= (gf_t)(bit << i);
    }

    /* Verifier que elem n'est pas deja dans L */
    int found = 0;
    for (size_t j = 0; j < filled; j++) {
      if (L[j] == elem) {
        found = 1;
        break;
      }
    }
    if (!found) {
      L[filled] = elem;
      filled++;
    }
  }

  return 0;
}

/* ===== Construction de H = Y * Z ===== */

void build_parity_check(gf_t *H, const gf_t *L, size_t n, const poly_t *g, int t, const gf_ctx_t *ctx) {
  if (!H || !L || !g || !ctx || !g->coeff) {
    return;
  }

  /*
   * Construire la matrice de parite H de taille t x n.
   * H_{i,j} = (g(L_j)^{-1}) * L_j^i pour i = 0..t-1, j = 0..n-1
   * 
   * Cela correspond a Y * Z ou :
   * Y_{i,j} = L_j^i
   * Z = diag(g(L_j)^{-1})
   */

  for (size_t j = 0; j < n; j++) {
    gf_t g_inv = gf_inv(ctx, poly_eval(g, L[j], ctx));
    gf_t power = 1; /* L_j^0 = 1 */

    for (int i = 0; i < t; i++) {
      gf_t H_ij = gf_mul(ctx, power, g_inv);
      H[i * n + j] = H_ij;
      power = gf_mul(ctx, power, L[j]);
    }
  }
}

/* ===== Expansion binaire H -> H' ===== */

void expand_parity_check(binmat_t *H_bin, const gf_t *H, size_t n, int t, int m) {
  if (!H_bin || !H || H_bin->data == NULL) {
    return;
  }

  /*
   * Expansion : chaque element H_{i,j} dans F2^m devient m bits.
   * H_bin a des dimensions (m*t) x n.
   * Le bit de position (i*m + k, j) dans H_bin correspond au k-eme bit
   * de H[i*n + j] dans F2^m.
   */

  for (size_t j = 0; j < n; j++) {
    for (int i = 0; i < t; i++) {
      gf_t coeff = H[i * n + j];
      for (int k = 0; k < m; k++) {
        int bit = (coeff >> k) & 1;
        binmat_set(H_bin, (size_t)(i * m + k), j, bit);
      }
    }
  }
}

/* ===== Gauss-Jordan mod 2 ===== */

size_t gauss_jordan_mod2(binmat_t *m, size_t *pivot_cols, size_t max_pivots) {
  if (!m || !m->data || !pivot_cols) {
    return 0;
  }

  size_t pivot_count = 0;
  size_t current_row = 0;

  for (size_t col = 0; col < m->cols && current_row < m->rows && pivot_count < max_pivots; col++) {
    /* Trouver un pivot dans cette colonne (ligne >= current_row) */
    int pivot_row = -1;
    for (size_t row = current_row; row < m->rows; row++) {
      if (binmat_get(m, row, col) == 1) {
        pivot_row = (int)row;
        break;
      }
    }

    if (pivot_row == -1) {
      continue; /* Pas de pivot, passer a la colonne suivante */
    }

    /* Swap les lignes */
    if ((size_t)pivot_row != current_row) {
      binmat_row_swap(m, current_row, (size_t)pivot_row);
    }

    /* Enregistrer la colonne pivot */
    pivot_cols[pivot_count] = col;
    pivot_count++;

    /* Eliminer les autres 1s dans cette colonne */
    for (size_t row = 0; row < m->rows; row++) {
      if (row != current_row && binmat_get(m, row, col) == 1) {
        binmat_row_xor(m, row, current_row);
      }
    }

    current_row++;
  }

  return pivot_count;
}

/* ===== Systematisation : H' -> H'_r = [A | I_r] ===== */

int systematize(binmat_t *h, size_t r, size_t n, size_t *col_perm) {
  if (!h || !col_perm) return 0;

  size_t pivot_cols_buf[512];
  size_t pivot_count = gauss_jordan_mod2(h, pivot_cols_buf, r);
  if (pivot_count < r) return 0;

  /* === COL_PERM : non-pivots en premier, pivots en dernier → [A | I] === */
  size_t k = n - r;
  size_t idx = 0;
  for (size_t c = 0; c < n; c++) {
    int is_pivot = 0;
    for (size_t p = 0; p < pivot_count; p++) {
      if (pivot_cols_buf[p] == c) { is_pivot = 1; break; }
    }
    if (!is_pivot) col_perm[idx++] = c;
  }
  for (size_t p = 0; p < pivot_count; p++) {
    col_perm[idx++] = pivot_cols_buf[p];
  }

  /* === APPLIQUER LA PERMUTATION DE COLONNES (H_bin devient [A | I]) === */
  binmat_t temp = binmat_alloc(h->rows, h->cols);
  for (size_t row = 0; row < h->rows; row++) {
    for (size_t j = 0; j < h->cols; j++) {
      int bit = binmat_get(h, row, col_perm[j]);
      binmat_set(&temp, row, j, bit);
    }
  }
  memcpy(h->data, temp.data, h->rows * h->words_per_row * sizeof(uint64_t));
  binmat_free(&temp);

  return 1;
}

/* ===== Construction de G = [I_k | A^T] ===== */

void build_generator(const binmat_t *A, binmat_t *G) {
  if (!A || !G || !A->data || !G->data) {
    return;
  }

  /*
   * Construire G de taille k x n a partir de A de taille (m*t) x k.
   * G = [I_k | A^T]
   * 
   * Dans ce code, on suppose que G est deja dimensionne correctement.
   * G.rows = k (dimension du code)
   * G.cols = n (longueur du code)
   */

  size_t k = A->cols;

  /* Remplir la partie I_k (identite) */
  for (size_t i = 0; i < k && i < G->rows; i++) {
    binmat_set(G, i, i, 1);
  }

  /* Remplir la partie A^T (transpose de A) */
  for (size_t i = 0; i < k && i < G->rows; i++) {
    for (size_t j = 0; j < A->rows && k + j < G->cols; j++) {
      int bit = binmat_get(A, j, i);
      binmat_set(G, i, k + j, bit);
    }
  }
}

/* ===== KeyGen complet (Algorithme 1) ===== */

int mc_keygen(const mc_params_t *params, mc_public_key_t *pk, mc_secret_key_t *sk, int (*u8rnd)()) {
  if (!params || !pk || !sk || !u8rnd) {
    return -1;
  }

  /* Initialiser le contexte GF */
  gf_ctx_t ctx = {0};
  if (gf_init(&ctx, params) != 0) {
    return -1;
  }

  /* 1. Construire le support L */
  gf_t *L = (gf_t *)malloc((size_t)params->n * sizeof(gf_t));
  if (!L) {
    gf_clear(&ctx);
    return -1;
  }

  if (build_support(L, (size_t)params->n, &ctx, u8rnd) != 0) {
    free(L);
    gf_clear(&ctx);
    return -1;
  }

  /* 2. Generer un polynome de Goppa irreductible */
  poly_t g = poly_alloc(params->t);
  if (!g.coeff) {
    free(L);
    gf_clear(&ctx);
    return -1;
  }

  if (poly_rand_irreducible(&g, params->t, &ctx, u8rnd) != 0) {
    poly_free(&g);
    free(L);
    gf_clear(&ctx);
    return -1;
  }

  /* 3. Construire la matrice de parite H */
  gf_t *H = (gf_t *)malloc((size_t)(params->t * params->n) * sizeof(gf_t));
  if (!H) {
    poly_free(&g);
    free(L);
    gf_clear(&ctx);
    return -1;
  }

  build_parity_check(H, L, (size_t)params->n, &g, params->t, &ctx);

  /* 4. Expansion binaire H -> H' */
  binmat_t H_bin = binmat_alloc((size_t)(params->m * params->t), (size_t)params->n);
  if (!H_bin.data) {
    free(H);
    poly_free(&g);
    free(L);
    gf_clear(&ctx);
    return -1;
  }

  expand_parity_check(&H_bin, H, (size_t)params->n, params->t, params->m);

  /* 5. Systematiser H' en [A | I_{m*t}] */
  size_t col_perm[512];
  memset(col_perm, 0, sizeof(col_perm));

  if (!systematize(&H_bin, (size_t)(params->m * params->t), (size_t)params->n, col_perm)) {
    binmat_free(&H_bin);
    free(H);
    poly_free(&g);
    free(L);
    gf_clear(&ctx);
    return -1;
  }

  /* === PERMUTER L SELON LA MÊME PERMUTATION === */
  gf_t *L_new = (gf_t *)malloc((size_t)params->n * sizeof(gf_t));
  if (!L_new) {
    binmat_free(&H_bin);
    free(H);
    poly_free(&g);
    free(L);
    gf_clear(&ctx);
    return -1;
  }
  for (size_t j = 0; j < (size_t)params->n; j++) {
    L_new[j] = L[col_perm[j]];
  }
  free(L);
  L = L_new;

  /* 6. Extraire A et construire G = [I_k | A^T] */
  /* A est H_bin[*, 0..k-1] */
  binmat_t A = binmat_alloc((size_t)(params->m * params->t), (size_t)params->k);
  if (!A.data) {
    binmat_free(&H_bin);
    free(H);
    poly_free(&g);
    free(L);
    gf_clear(&ctx);
    return -1;
  }

  binmat_copy_block(&A, 0, 0, &H_bin, 0, 0, (size_t)(params->m * params->t), (size_t)params->k);

  pk->G_pub = binmat_alloc((size_t)params->k, (size_t)params->n);
  if (!pk->G_pub.data) {
    binmat_free(&A);
    binmat_free(&H_bin);
    free(H);
    poly_free(&g);
    free(L);
    gf_clear(&ctx);
    return -1;
  }

  build_generator(&A, &pk->G_pub);

  /* 7. Remplir les structures de cles */
  pk->params = *params;
  pk->t = params->t;

  sk->params = *params;
  sk->L = L;
  sk->g = g;

  /* Cleanup */
  binmat_free(&A);
  binmat_free(&H_bin);
  free(H);
  gf_clear(&ctx);

  return 0;
}
