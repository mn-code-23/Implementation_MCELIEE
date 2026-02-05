#include <stdlib.h>
#include <string.h>
#include "matrix.h"

/* ===== Allocation/Liberation ===== */

binmat_t binmat_alloc(size_t rows, size_t cols) {
  binmat_t m;
  m.rows = rows;
  m.cols = cols;
  m.words_per_row = (cols + 63) / 64;
  size_t total_words = rows * m.words_per_row;
  m.data = (uint64_t *)calloc(total_words, sizeof(uint64_t));
  return m;
}

void binmat_free(binmat_t *m) {
  if (!m) {
    return;
  }
  if (m->data) {
    free(m->data);
    m->data = NULL;
  }
  m->rows = 0;
  m->cols = 0;
  m->words_per_row = 0;
}

/* ===== Acces aux bits ===== */

int binmat_get(const binmat_t *m, size_t r, size_t c) {
  if (!m || !m->data || r >= m->rows || c >= m->cols) {
    return 0;
  }
  size_t word_idx = r * m->words_per_row + (c / 64);
  size_t bit_idx = c % 64;
  return (int)((m->data[word_idx] >> bit_idx) & 1);
}

void binmat_set(binmat_t *m, size_t r, size_t c, int value) {
  if (!m || !m->data || r >= m->rows || c >= m->cols) {
    return;
  }
  size_t word_idx = r * m->words_per_row + (c / 64);
  size_t bit_idx = c % 64;
  uint64_t mask = 1ULL << bit_idx;
  if (value) {
    m->data[word_idx] |= mask;
  } else {
    m->data[word_idx] &= ~mask;
  }
}

/* ===== Operations sur les lignes ===== */

void binmat_row_xor(binmat_t *m, size_t dst, size_t src) {
  if (!m || !m->data || dst >= m->rows || src >= m->rows) {
    return;
  }
  size_t src_offset = src * m->words_per_row;
  size_t dst_offset = dst * m->words_per_row;
  for (size_t w = 0; w < m->words_per_row; w++) {
    m->data[dst_offset + w] ^= m->data[src_offset + w];
  }
}

void binmat_row_swap(binmat_t *m, size_t r1, size_t r2) {
  if (!m || !m->data || r1 >= m->rows || r2 >= m->rows || r1 == r2) {
    return;
  }
  size_t off1 = r1 * m->words_per_row;
  size_t off2 = r2 * m->words_per_row;
  for (size_t w = 0; w < m->words_per_row; w++) {
    uint64_t tmp = m->data[off1 + w];
    m->data[off1 + w] = m->data[off2 + w];
    m->data[off2 + w] = tmp;
  }
}

/* ===== Operations sur les colonnes ===== */

void binmat_col_swap(binmat_t *m, size_t c1, size_t c2) {
  if (!m || !m->data || c1 >= m->cols || c2 >= m->cols || c1 == c2) {
    return;
  }
  for (size_t r = 0; r < m->rows; r++) {
    int b1 = binmat_get(m, r, c1);
    int b2 = binmat_get(m, r, c2);
    binmat_set(m, r, c1, b2);
    binmat_set(m, r, c2, b1);
  }
}

/* ===== Copie de bloc ===== */

void binmat_copy_block(
  binmat_t *dst,
  size_t dst_r,
  size_t dst_c,
  const binmat_t *src,
  size_t src_r,
  size_t src_c,
  size_t rows,
  size_t cols
) {
  if (!dst || !src || !dst->data || !src->data) {
    return;
  }
  if (dst_r + rows > dst->rows || src_r + rows > src->rows ||
      dst_c + cols > dst->cols || src_c + cols > src->cols) {
    return;
  }

  for (size_t r = 0; r < rows; r++) {
    for (size_t c = 0; c < cols; c++) {
      int bit = binmat_get(src, src_r + r, src_c + c);
      binmat_set(dst, dst_r + r, dst_c + c, bit);
    }
  }
}

/* ===== Multiplication matricielle (dans F2) ===== */

void binmat_mul(binmat_t *c, const binmat_t *a, const binmat_t *b) {
  if (!c || !a || !b || !c->data || !a->data || !b->data) {
    return;
  }
  if (a->cols != b->rows || a->rows != c->rows || b->cols != c->cols) {
    return;
  }

  /* Initialiser C a 0 */
  memset(c->data, 0, c->rows * c->words_per_row * sizeof(uint64_t));

  /* Multiplication C = A * B dans F2 */
  for (size_t i = 0; i < a->rows; i++) {
    for (size_t j = 0; j < b->cols; j++) {
      int sum = 0;
      for (size_t k = 0; k < a->cols; k++) {
        int a_bit = binmat_get(a, i, k);
        int b_bit = binmat_get(b, k, j);
        sum ^= (a_bit & b_bit);
      }
      binmat_set(c, i, j, sum);
    }
  }
}

/* ===== Multiplication vecteur-matrice (dans F2) ===== */

void binvec_mul_mat(uint8_t *y, const uint8_t *x, const binmat_t *m) {
  if (!y || !x || !m || !m->data) {
    return;
  }

  /* y = x * M, ou x est un vecteur ligne de m->rows bits */
  /* et y est un vecteur ligne de m->cols bits */

  for (size_t j = 0; j < m->cols; j++) {
    int sum = 0;
    for (size_t i = 0; i < m->rows; i++) {
      int x_bit = (x[i / 8] >> (i % 8)) & 1;
      int m_bit = binmat_get(m, i, j);
      sum ^= (x_bit & m_bit);
    }
    y[j / 8] = (uint8_t)((y[j / 8] & ~(1u << (j % 8))) | (sum << (j % 8)));
  }
}
