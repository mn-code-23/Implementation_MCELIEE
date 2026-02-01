#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>
#include <stdint.h>

/*
 * Matrice binaire stockee en bitset.
 * - rows, cols : dimensions en bits.
 * - words_per_row : nombre de mots 64-bits par ligne.
 * - data : buffer contigu de size rows * words_per_row.
 */
typedef struct {
  size_t rows;
  size_t cols;
  size_t words_per_row;
  uint64_t *data;
} binmat_t;

/*
 * Alloue une matrice binaire rows x cols et initialise a 0.
 * Operation mathematique : creation de M in F2^{rows x cols}.
 */
binmat_t binmat_alloc(size_t rows, size_t cols);

/* Libere la memoire associee a la matrice. */
void binmat_free(binmat_t *m);

/*
 * Acces a un bit : retourne M[r,c] in {0,1}.
 */
int binmat_get(const binmat_t *m, size_t r, size_t c);

/*
 * Ecriture d'un bit : fixe M[r,c] <- value (0 ou 1).
 */
void binmat_set(binmat_t *m, size_t r, size_t c, int value);

/*
 * XOR de lignes : L_dst <- L_dst XOR L_src.
 * Operation mathematique : addition de lignes dans F2.
 */
void binmat_row_xor(binmat_t *m, size_t dst, size_t src);

/*
 * Echange de lignes : swap(L_r1, L_r2).
 */
void binmat_row_swap(binmat_t *m, size_t r1, size_t r2);

/*
 * Echange de colonnes : swap(C_c1, C_c2).
 */
void binmat_col_swap(binmat_t *m, size_t c1, size_t c2);

/*
 * Copie d'un bloc : dst[dst_r..] <- src[src_r..] pour un bloc rows x cols.
 */
void binmat_copy_block(
  binmat_t *dst,
  size_t dst_r,
  size_t dst_c,
  const binmat_t *src,
  size_t src_r,
  size_t src_c,
  size_t rows,
  size_t cols
);

/*
 * Multiplication binaire : C = A * B (dans F2).
 * Dimensions : A (r x k), B (k x n), C (r x n).
 */
void binmat_mul(binmat_t *c, const binmat_t *a, const binmat_t *b);

/*
 * Multiplication vecteur-matrice : y = x * M (dans F2).
 * x : vecteur binaire de taille M.rows, y : taille M.cols.
 */
void binvec_mul_mat(uint8_t *y, const uint8_t *x, const binmat_t *m);

#endif /* MATRIX_H */
