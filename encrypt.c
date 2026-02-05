#include <stdlib.h>
#include <string.h>
#include "encrypt.h"
#include "util.h"

/* ===== Generation du vecteur d'erreur ===== */

int generate_error_vector(uint8_t *e, size_t n, int t, int (*u8rnd)()) {
  if (!e || !u8rnd || t < 0 || t > (int)n) {
    return -1;
  }

  /* Initialiser e a 0 */
  memset(e, 0, (n + 7) / 8);

  /* Generer t positions aleatoires distinctes et fixer les bits correspondants */
  int count = 0;
  while (count < t) {
    size_t pos = 0;
    for (int i = 0; i < (int)sizeof(size_t); i++) {
      int byte = u8rnd();
      if (byte < 0) {
        byte = 0;
      }
      pos = (pos << 8) | (byte & 0xFF);
    }
    pos = pos % n;

    /* Verifier que ce bit n'est pas deja fixe */
    int bit_idx = pos % 8;
    int byte_idx = pos / 8;
    if ((e[byte_idx] & (1u << bit_idx)) == 0) {
      e[byte_idx] |= (1u << bit_idx);
      count++;
    }
  }

  return 0;
}

/* ===== Chiffrement (Algorithme 2) ===== */

int mc_encrypt(uint8_t *c, const uint8_t *m, const mc_public_key_t *pk, int (*u8rnd)()) {
  if (!c || !m || !pk || !u8rnd) {
    return -1;
  }
  if (!pk->G_pub.data) {
    return -1;
  }

  /*
   * Chiffrement :
   * 1. y = m * G_pub
   * 2. e <- vecteur d'erreur de poids t
   * 3. c = y XOR e
   */

  /* 1. Calculer y = m * G_pub */
  memset(c, 0, (pk->G_pub.cols + 7) / 8);
  binvec_mul_mat(c, m, &pk->G_pub);

  /* 2. Generer le vecteur d'erreur e */
  uint8_t *e = (uint8_t *)calloc((pk->G_pub.cols + 7) / 8, sizeof(uint8_t));
  if (!e) {
    return -1;
  }

  if (generate_error_vector(e, pk->G_pub.cols, pk->t, u8rnd) != 0) {
    free(e);
    return -1;
  }

  /* 3. c = y XOR e */
  for (size_t i = 0; i < (pk->G_pub.cols + 7) / 8; i++) {
    c[i] ^= e[i];
  }

  free(e);
  return 0;
}
