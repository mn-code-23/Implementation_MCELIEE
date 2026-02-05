#include "util.h"

/* Poids de Hamming : compte le nombre de bits a 1 */
size_t hamming_weight(const uint8_t *v, size_t n) {
  if (!v) {
    return 0;
  }
  size_t weight = 0;
  for (size_t i = 0; i < n; i++) {
    uint8_t byte = v[i];
    while (byte) {
      weight += (byte & 1);
      byte >>= 1;
    }
  }
  return weight;
}

/* XOR de deux vecteurs binaires */
void xor_vec(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n) {
  if (!out || !a || !b) {
    return;
  }
  for (size_t i = 0; i < n; i++) {
    out[i] = a[i] ^ b[i];
  }
}
