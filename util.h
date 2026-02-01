#ifndef UTIL_H
#define UTIL_H

#include <stddef.h>
#include <stdint.h>

/*
 * Poids de Hamming d'un vecteur binaire de taille n.
 */
size_t hamming_weight(const uint8_t *v, size_t n);

/*
 * XOR de deux vecteurs binaires : out = a XOR b (taille n).
 */
void xor_vec(uint8_t *out, const uint8_t *a, const uint8_t *b, size_t n);

#endif /* UTIL_H */
