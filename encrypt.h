#ifndef ENCRYPT_H
#define ENCRYPT_H

#include <stddef.h>
#include <stdint.h>
#include "key_gen.h"

/*
 * Genere un vecteur d'erreur e de poids t.
 * Operation mathematique : e in F2^n, wt(e)=t.
 */
int generate_error_vector(uint8_t *e, size_t n, int t, int (*u8rnd)());

/*
 * Chiffrement (Alg. 2 du guide) :
 *   y = m * G_pub
 *   c = y XOR e
 */
int mc_encrypt(uint8_t *c, const uint8_t *m, const mc_public_key_t *pk, int (*u8rnd)());

#endif /* ENCRYPT_H */
