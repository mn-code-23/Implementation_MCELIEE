#include "encrypt.h"

int mc_encrypt(uint8_t *c, const uint8_t *m, const mc_public_key_t *pk, int (*u8rnd)()) {
  (void)c;
  (void)m;
  (void)pk;
  (void)u8rnd;
  return -1; /* TODO: implement Algorithm 2 */
}
