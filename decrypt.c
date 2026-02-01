#include "decrypt.h"

int mc_decode(uint8_t *y, const uint8_t *c, const mc_secret_key_t *sk, const gf_ctx_t *ctx, int (*u8rnd)()) {
  (void)y;
  (void)c;
  (void)sk;
  (void)ctx;
  (void)u8rnd;
  return -1; /* TODO: implement Patterson (Alg. 4-7) */
}

int mc_decrypt(uint8_t *m, const uint8_t *c, const mc_secret_key_t *sk, const gf_ctx_t *ctx, int (*u8rnd)()) {
  (void)m;
  (void)c;
  (void)sk;
  (void)ctx;
  (void)u8rnd;
  return -1; /* TODO: implement Algorithm 3 */
}
