#include "key_gen.h"

int mc_keygen(const mc_params_t *params, mc_public_key_t *pk, mc_secret_key_t *sk, int (*u8rnd)()) {
  (void)params;
  (void)pk;
  (void)sk;
  (void)u8rnd;
  return -1; /* TODO: implement Algorithm 1 */
}
