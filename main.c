#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "param.h"
#include "gf.h"
#include "key_gen.h"
#include "encrypt.h"
#include "decrypt.h"
#include "rng.h"

static int u8rnd(void) {
  unsigned char b = 0;
  if (randombytes(&b, 1) != RNG_SUCCESS) {
    return -1;
  }
  return (int)b;
}

static void rand_vec(uint8_t *v, size_t num_bits) {
  size_t bytes = (num_bits + 7) / 8;
  for (size_t i = 0; i < bytes; i++) {
    v[i] = (uint8_t)u8rnd();
  }
  if (num_bits % 8) {
    v[bytes - 1] &= (1u << (num_bits % 8)) - 1;
  }
}

int main(void) {
  printf("=== Démarrage du test McEliece ===\n");
  fflush(stdout);
  
  mc_params_t params = MC_PARAMS_TOY;
  gf_ctx_t gf = {0};
  mc_public_key_t pk = {0};
  mc_secret_key_t sk = {0};

  printf("Initialisation du corps de Galois...\n");
  fflush(stdout);
  
  if (gf_init(&gf, &params) != 0) {
    fprintf(stderr, "gf_init failed\n");
    fflush(stderr);
    return 1;
  }
  
  printf("Corps de Galois initialisé avec succès\n");
  fflush(stdout);

  printf("Génération des clés...\n");
  fflush(stdout);
  
  if (mc_keygen(&params, &pk, &sk, u8rnd) != 0) {
    fprintf(stderr, "mc_keygen failed\n");
    fflush(stderr);
    gf_clear(&gf);
    return 1;
  }
  
  printf("Clés générées avec succès\n");
  fflush(stdout);

  size_t k_bytes = ((size_t)params.k + 7) / 8;
  size_t n_bytes = ((size_t)params.n + 7) / 8;

  uint8_t *m  = calloc(k_bytes, 1);
  uint8_t *c  = calloc(n_bytes, 1);
  uint8_t *m2 = calloc(k_bytes, 1);
  if (!m || !c || !m2) {
    fprintf(stderr, "alloc failed\n");
    fflush(stderr);
    free(m);
    free(c);
    free(m2);
    gf_clear(&gf);
    return 1;
  }

  printf("Génération du message aléatoire...\n");
  fflush(stdout);
  
  rand_vec(m, (size_t)params.k);
  
  printf("Chiffrement du message...\n");
  fflush(stdout);

  if (mc_encrypt(c, m, &pk, u8rnd) != 0) {
    fprintf(stderr, "mc_encrypt failed\n");
    fflush(stderr);
    free(m);
    free(c);
    free(m2);
    gf_clear(&gf);
    return 1;
  }
  
  printf("Déchiffrement du cryptogramme...\n");
  fflush(stdout);

  if (mc_decrypt(m2, c, &sk, &gf, u8rnd) != 0) {
    fprintf(stderr, "mc_decrypt failed\n");
    fflush(stderr);
    free(m);
    free(c);
    free(m2);
    gf_clear(&gf);
    return 1;
  }
  
  printf("Vérification du résultat...\n");
  fflush(stdout);

  if (memcmp(m, m2, k_bytes) == 0) {
    printf("✓ OK: decrypt(encrypt(m)) == m\n");
  } else {
    fprintf(stderr, "decrypt mismatch\n");
    fflush(stderr);
  }
  
  free(m);
  free(c);
  free(m2);
  gf_clear(&gf);
  
  printf("=== Test réussi ===\n");
  fflush(stdout);
  
  return 0;
}
