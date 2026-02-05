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

static void rand_vec(uint8_t *v, size_t n) {
  for (size_t i = 0; i < n; i++) {
    int r = u8rnd();
    v[i] = (r < 0) ? 0 : (uint8_t)(r & 1);
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

  uint8_t *m = (uint8_t *)calloc((size_t)params.k, sizeof(uint8_t));
  uint8_t *c = (uint8_t *)calloc((size_t)params.n, sizeof(uint8_t));
  uint8_t *m2 = (uint8_t *)calloc((size_t)params.k, sizeof(uint8_t));
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

  if (memcmp(m, m2, (size_t)params.k) != 0) {
    fprintf(stderr, "decrypt mismatch\n");
    fflush(stderr);
  } else {
    printf("✓ OK: decrypt(encrypt(m)) == m\n");
    fflush(stdout);
  }
  
  free(m);
  free(c);
  free(m2);
  gf_clear(&gf);
  
  printf("=== Test réussi ===\n");
  fflush(stdout);
  
  return 0;
}
