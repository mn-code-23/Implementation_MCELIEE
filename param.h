#ifndef PARAM_H
#define PARAM_H

#include <stdint.h>

/*
 * Parametres McEliece / Goppa.
 * - m : degre du corps F2^m (n = 2^m)
 * - t : capacite de correction (deg(g) = t)
 * - n : longueur du code
 * - k : dimension du code (k = n - m*t)
 * - field_poly : polynome irreductible de F2^m (bitmask)
 */
typedef struct {
  int m;
  int t;
  int n;
  int k;
  uint32_t field_poly;
} mc_params_t;

/*
 * Initialise les champs derives (n, k) a partir de m et t.
 * Operation mathematique : n = 2^m, k = n - m*t.
 */
void mc_params_init(mc_params_t *p, int m, int t, uint32_t field_poly);

/*
 * Exemple de jeu de parametres (petit, pour tests rapides).
 * A adapter selon le cours.
 */
extern const mc_params_t MC_PARAMS_TOY;

#endif /* PARAM_H */
