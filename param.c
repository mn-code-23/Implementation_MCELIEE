#include "param.h"

void mc_params_init(mc_params_t *p, int m, int t, uint32_t field_poly) {
  if (!p) {
    return;
  }
  p->m = m;
  p->t = t;
  p->n = 1 << m;
  p->k = p->n - m * t;
  p->field_poly = field_poly;
}

const mc_params_t MC_PARAMS_TOY = {
  .m = 3,
  .t = 2,
  .n = 8,
  .k = 2,
  .field_poly = 0x43 /* x^6 + x + 1 */
};

/*
  .m = 6,
  .t = 5,
  .n = 64,
  .k = 34,
*/