#include <stdlib.h>
#include <string.h>

#include "gf.h"

static gf_t gf_mul_slow(const gf_ctx_t *ctx, gf_t a, gf_t b) {
  gf_t res = 0;
  gf_t aa = a;
  gf_t bb = b;
  int m = ctx->m;
  gf_t mod = (gf_t)ctx->prim_poly;
  gf_t hi_bit = (gf_t)(1u << (m - 1));

  while (bb) {
    if (bb & 1u) {
      res ^= aa;
    }
    bb >>= 1;
    if (aa & hi_bit) {
      aa <<= 1;
      aa ^= mod;
    } else {
      aa <<= 1;
    }
  }
  return res;
}

int gf_init(gf_ctx_t *ctx, const mc_params_t *p) {
  if (!ctx || !p) {
    return -1;
  }
  memset(ctx, 0, sizeof(*ctx));
  ctx->m = p->m;
  ctx->q = 1 << p->m;
  ctx->q_minus_1 = ctx->q - 1;
  ctx->prim_poly = p->field_poly;

  ctx->exp = (gf_t *)calloc((size_t)ctx->q, sizeof(gf_t));
  ctx->log = (gf_t *)calloc((size_t)ctx->q, sizeof(gf_t));
  if (!ctx->exp || !ctx->log) {
    gf_clear(ctx);
    return -1;
  }

  ctx->exp[0] = 1;
  for (int i = 1; i < ctx->q_minus_1; ++i) {
    gf_t x = (gf_t)(ctx->exp[i - 1] << 1);
    if (ctx->exp[i - 1] & (1u << (ctx->m - 1))) {
      x ^= (gf_t)ctx->prim_poly;
    }
    ctx->exp[i] = (gf_t)(x & ctx->q_minus_1);
  }
  ctx->exp[ctx->q_minus_1] = 1;

  ctx->log[0] = (gf_t)ctx->q_minus_1;
  for (int i = 0; i < ctx->q_minus_1; ++i) {
    ctx->log[ctx->exp[i]] = (gf_t)i;
  }

  return 0;
}

void gf_clear(gf_ctx_t *ctx) {
  if (!ctx) {
    return;
  }
  free(ctx->exp);
  free(ctx->log);
  ctx->exp = NULL;
  ctx->log = NULL;
}

gf_t gf_mul(const gf_ctx_t *ctx, gf_t a, gf_t b) {
  if (!a || !b) {
    return 0;
  }
  if (ctx && ctx->exp && ctx->log) {
    int idx = ctx->log[a] + ctx->log[b];
    if (idx >= ctx->q_minus_1) {
      idx -= ctx->q_minus_1;
    }
    return ctx->exp[idx];
  }
  return gf_mul_slow(ctx, a, b);
}

gf_t gf_square(const gf_ctx_t *ctx, gf_t a) {
  return gf_mul(ctx, a, a);
}

gf_t gf_inv(const gf_ctx_t *ctx, gf_t a) {
  if (!a || !ctx || !ctx->exp || !ctx->log) {
    return 0;
  }
  int idx = ctx->q_minus_1 - ctx->log[a];
  if (idx == ctx->q_minus_1) {
    idx = 0;
  }
  return ctx->exp[idx];
}

gf_t gf_pow(const gf_ctx_t *ctx, gf_t a, int e) {
  if (e == 0) {
    return 1;
  }
  if (a == 0) {
    return 0;
  }
  if (ctx && ctx->exp && ctx->log) {
    int exp = ctx->log[a];
    long long ee = (long long)exp * (long long)e;
    long long mod = ctx->q_minus_1;
    long long r = ee % mod;
    if (r < 0) {
      r += mod;
    }
    return ctx->exp[(int)r];
  }

  gf_t base = a;
  gf_t res = 1;
  int exp = e;
  while (exp > 0) {
    if (exp & 1) {
      res = gf_mul(ctx, res, base);
    }
    base = gf_mul(ctx, base, base);
    exp >>= 1;
  }
  return res;
}

gf_t gf_sqrt(const gf_ctx_t *ctx, gf_t a) {
  if (!ctx) {
    return 0;
  }
  int e = 1 << (ctx->m - 1);
  return gf_pow(ctx, a, e);
}
