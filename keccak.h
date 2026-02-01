#define FOR(i,n) for(i=0; i<n; ++i)
typedef unsigned char u8;
typedef unsigned long long int u64;
typedef unsigned int ui;
#define ROL(a,o) ((((u64)a)<<o)^(((u64)a)>>(64-o)))
#define rL(x,y) load64((u8*)s+8*(x+5*y))
#define wL(x,y,l) store64((u8*)s+8*(x+5*y),l)
#define XL(x,y,l) xor64((u8*)s+8*(x+5*y),l)
void Keccak(ui r, ui c, const u8 *in, u64 inLen, u8 sfx, u8 *out, u64 outLen);
void FIPS202_SHAKE128(const u8 *in, u64 inLen, u8 *out, u64 outLen);
void FIPS202_SHAKE256(const u8 *in, u64 inLen, u8 *out, u64 outLen);
void FIPS202_SHA3_224(const u8 *in, u64 inLen, u8 *out);
void FIPS202_SHA3_256(const u8 *in, u64 inLen, u8 *out); 
void FIPS202_SHA3_384(const u8 *in, u64 inLen, u8 *out);
void FIPS202_SHA3_512(const u8 *in, u64 inLen, u8 *out);
int LFSR86540(u8 *R);
void KeccakF1600(void *s);