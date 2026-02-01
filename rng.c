#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "rng.h"

#if defined(__APPLE__)
#include <stdlib.h>
#else
#include <fcntl.h>
#include <unistd.h>
#if defined(__linux__)
#include <sys/random.h>
#endif
#endif

int randombytes(unsigned char *x, unsigned long long xlen) {
  if (!x && xlen > 0) {
    return RNG_BAD_OUTBUF;
  }

#if defined(__APPLE__)
  arc4random_buf(x, (size_t)xlen);
  return RNG_SUCCESS;
#else
#if defined(__linux__)
  ssize_t r = getrandom(x, (size_t)xlen, 0);
  if (r == (ssize_t)xlen) {
    return RNG_SUCCESS;
  }
#endif
  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) {
    return RNG_BAD_OUTBUF;
  }
  size_t off = 0;
  while (off < xlen) {
    ssize_t n = read(fd, x + off, (size_t)(xlen - off));
    if (n < 0) {
      if (errno == EINTR) {
        continue;
      }
      close(fd);
      return RNG_BAD_OUTBUF;
    }
    off += (size_t)n;
  }
  close(fd);
  return RNG_SUCCESS;
#endif
}
