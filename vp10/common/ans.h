#ifndef VP10_COMMON_ANS_H_
#define VP10_COMMON_ANS_H_
// An implementation of Asymmetric Numeral Systems
// http://arxiv.org/abs/1311.2540v2

#include <stdint.h>
#include "vpx_ports/mem_ops.h"

#define ANS_DIVIDE_BY_MULTIPLY 0
#if ANS_DIVIDE_BY_MULTIPLY
#include "divide.h"
#define ANS_INIT_DIVIDE init_fastdiv()
#define ANS_DIVREM(quotient, remainder, dividend, divisor) \
  do { \
    quotient = FASTDIV(dividend, divisor); \
    remainder = dividend - quotient * divisor; \
  } while (0)
#else
#define ANS_INIT_DIVIDE
#define ANS_DIVREM(quotient, remainder, dividend, divisor) \
  do { \
    quotient = dividend / divisor; \
    remainder = dividend % divisor; \
  } while (0)
#endif

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

struct AnsCoder {
  uint8_t *buf;
  int buf_offset;
  uint32_t state;
};

struct AnsDecoder {
  const uint8_t *buf;
  int buf_offset;
  uint32_t state;
};

typedef uint8_t AnsP8;
#define ans_p8_precision 256
#define ans_p8_shift 8
#define l_base (ans_p8_precision * 4)  // l_base % precision must be 0
#define io_base 256
// Range I = { l_base, l_base + 1, ..., l_base * io_base - 1 }

static inline void ans_write_init(struct AnsCoder *const ans,
                                  uint8_t *const buf) {
  ans->buf = buf;
  ans->buf_offset = 0;
  ans->state = l_base;
}

static inline int ans_write_end(struct AnsCoder *const ans) {
  mem_put_le24(ans->buf + ans->buf_offset, ans->state);
  return ans->buf_offset + 3;
}

// rABS with normalization
// p or p0 takes the place of l_s from the paper
// ans_p8_precision is m
static inline void rabs_write(struct AnsCoder *ans, int val, AnsP8 p0) {
  const AnsP8 p = ans_p8_precision - p0;
  const unsigned l_s = val ? p : p0;
  unsigned quot, rem;
  if ((!val && ans->state >= l_base / ans_p8_precision * io_base * p0) ||
      (val && ans->state >= l_base / ans_p8_precision * io_base * p)) {
    ans->buf[ans->buf_offset++] = ans->state % io_base;
    ans->state /= io_base;
  }
  ANS_DIVREM(quot, rem, ans->state, l_s);
  ans->state = quot * ans_p8_precision + rem + (val ? 0 : p);
}

#define ANS_IMPL1 0
#define UNPREDICTABLE(x) x
static inline int rabs_read(struct AnsDecoder *ans, AnsP8 p0) {
  int val;
#if ANS_IMPL1
  unsigned l_s;
#else
  unsigned quot, rem, x, xn;
#endif
  const AnsP8 p = ans_p8_precision - p0;
  if (ans->state < l_base) {
    ans->state = ans->state * io_base + ans->buf[--ans->buf_offset];
  }
#if ANS_IMPL1
  val = ans->state % ans_p8_precision < p;
  l_s = val ? p : p0;
  ans->state = (ans->state / ans_p8_precision) * l_s +
               ans->state % ans_p8_precision - (!val * p);
#else
  x = ans->state;
  quot = x / ans_p8_precision;
  rem = x % ans_p8_precision;
  xn = quot * p;
  val = rem < p;
  if (UNPREDICTABLE(val)) {
  ans->state = xn + rem;
  } else {
  //ans->state = quot * p0 + rem - p;
  ans->state = x - xn - p;
  }
#endif
  return val;
}

struct rans_sym {
  AnsP8 prob;
  AnsP8 cum_prob;  // not-inclusive
};

struct rans_dec_sym {
  uint8_t val;
  AnsP8 prob;
  AnsP8 cum_prob;  // not-inclusive
};

static inline void rans_build_dec_tab(const AnsP8 token_probs[], struct rans_dec_sym dec_tab[]) {
  int val = 0;
  int cum_prob = 0;
  int sym_end = token_probs[0];
  int i;
  for (i = 0; i < 256; ++i) {
    if (i == sym_end) {
      ++val;
      cum_prob = sym_end;
      sym_end += token_probs[val];
    }
    dec_tab[i].val = val;
    dec_tab[i].prob = token_probs[val];
    dec_tab[i].cum_prob = cum_prob;
  }
}

#define DBG_RANS 0
// rANS with normalization
// sym->prob takes the place of l_s from the paper
// ans_p8_precision is m
static inline void rans_stream_encode(struct AnsCoder *ans, const struct rans_sym *const sym) {
  const AnsP8 p = sym->prob;
//  const unsigned int s0 = ans->state;
  while (ans->state >= l_base / ans_p8_precision * io_base * p) {
    ans->buf[ans->buf_offset++] = ans->state % io_base;
    ans->state /= io_base;
  }
#if DBG_RANS
  unsigned state0 = ans->state;
#endif
  ans->state = (ans->state / p) * ans_p8_precision + ans->state % p + sym->cum_prob;
#if DBG_RANS
  fprintf(stderr, "C(val = [%02x %02x], %x) = %x\n", sym->cum_prob, sym->prob, state0, ans->state);
#endif
}

static inline int rans_stream_decode(struct AnsDecoder *ans,
                                     const struct rans_dec_sym tab[]) {
  unsigned rem;
  unsigned quo;
  int val;
//  const unsigned int s0 = ans->state;
  if (ans->state < l_base) {
    ans->state = ans->state * io_base + ans->buf[--ans->buf_offset];
  }
#if DBG_RANS
  unsigned state0 = ans->state;
#endif
  quo = ans->state / ans_p8_precision;
  rem = ans->state % ans_p8_precision;
  val = tab[rem].val;

  ans->state = quo * tab[rem].prob + rem - tab[rem].cum_prob;
#if DBG_RANS
  fprintf(stderr, "D(%x) = (%d = [%02x %02x], %x)\n", state0, val, tab[rem].cum_prob, tab[rem].prob, ans->state);
#endif
  return val;
}

static inline int ans_read_init(struct AnsDecoder *const ans,
                                 const uint8_t *const buf,
                                 int offset) {
  if (offset < 3)
    return 1;
  ans->buf = buf;
  ans->buf_offset = offset - 3;
  ans->state = mem_get_le24(buf + offset - 3);
  return 0;
}

static inline int ans_read_end(struct AnsDecoder *const ans) {
  return ans->state == l_base;
}
#undef ANS_DIVREM
#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
#endif  // VP10_COMMON_ANS_H_
