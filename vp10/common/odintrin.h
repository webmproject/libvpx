#ifndef VP10_COMMON_ODINTRIN_H_
#define VP10_COMMON_ODINTRIN_H_

#include "vp10/common/enums.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_ports/bitops.h"

/*Smallest blocks are 4x4*/
# define OD_LOG_BSIZE0 (2)
/*There are 5 block sizes total (4x4, 8x8, 16x16, 32x32 and 64x64).*/
# define OD_NBSIZES    (5)
/*The log of the maximum length of the side of a block.*/
# define OD_LOG_BSIZE_MAX (OD_LOG_BSIZE0 + OD_NBSIZES - 1)
/*The maximum length of the side of a block.*/
# define OD_BSIZE_MAX     (1 << OD_LOG_BSIZE_MAX)

typedef int od_coeff;

typedef int16_t od_dering_in;

# define OD_DIVU_DMAX (1024)

extern uint32_t OD_DIVU_SMALL_CONSTS[OD_DIVU_DMAX][2];

/*Translate unsigned division by small divisors into multiplications.*/
# define OD_DIVU_SMALL(_x, _d) \
  ((uint32_t)((OD_DIVU_SMALL_CONSTS[(_d)-1][0]* \
  (uint64_t)(_x)+OD_DIVU_SMALL_CONSTS[(_d)-1][1])>>32)>> \
  (OD_ILOG(_d)-1))

# define OD_DIVU(_x, _d) \
  (((_d) < OD_DIVU_DMAX)?(OD_DIVU_SMALL((_x), (_d))):((_x)/(_d)))

#define OD_MINI VPXMIN
#define OD_CLAMPI(min, val, max) clamp((val), (min), (max))

# define OD_CLZ0 (1)
# define OD_CLZ(x) (-get_msb(x))
# define OD_ILOG_NZ(x) (OD_CLZ0 - OD_CLZ(x))
/*Note that __builtin_clz is not defined when x == 0, according to the gcc
   documentation (and that of the x86 BSR instruction that implements it), so
   we have to special-case it.
  We define a special version of the macro to use when x can be zero.*/
# define OD_ILOG(x) ((x) ? OD_ILOG_NZ(x) : 0)

#endif
