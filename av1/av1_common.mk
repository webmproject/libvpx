##
## Copyright (c) 2016, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##

AV1_COMMON_SRCS-yes += av1_common.mk
AV1_COMMON_SRCS-yes += av1_iface_common.h
AV1_COMMON_SRCS-yes += common/alloccommon.c
AV1_COMMON_SRCS-yes += common/blockd.c
AV1_COMMON_SRCS-yes += common/debugmodes.c
AV1_COMMON_SRCS-yes += common/entropy.c
AV1_COMMON_SRCS-yes += common/entropymode.c
AV1_COMMON_SRCS-yes += common/entropymv.c
AV1_COMMON_SRCS-yes += common/frame_buffers.c
AV1_COMMON_SRCS-yes += common/frame_buffers.h
AV1_COMMON_SRCS-yes += common/alloccommon.h
AV1_COMMON_SRCS-yes += common/blockd.h
AV1_COMMON_SRCS-yes += common/common.h
AV1_COMMON_SRCS-yes += common/entropy.h
AV1_COMMON_SRCS-yes += common/entropymode.h
AV1_COMMON_SRCS-yes += common/entropymv.h
AV1_COMMON_SRCS-yes += common/enums.h
AV1_COMMON_SRCS-yes += common/filter.h
AV1_COMMON_SRCS-yes += common/filter.c
AV1_COMMON_SRCS-yes += common/idct.h
AV1_COMMON_SRCS-yes += common/idct.c
AV1_COMMON_SRCS-yes += common/loopfilter.h
AV1_COMMON_SRCS-yes += common/thread_common.h
AV1_COMMON_SRCS-yes += common/mv.h
AV1_COMMON_SRCS-yes += common/onyxc_int.h
AV1_COMMON_SRCS-yes += common/pred_common.h
AV1_COMMON_SRCS-yes += common/pred_common.c
AV1_COMMON_SRCS-yes += common/quant_common.h
AV1_COMMON_SRCS-yes += common/reconinter.h
AV1_COMMON_SRCS-yes += common/reconintra.h
AV1_COMMON_SRCS-yes += common/av1_rtcd.c
AV1_COMMON_SRCS-yes += common/av1_rtcd_defs.pl
AV1_COMMON_SRCS-yes += common/scale.h
AV1_COMMON_SRCS-yes += common/scale.c
AV1_COMMON_SRCS-yes += common/seg_common.h
AV1_COMMON_SRCS-yes += common/seg_common.c
AV1_COMMON_SRCS-yes += common/tile_common.h
AV1_COMMON_SRCS-yes += common/tile_common.c
AV1_COMMON_SRCS-yes += common/loopfilter.c
AV1_COMMON_SRCS-yes += common/thread_common.c
AV1_COMMON_SRCS-yes += common/mvref_common.c
AV1_COMMON_SRCS-yes += common/mvref_common.h
AV1_COMMON_SRCS-yes += common/quant_common.c
AV1_COMMON_SRCS-yes += common/reconinter.c
AV1_COMMON_SRCS-yes += common/reconintra.c
AV1_COMMON_SRCS-yes += common/restoration.h
AV1_COMMON_SRCS-yes += common/common_data.h
AV1_COMMON_SRCS-yes += common/scan.c
AV1_COMMON_SRCS-yes += common/scan.h
# TODO(angiebird) the forward transform belongs under encoder/
AV1_COMMON_SRCS-yes += common/av1_txfm.h
AV1_COMMON_SRCS-yes += common/av1_fwd_txfm1d.h
AV1_COMMON_SRCS-yes += common/av1_fwd_txfm1d.c
AV1_COMMON_SRCS-yes += common/av1_inv_txfm1d.h
AV1_COMMON_SRCS-yes += common/av1_inv_txfm1d.c
AV1_COMMON_SRCS-yes += common/av1_fwd_txfm2d.c
AV1_COMMON_SRCS-yes += common/av1_fwd_txfm2d_cfg.h
AV1_COMMON_SRCS-yes += common/av1_inv_txfm2d.c
AV1_COMMON_SRCS-yes += common/av1_inv_txfm2d_cfg.h
AV1_COMMON_SRCS-$(HAVE_SSSE3) += common/x86/av1_convolve_ssse3.c
AV1_COMMON_SRCS-$(HAVE_SSSE3) += common/x86/av1_convolve_filters_ssse3.h
ifeq ($(CONFIG_AOM_HIGHBITDEPTH),yes)
AV1_COMMON_SRCS-$(HAVE_SSE4_1) += common/x86/av1_highbd_convolve_sse4.c
AV1_COMMON_SRCS-$(HAVE_SSE4_1) += common/x86/av1_highbd_convolve_filters_sse4.h
endif
AV1_COMMON_SRCS-yes += common/convolve.c
AV1_COMMON_SRCS-yes += common/convolve.h
AV1_COMMON_SRCS-$(CONFIG_LOOP_RESTORATION) += common/restoration.h
AV1_COMMON_SRCS-$(CONFIG_LOOP_RESTORATION) += common/restoration.c
ifeq (yes,$(filter $(CONFIG_GLOBAL_MOTION) $(CONFIG_WARPED_MOTION),yes))
AV1_COMMON_SRCS-yes += common/warped_motion.h
AV1_COMMON_SRCS-yes += common/warped_motion.c
endif
ifeq ($(CONFIG_CLPF),yes)
AV1_COMMON_SRCS-yes += common/clpf.c
AV1_COMMON_SRCS-yes += common/clpf.h
AV1_COMMON_SRCS-yes += common/clpf_simd.h
AV1_COMMON_SRCS-$(HAVE_SSE2) += common/clpf_sse2.c
AV1_COMMON_SRCS-$(HAVE_SSSE3) += common/clpf_ssse3.c
AV1_COMMON_SRCS-$(HAVE_SSE4_1) += common/clpf_sse4.c
AV1_COMMON_SRCS-$(HAVE_NEON) += common/clpf_neon.c
endif
ifeq ($(CONFIG_DERING),yes)
AV1_COMMON_SRCS-yes += common/od_dering.c
AV1_COMMON_SRCS-yes += common/od_dering.h
AV1_COMMON_SRCS-$(HAVE_SSE4_1) += common/x86/od_dering_sse4.c
AV1_COMMON_SRCS-$(HAVE_SSE4_1) += common/x86/od_dering_sse4.h
AV1_COMMON_SRCS-yes += common/dering.c
AV1_COMMON_SRCS-yes += common/dering.h
endif
ifeq ($(CONFIG_ACCOUNTING),yes)
AV1_COMMON_SRCS-yes += common/accounting.h
AV1_COMMON_SRCS-yes += common/accounting.c
endif
AV1_COMMON_SRCS-yes += common/odintrin.c
AV1_COMMON_SRCS-yes += common/odintrin.h

ifeq ($(CONFIG_PVQ),yes)
# PVQ from daala
AV1_COMMON_SRCS-yes += common/pvq.c
AV1_COMMON_SRCS-yes += common/pvq.h
AV1_COMMON_SRCS-yes += common/partition.c
AV1_COMMON_SRCS-yes += common/partition.h
AV1_COMMON_SRCS-yes += common/zigzag4.c
AV1_COMMON_SRCS-yes += common/zigzag8.c
AV1_COMMON_SRCS-yes += common/zigzag16.c
AV1_COMMON_SRCS-yes += common/zigzag32.c
AV1_COMMON_SRCS-yes += common/zigzag.h
AV1_COMMON_SRCS-yes += common/generic_code.c
AV1_COMMON_SRCS-yes += common/generic_code.h
AV1_COMMON_SRCS-yes += common/pvq_state.c
AV1_COMMON_SRCS-yes += common/pvq_state.h
AV1_COMMON_SRCS-yes += common/laplace_tables.c
endif

ifneq ($(CONFIG_AOM_HIGHBITDEPTH),yes)
AV1_COMMON_SRCS-$(HAVE_DSPR2)  += common/mips/dspr2/itrans4_dspr2.c
AV1_COMMON_SRCS-$(HAVE_DSPR2)  += common/mips/dspr2/itrans8_dspr2.c
AV1_COMMON_SRCS-$(HAVE_DSPR2)  += common/mips/dspr2/itrans16_dspr2.c
endif

# common (msa)
AV1_COMMON_SRCS-$(HAVE_MSA) += common/mips/msa/idct4x4_msa.c
AV1_COMMON_SRCS-$(HAVE_MSA) += common/mips/msa/idct8x8_msa.c
AV1_COMMON_SRCS-$(HAVE_MSA) += common/mips/msa/idct16x16_msa.c

AV1_COMMON_SRCS-$(HAVE_SSE2) += common/x86/idct_intrin_sse2.c
AV1_COMMON_SRCS-$(HAVE_AVX2) += common/x86/hybrid_inv_txfm_avx2.c

ifeq ($(CONFIG_AV1_ENCODER),yes)
AV1_COMMON_SRCS-$(HAVE_SSE4_1) += common/x86/av1_txfm1d_sse4.h
AV1_COMMON_SRCS-$(HAVE_SSE4_1) += common/x86/av1_fwd_txfm1d_sse4.c
AV1_COMMON_SRCS-$(HAVE_SSE4_1) += common/x86/av1_fwd_txfm2d_sse4.c
endif
ifeq ($(CONFIG_AOM_HIGHBITDEPTH),yes)
AV1_COMMON_SRCS-$(HAVE_SSE4_1) += common/x86/highbd_txfm_utility_sse4.h
endif

ifneq ($(CONFIG_AOM_HIGHBITDEPTH),yes)
AV1_COMMON_SRCS-$(HAVE_NEON) += common/arm/neon/iht4x4_add_neon.c
AV1_COMMON_SRCS-$(HAVE_NEON) += common/arm/neon/iht8x8_add_neon.c
endif

ifeq ($(CONFIG_FILTER_INTRA),yes)
AV1_COMMON_SRCS-$(HAVE_SSE4_1) += common/x86/filterintra_sse4.c
endif

$(eval $(call rtcd_h_template,av1_rtcd,av1/common/av1_rtcd_defs.pl))
