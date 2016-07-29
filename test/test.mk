LIBVPX_TEST_SRCS-yes += acm_random.h
LIBVPX_TEST_SRCS-yes += clear_system_state.h
LIBVPX_TEST_SRCS-yes += codec_factory.h
LIBVPX_TEST_SRCS-yes += md5_helper.h
LIBVPX_TEST_SRCS-yes += register_state_check.h
LIBVPX_TEST_SRCS-yes += test.mk
LIBVPX_TEST_SRCS-yes += test_libvpx.cc
LIBVPX_TEST_SRCS-yes += util.h
LIBVPX_TEST_SRCS-yes += video_source.h
LIBVPX_TEST_SRCS-yes += transform_test_base.h
LIBVPX_TEST_SRCS-yes += function_equivalence_test.h

##
## BLACK BOX TESTS
##
## Black box tests only use the public API.
##
LIBVPX_TEST_SRCS-yes                   += ../md5_utils.h ../md5_utils.c
LIBVPX_TEST_SRCS-$(CONFIG_DECODERS)    += ivf_video_source.h
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += ../y4minput.h ../y4minput.c
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += altref_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += aq_segment_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += datarate_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += encode_api_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += error_resilience_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += i420_video_source.h
#LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += realtime_test.cc
#LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += resize_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += y4m_video_source.h
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += yuv_video_source.h

#LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += level_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += active_map_refresh_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += active_map_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += borders_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += cpu_speed_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += frame_size_tests.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += lossless_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += ethread_test.cc

LIBVPX_TEST_SRCS-yes                   += decode_test_driver.cc
LIBVPX_TEST_SRCS-yes                   += decode_test_driver.h
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += encode_test_driver.cc
LIBVPX_TEST_SRCS-yes                   += encode_test_driver.h

## IVF writing.
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += ../ivfenc.c ../ivfenc.h

## Y4m parsing.
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS)    += y4m_test.cc ../y4menc.c ../y4menc.h

## WebM Parsing
ifeq ($(CONFIG_WEBM_IO), yes)
LIBWEBM_PARSER_SRCS += ../third_party/libwebm/mkvparser/mkvparser.cc
LIBWEBM_PARSER_SRCS += ../third_party/libwebm/mkvparser/mkvreader.cc
LIBWEBM_PARSER_SRCS += ../third_party/libwebm/mkvparser/mkvparser.h
LIBWEBM_PARSER_SRCS += ../third_party/libwebm/mkvparser/mkvreader.h
LIBVPX_TEST_SRCS-$(CONFIG_DECODERS)    += $(LIBWEBM_PARSER_SRCS)
LIBVPX_TEST_SRCS-$(CONFIG_DECODERS)    += ../tools_common.h
LIBVPX_TEST_SRCS-$(CONFIG_DECODERS)    += ../webmdec.cc
LIBVPX_TEST_SRCS-$(CONFIG_DECODERS)    += ../webmdec.h
LIBVPX_TEST_SRCS-$(CONFIG_DECODERS)    += webm_video_source.h
endif

LIBVPX_TEST_SRCS-$(CONFIG_DECODERS)    += decode_api_test.cc

# Currently we only support decoder perf tests for vp9. Also they read from WebM
# files, so WebM IO is required.
ifeq ($(CONFIG_DECODE_PERF_TESTS)$(CONFIG_VP10_DECODER)$(CONFIG_WEBM_IO), \
      yesyesyes)
LIBVPX_TEST_SRCS-yes                   += decode_perf_test.cc
endif

# encode perf tests are vp9 only
ifeq ($(CONFIG_ENCODE_PERF_TESTS)$(CONFIG_VP10_ENCODER), yesyes)
LIBVPX_TEST_SRCS-yes += encode_perf_test.cc
endif

## Multi-codec / unconditional black box tests.
ifeq ($(findstring yes,$(CONFIG_VP10_ENCODER)),yes)
LIBVPX_TEST_SRCS-yes += active_map_refresh_test.cc
LIBVPX_TEST_SRCS-yes += active_map_test.cc
LIBVPX_TEST_SRCS-yes += end_to_end_test.cc
endif

##
## WHITE BOX TESTS
##
## Whitebox tests invoke functions not exposed via the public API. Certain
## shared library builds don't make these functions accessible.
##
ifeq ($(CONFIG_SHARED),)

## VP10
ifeq ($(CONFIG_VP10),yes)

# These tests require both the encoder and decoder to be built.
ifeq ($(CONFIG_VP10_ENCODER)$(CONFIG_VP10_DECODER),yesyes)
# IDCT test currently depends on FDCT function
LIBVPX_TEST_SRCS-yes                   += idct8x8_test.cc
LIBVPX_TEST_SRCS-yes                   += partial_idct_test.cc
LIBVPX_TEST_SRCS-yes                   += superframe_test.cc
LIBVPX_TEST_SRCS-yes                   += tile_independence_test.cc
LIBVPX_TEST_SRCS-yes                   += boolcoder_test.cc
#LIBVPX_TEST_SRCS-yes                   += encoder_parms_get_to_decoder.cc
endif

#LIBVPX_TEST_SRCS-yes                   += convolve_test.cc
LIBVPX_TEST_SRCS-yes                   += lpf_8_test.cc
LIBVPX_TEST_SRCS-yes                   += intrapred_test.cc
#LIBVPX_TEST_SRCS-$(CONFIG_VP10_DECODER) += vp9_thread_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += dct16x16_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += dct32x32_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += fdct4x4_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += fdct8x8_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += hadamard_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += minmax_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += variance_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += error_block_test.cc
#LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += vp9_quantize_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += subtract_test.cc

ifeq ($(CONFIG_VP10_ENCODER),yes)
LIBVPX_TEST_SRCS-$(CONFIG_SPATIAL_SVC) += svc_test.cc
endif

ifeq ($(CONFIG_VP10_ENCODER)$(CONFIG_VP10_TEMPORAL_DENOISING),yesyes)
LIBVPX_TEST_SRCS-$(HAVE_SSE2) += denoiser_sse2_test.cc
endif
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += arf_freq_test.cc


LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += vp10_inv_txfm_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += vp10_dct_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += vp10_fht4x4_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += vp10_fht8x8_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += vp10_fht16x16_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_ANS)          += vp10_ans_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_EXT_TILE)     += vp10_ext_tile_test.cc

LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += sum_squares_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += subtract_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += blend_a64_mask_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += blend_a64_mask_1d_test.cc

ifeq ($(CONFIG_EXT_INTER),yes)
LIBVPX_TEST_SRCS-$(HAVE_SSSE3) += masked_variance_test.cc
LIBVPX_TEST_SRCS-$(HAVE_SSSE3) += masked_sad_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += vp10_wedge_utils_test.cc
endif

ifeq ($(CONFIG_OBMC),yes)
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += obmc_sad_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10_ENCODER) += obmc_variance_test.cc
endif

ifeq ($(CONFIG_VPX_HIGHBITDEPTH),yes)
LIBVPX_TEST_SRCS-$(HAVE_SSE4_1) += vp10_quantize_test.cc
LIBVPX_TEST_SRCS-$(HAVE_SSE4_1) += vp10_highbd_iht_test.cc
endif # CONFIG_VPX_HIGHBITDEPTH
endif # VP10

## Multi-codec / unconditional whitebox tests.

ifeq ($(CONFIG_VP10_ENCODER),yes)
LIBVPX_TEST_SRCS-yes += avg_test.cc
endif
ifeq ($(CONFIG_INTERNAL_STATS),yes)
LIBVPX_TEST_SRCS-$(CONFIG_VPX_HIGHBITDEPTH) += hbd_metrics_test.cc
endif
LIBVPX_TEST_SRCS-$(CONFIG_ENCODERS) += sad_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10) += vp10_txfm_test.h
LIBVPX_TEST_SRCS-$(CONFIG_VP10) += vp10_txfm_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10) += vp10_fwd_txfm1d_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10) += vp10_inv_txfm1d_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10) += vp10_fwd_txfm2d_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10) += vp10_inv_txfm2d_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10) += vp10_convolve_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP10) += vp10_convolve_optimz_test.cc

TEST_INTRA_PRED_SPEED_SRCS-yes := test_intra_pred_speed.cc
TEST_INTRA_PRED_SPEED_SRCS-yes += ../md5_utils.h ../md5_utils.c

endif # CONFIG_SHARED

include $(SRC_PATH_BARE)/test/test-data.mk
