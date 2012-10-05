LIBVPX_TEST_SRCS-yes += acm_random.h
LIBVPX_TEST_SRCS-yes += test.mk
LIBVPX_TEST_SRCS-yes += test_libvpx.cc
LIBVPX_TEST_SRCS-yes += util.h

##
## BLACK BOX TESTS
##
## Black box tests only use the public API.
##
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += altref_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += config_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += cq_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += encode_test_driver.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += encode_test_driver.h
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += error_resilience_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += i420_video_source.h
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += keyframe_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += resize_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += video_source.h
LIBVPX_TEST_SRCS-$(CONFIG_VP8_DECODER) += decode_test_driver.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP8_DECODER) += decode_test_driver.h

##
## WHITE BOX TESTS
##
## Whitebox tests invoke functions not exposed via the public API. Certain
## shared library builds don't make these functions accessible.
##
ifeq ($(CONFIG_SHARED),)

# These tests require both the encoder and decoder to be built.
ifeq ($(CONFIG_VP8_ENCODER)$(CONFIG_VP8_DECODER),yesyes)
LIBVPX_TEST_SRCS-yes                   += boolcoder_test.cc
endif

LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += fdct4x4_test.cc
LIBVPX_TEST_SRCS-yes                   += idctllm_test.cc
LIBVPX_TEST_SRCS-yes                   += intrapred_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_POSTPROC)    += pp_filter_test.cc
LIBVPX_TEST_SRCS-yes                   += sad_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += set_roi.cc
LIBVPX_TEST_SRCS-yes                   += sixtap_predict_test.cc
LIBVPX_TEST_SRCS-$(CONFIG_VP8_ENCODER) += subtract_test.cc

endif


##
## TEST DATA
##
LIBVPX_TEST_DATA-$(CONFIG_VP8_ENCODER) += hantro_collage_w352h288.yuv
