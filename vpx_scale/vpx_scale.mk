SCALE_SRCS-yes += vpx_scale.mk
SCALE_SRCS-yes += scale_mode.h
SCALE_SRCS-yes += yv12config.h
SCALE_SRCS-yes += vpxscale.h
SCALE_SRCS-yes += generic/vpxscale.c
SCALE_SRCS-yes += generic/yv12config.c
SCALE_SRCS-yes += generic/yv12extend.c
SCALE_SRCS-$(CONFIG_SPATIAL_RESAMPLING) += generic/gen_scalers.c
SCALE_SRCS-yes += rtcd.c

#neon
SCALE_SRCS-$(HAVE_NEON)  += arm/neon/vp8_vpxyv12_copyframe_func_neon$(ASM)
SCALE_SRCS-$(HAVE_NEON)  += arm/neon/vp8_vpxyv12_copy_y_neon$(ASM)
SCALE_SRCS-$(HAVE_NEON)  += arm/neon/vp8_vpxyv12_copysrcframe_func_neon$(ASM)
SCALE_SRCS-$(HAVE_NEON)  += arm/neon/vp8_vpxyv12_extendframeborders_neon$(ASM)
SCALE_SRCS-$(HAVE_NEON)  += arm/neon/yv12extend_arm.c

SCALE_SRCS-no += $(SCALE_SRCS_REMOVE-yes)

$(eval $(call rtcd_h_template,vpx_scale_rtcd,vpx_scale/vpx_scale_rtcd.sh))
