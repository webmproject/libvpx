##
##  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##


# ARM assembly files are written in RVCT-style. We use some make magic to
# filter those files to allow GCC compilation
ifeq ($(ARCH_ARM),yes)
  ASM:=$(if $(filter yes,$(CONFIG_GCC)),.asm.s,.asm)
else
  ASM:=.asm
endif

CODEC_SRCS-yes += libs.mk

include $(SRC_PATH_BARE)/vpx/vpx_codec.mk
CODEC_SRCS-yes += $(addprefix vpx/,$(call enabled,API_SRCS))

include $(SRC_PATH_BARE)/vpx_mem/vpx_mem.mk
CODEC_SRCS-yes += $(addprefix vpx_mem/,$(call enabled,MEM_SRCS))

include $(SRC_PATH_BARE)/vpx_scale/vpx_scale.mk
CODEC_SRCS-yes += $(addprefix vpx_scale/,$(call enabled,SCALE_SRCS))


ifeq ($(CONFIG_VP8_ENCODER),yes)
  VP8_PREFIX=vp8/
  include $(SRC_PATH_BARE)/$(VP8_PREFIX)vp8cx.mk
  CODEC_SRCS-yes += $(addprefix $(VP8_PREFIX),$(call enabled,VP8_CX_SRCS))
  CODEC_EXPORTS-yes += $(addprefix $(VP8_PREFIX),$(VP8_CX_EXPORTS))
  CODEC_SRCS-yes += $(VP8_PREFIX)vp8cx.mk vpx/vp8.h vpx/vp8cx.h vpx/vp8e.h
  INSTALL-LIBS-yes += include/vpx/vp8.h include/vpx/vp8e.h include/vpx/vp8cx.h
  INSTALL_MAPS += include/vpx/% $(SRC_PATH_BARE)/$(VP8_PREFIX)/%
  CODEC_DOC_SRCS += vpx/vp8.h vpx/vp8cx.h
  CODEC_DOC_SECTIONS += vp8 vp8_encoder
endif

ifeq ($(CONFIG_VP8_DECODER),yes)
  VP8_PREFIX=vp8/
  include $(SRC_PATH_BARE)/$(VP8_PREFIX)vp8dx.mk
  CODEC_SRCS-yes += $(addprefix $(VP8_PREFIX),$(call enabled,VP8_DX_SRCS))
  CODEC_EXPORTS-yes += $(addprefix $(VP8_PREFIX),$(VP8_DX_EXPORTS))
  CODEC_SRCS-yes += $(VP8_PREFIX)vp8dx.mk vpx/vp8.h vpx/vp8dx.h
  INSTALL-LIBS-yes += include/vpx/vp8.h include/vpx/vp8dx.h
  INSTALL_MAPS += include/vpx/% $(SRC_PATH_BARE)/$(VP8_PREFIX)/%
  CODEC_DOC_SRCS += vpx/vp8.h vpx/vp8dx.h
  CODEC_DOC_SECTIONS += vp8 vp8_decoder
endif


ifeq ($(CONFIG_ENCODERS),yes)
  CODEC_DOC_SECTIONS += encoder
endif
ifeq ($(CONFIG_DECODERS),yes)
  CODEC_DOC_SECTIONS += decoder
endif


ifeq ($(CONFIG_MSVS),yes)
CODEC_LIB=$(if $(CONFIG_STATIC_MSVCRT),vpxmt,vpxmd)
# This variable uses deferred expansion intentionally, since the results of
# $(wildcard) may change during the course of the Make.
VS_PLATFORMS = $(foreach d,$(wildcard */Release/$(CODEC_LIB).lib),$(word 1,$(subst /, ,$(d))))
endif

# The following pairs define a mapping of locations in the distribution
# tree to locations in the source/build trees.
INSTALL_MAPS += include/vpx/% $(SRC_PATH_BARE)/vpx/%
INSTALL_MAPS += include/vpx/% $(SRC_PATH_BARE)/vpx_ports/%
INSTALL_MAPS += $(LIBSUBDIR)/%     %
INSTALL_MAPS += src/%     $(SRC_PATH_BARE)/%
ifeq ($(CONFIG_MSVS),yes)
INSTALL_MAPS += $(foreach p,$(VS_PLATFORMS),$(LIBSUBDIR)/$(p)/%  $(p)/Release/%)
INSTALL_MAPS += $(foreach p,$(VS_PLATFORMS),$(LIBSUBDIR)/$(p)/%  $(p)/Debug/%)
endif

# If this is a universal (fat) binary, then all the subarchitectures have
# already been built and our job is to stitch them together. The
# BUILD_LIBVPX variable indicates whether we should be building
# (compiling, linking) the library. The LIPO_LIBVPX variable indicates
# that we're stitching.
$(eval $(if $(filter universal%,$(TOOLCHAIN)),LIPO_LIBVPX,BUILD_LIBVPX):=yes)

CODEC_SRCS-$(BUILD_LIBVPX) += build/make/version.sh
CODEC_SRCS-$(BUILD_LIBVPX) += vpx/vpx_integer.h
CODEC_SRCS-$(BUILD_LIBVPX) += vpx_ports/vpx_timer.h
CODEC_SRCS-$(BUILD_LIBVPX) += vpx_ports/mem.h
CODEC_SRCS-$(BUILD_LIBVPX) += $(BUILD_PFX)vpx_config.c
INSTALL-SRCS-no += $(BUILD_PFX)vpx_config.c
ifeq ($(ARCH_X86)$(ARCH_X86_64),yes)
CODEC_SRCS-$(BUILD_LIBVPX) += vpx_ports/emms.asm
CODEC_SRCS-$(BUILD_LIBVPX) += vpx_ports/x86.h
CODEC_SRCS-$(BUILD_LIBVPX) += vpx_ports/x86_abi_support.asm
CODEC_SRCS-$(BUILD_LIBVPX) += vpx_ports/x86_cpuid.c
endif
CODEC_SRCS-$(ARCH_ARM) += vpx_ports/arm_cpudetect.c
CODEC_SRCS-$(ARCH_ARM) += $(BUILD_PFX)vpx_config.asm
CODEC_EXPORTS-$(BUILD_LIBVPX) += vpx/exports_com
CODEC_EXPORTS-$(CONFIG_ENCODERS) += vpx/exports_enc
CODEC_EXPORTS-$(CONFIG_DECODERS) += vpx/exports_dec

INSTALL-LIBS-yes += include/vpx/vpx_codec.h
INSTALL-LIBS-yes += include/vpx/vpx_image.h
INSTALL-LIBS-yes += include/vpx/vpx_integer.h
INSTALL-LIBS-yes += include/vpx/vpx_codec_impl_top.h
INSTALL-LIBS-yes += include/vpx/vpx_codec_impl_bottom.h
INSTALL-LIBS-$(CONFIG_DECODERS) += include/vpx/vpx_decoder.h
INSTALL-LIBS-$(CONFIG_DECODERS) += include/vpx/vpx_decoder_compat.h
INSTALL-LIBS-$(CONFIG_ENCODERS) += include/vpx/vpx_encoder.h
ifeq ($(CONFIG_EXTERNAL_BUILD),yes)
ifeq ($(CONFIG_MSVS),yes)
INSTALL-LIBS-yes                  += $(foreach p,$(VS_PLATFORMS),$(LIBSUBDIR)/$(p)/$(CODEC_LIB).lib)
INSTALL-LIBS-$(CONFIG_DEBUG_LIBS) += $(foreach p,$(VS_PLATFORMS),$(LIBSUBDIR)/$(p)/$(CODEC_LIB)d.lib)
INSTALL-LIBS-$(CONFIG_SHARED) += $(foreach p,$(VS_PLATFORMS),$(LIBSUBDIR)/$(p)/vpx.dll)
INSTALL-LIBS-$(CONFIG_SHARED) += $(foreach p,$(VS_PLATFORMS),$(LIBSUBDIR)/$(p)/vpx.exp)
endif
else
INSTALL-LIBS-yes += $(LIBSUBDIR)/libvpx.a
INSTALL-LIBS-$(CONFIG_DEBUG_LIBS) += $(LIBSUBDIR)/libvpx_g.a
endif

CODEC_SRCS=$(call enabled,CODEC_SRCS)
INSTALL-SRCS-$(CONFIG_CODEC_SRCS) += $(CODEC_SRCS)
INSTALL-SRCS-$(CONFIG_CODEC_SRCS) += $(call enabled,CODEC_EXPORTS)

ifeq ($(CONFIG_EXTERNAL_BUILD),yes)
ifeq ($(CONFIG_MSVS),yes)

obj_int_extract.vcproj: $(SRC_PATH_BARE)/build/make/obj_int_extract.c
	@cp $(SRC_PATH_BARE)/build/x86-msvs/obj_int_extract.bat .
	@echo "    [CREATE] $@"
	$(SRC_PATH_BARE)/build/make/gen_msvs_proj.sh \
    --exe \
    --target=$(TOOLCHAIN) \
    --name=obj_int_extract \
    --ver=$(CONFIG_VS_VERSION) \
    --proj-guid=E1360C65-D375-4335-8057-7ED99CC3F9B2 \
    $(if $(CONFIG_STATIC_MSVCRT),--static-crt) \
    --out=$@ $^ \
    -I. \
    -I"$(SRC_PATH_BARE)" \

PROJECTS-$(BUILD_LIBVPX) += obj_int_extract.vcproj
PROJECTS-$(BUILD_LIBVPX) += obj_int_extract.bat

vpx.def: $(call enabled,CODEC_EXPORTS)
	@echo "    [CREATE] $@"
	$(SRC_PATH_BARE)/build/make/gen_msvs_def.sh\
            --name=vpx\
            --out=$@ $^
CLEAN-OBJS += vpx.def

vpx.vcproj: $(CODEC_SRCS) vpx.def
	@echo "    [CREATE] $@"
	$(SRC_PATH_BARE)/build/make/gen_msvs_proj.sh \
			--lib \
			--target=$(TOOLCHAIN) \
            $(if $(CONFIG_STATIC_MSVCRT),--static-crt) \
            --name=vpx \
            --proj-guid=DCE19DAF-69AC-46DB-B14A-39F0FAA5DB74 \
            --module-def=vpx.def \
            --ver=$(CONFIG_VS_VERSION) \
            --out=$@ $(CFLAGS) $^ \
            --src-path-bare="$(SRC_PATH_BARE)" \

PROJECTS-$(BUILD_LIBVPX) += vpx.vcproj

vpx.vcproj: vpx_config.asm

endif
else
LIBVPX_OBJS=$(call objs,$(CODEC_SRCS))
OBJS-$(BUILD_LIBVPX) += $(LIBVPX_OBJS)
LIBS-$(BUILD_LIBVPX) += $(BUILD_PFX)libvpx.a $(BUILD_PFX)libvpx_g.a
$(BUILD_PFX)libvpx_g.a: $(LIBVPX_OBJS)

BUILD_LIBVPX_SO         := $(if $(BUILD_LIBVPX),$(CONFIG_SHARED))
LIBVPX_SO               := libvpx.so.$(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)
LIBS-$(BUILD_LIBVPX_SO) += $(BUILD_PFX)$(LIBVPX_SO)
$(BUILD_PFX)$(LIBVPX_SO): $(LIBVPX_OBJS) libvpx.ver
$(BUILD_PFX)$(LIBVPX_SO): extralibs += -lm -pthread
$(BUILD_PFX)$(LIBVPX_SO): SONAME = libvpx.so.$(VERSION_MAJOR)
$(BUILD_PFX)$(LIBVPX_SO): SO_VERSION_SCRIPT = libvpx.ver
LIBVPX_SO_SYMLINKS      := $(addprefix $(LIBSUBDIR)/, \
                             libvpx.so libvpx.so.$(VERSION_MAJOR) \
                             libvpx.so.$(VERSION_MAJOR).$(VERSION_MINOR))

libvpx.ver: $(call enabled,CODEC_EXPORTS)
	@echo "    [CREATE] $@"
	$(qexec)echo "{ global:" > $@
	$(qexec)for f in $?; do awk '{print $$2";"}' < $$f >>$@; done
	$(qexec)echo "local: *; };" >> $@
CLEAN-OBJS += libvpx.ver

$(addprefix $(DIST_DIR)/,$(LIBVPX_SO_SYMLINKS)):
	@echo "    [LN]      $@"
	$(qexec)ln -sf $(LIBVPX_SO) $@

INSTALL-LIBS-$(CONFIG_SHARED) += $(LIBVPX_SO_SYMLINKS)
INSTALL-LIBS-$(CONFIG_SHARED) += $(LIBSUBDIR)/$(LIBVPX_SO)

LIBS-$(BUILD_LIBVPX) += vpx.pc
vpx.pc: config.mk libs.mk
	@echo "    [CREATE] $@"
	$(qexec)echo '# pkg-config file from libvpx $(VERSION_STRING)' > $@
	$(qexec)echo 'prefix=$(PREFIX)' >> $@
	$(qexec)echo 'exec_prefix=$${prefix}' >> $@
	$(qexec)echo 'libdir=$${prefix}/lib' >> $@
	$(qexec)echo 'includedir=$${prefix}/include' >> $@
	$(qexec)echo '' >> $@
	$(qexec)echo 'Name: vpx' >> $@
	$(qexec)echo 'Description: WebM Project VPx codec implementation' >> $@
	$(qexec)echo 'Version: $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)' >> $@
	$(qexec)echo 'Requires:' >> $@
	$(qexec)echo 'Conflicts:' >> $@
	$(qexec)echo 'Libs: -L$${libdir} -lvpx' >> $@
	$(qexec)echo 'Cflags: -I$${includedir}' >> $@
INSTALL-LIBS-yes += $(LIBSUBDIR)/pkgconfig/vpx.pc
INSTALL_MAPS += $(LIBSUBDIR)/pkgconfig/%.pc %.pc
CLEAN-OBJS += vpx.pc
endif

LIBS-$(LIPO_LIBVPX) += libvpx.a
$(eval $(if $(LIPO_LIBVPX),$(call lipo_lib_template,libvpx.a)))

#
# Rule to make assembler configuration file from C configuration file
#
ifeq ($(ARCH_X86)$(ARCH_X86_64),yes)
# YASM
$(BUILD_PFX)vpx_config.asm: $(BUILD_PFX)vpx_config.h
	@echo "    [CREATE] $@"
	@egrep "#define [A-Z0-9_]+ [01]" $< \
	    | awk '{print $$2 " equ " $$3}' > $@
else
ADS2GAS=$(if $(filter yes,$(CONFIG_GCC)),| $(ASM_CONVERSION))
$(BUILD_PFX)vpx_config.asm: $(BUILD_PFX)vpx_config.h
	@echo "    [CREATE] $@"
	@egrep "#define [A-Z0-9_]+ [01]" $< \
	    | awk '{print $$2 " EQU " $$3}' $(ADS2GAS) > $@
	@echo "        END" $(ADS2GAS) >> $@
CLEAN-OBJS += $(BUILD_PFX)vpx_config.asm
endif

#
# Add assembler dependencies for configuration and offsets
#
$(filter %.s.o,$(OBJS-yes)):     $(BUILD_PFX)vpx_config.asm
$(filter %$(ASM).o,$(OBJS-yes)): $(BUILD_PFX)vpx_config.asm

#
# Calculate platform- and compiler-specific offsets for hand coded assembly
#
ifeq ($(CONFIG_EXTERNAL_BUILD),) # Visual Studio uses obj_int_extract.bat
  ifeq ($(ARCH_ARM), yes)
    asm_com_offsets.asm: obj_int_extract
    asm_com_offsets.asm: $(VP8_PREFIX)common/asm_com_offsets.c.o
	./obj_int_extract rvds $< $(ADS2GAS) > $@
    OBJS-yes += $(VP8_PREFIX)common/asm_com_offsets.c.o
    CLEAN-OBJS += asm_com_offsets.asm
    $(filter %$(ASM).o,$(OBJS-yes)): $(BUILD_PFX)asm_com_offsets.asm
  endif

  ifeq ($(ARCH_ARM)$(ARCH_X86)$(ARCH_X86_64), yes)
    ifeq ($(CONFIG_VP8_ENCODER), yes)
      asm_enc_offsets.asm: obj_int_extract
      asm_enc_offsets.asm: $(VP8_PREFIX)encoder/asm_enc_offsets.c.o
	./obj_int_extract rvds $< $(ADS2GAS) > $@
      OBJS-yes += $(VP8_PREFIX)encoder/asm_enc_offsets.c.o
      CLEAN-OBJS += asm_enc_offsets.asm
      $(filter %$(ASM).o,$(OBJS-yes)): $(BUILD_PFX)asm_enc_offsets.asm
    endif
  endif

  ifeq ($(ARCH_ARM), yes)
    ifeq ($(CONFIG_VP8_DECODER), yes)
      asm_dec_offsets.asm: obj_int_extract
      asm_dec_offsets.asm: $(VP8_PREFIX)decoder/asm_dec_offsets.c.o
	./obj_int_extract rvds $< $(ADS2GAS) > $@
      OBJS-yes += $(VP8_PREFIX)decoder/asm_dec_offsets.c.o
      CLEAN-OBJS += asm_dec_offsets.asm
      $(filter %$(ASM).o,$(OBJS-yes)): $(BUILD_PFX)asm_dec_offsets.asm
    endif
  endif
endif

$(shell $(SRC_PATH_BARE)/build/make/version.sh "$(SRC_PATH_BARE)" $(BUILD_PFX)vpx_version.h)
CLEAN-OBJS += $(BUILD_PFX)vpx_version.h

CODEC_DOC_SRCS += vpx/vpx_codec.h \
                  vpx/vpx_decoder.h \
                  vpx/vpx_encoder.h \
                  vpx/vpx_image.h

CLEAN-OBJS += libs.doxy
DOCS-yes += libs.doxy
libs.doxy: $(CODEC_DOC_SRCS)
	@echo "    [CREATE] $@"
	@rm -f $@
	@echo "INPUT += $^" >> $@
	@echo "PREDEFINED = VPX_CODEC_DISABLE_COMPAT" >> $@
	@echo "INCLUDE_PATH += ." >> $@;
	@echo "ENABLED_SECTIONS += $(sort $(CODEC_DOC_SECTIONS))" >> $@
