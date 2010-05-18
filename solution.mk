##
##  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license and patent
##  grant that can be found in the LICENSE file in the root of the source
##  tree. All contributing project authors may be found in the AUTHORS
##  file in the root of the source tree.
##


ifeq ($(ARCH_ARM),yes)
ARM_DEVELOP=no
ARM_DEVELOP:=$(if $(filter %vpx.vcproj,$(wildcard *.vcproj)),yes)

ifeq ($(ARM_DEVELOP),yes)
vpx.sln:
	@echo "    [COPY] $@"
	@cp $(SRC_PATH_BARE)/build/arm-wince-vs8/vpx.sln .
PROJECTS-yes += vpx.sln
else
vpx.sln: $(wildcard *.vcproj)
	@echo "    [CREATE] $@"
	$(SRC_PATH_BARE)/build/make/gen_msvs_sln.sh \
            $(if $(filter %vpx.vcproj,$^),--dep=ivfdec:vpx) \
            $(if $(filter %vpx.vcproj,$^),--dep=xma:vpx) \
            --ver=$(CONFIG_VS_VERSION)\
            --target=$(TOOLCHAIN)\
            --out=$@ $^
vpx.sln.mk: vpx.sln
	@true

PROJECTS-yes += vpx.sln vpx.sln.mk
-include vpx.sln.mk
endif

else
vpx.sln: $(wildcard *.vcproj)
	@echo "    [CREATE] $@"
	$(SRC_PATH_BARE)/build/make/gen_msvs_sln.sh \
            $(if $(filter %vpx.vcproj,$^),\
                $(foreach vcp,$(filter-out %vpx.vcproj,$^),\
                  --dep=$(vcp:.vcproj=):vpx)) \
            --ver=$(CONFIG_VS_VERSION)\
            --out=$@ $^
vpx.sln.mk: vpx.sln
	@true

PROJECTS-yes += vpx.sln vpx.sln.mk
-include vpx.sln.mk
endif

# Always install this file, as it is an unconditional post-build rule.
INSTALL_MAPS += src/%     $(SRC_PATH_BARE)/%
INSTALL-SRCS-yes            += $(target).mk
