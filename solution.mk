##
##  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##


# libaom reverse dependencies (targets that depend on libaom)
AOM_NONDEPS=$(addsuffix .$(VCPROJ_SFX),aom gtest)
AOM_RDEPS=$(foreach vcp,\
              $(filter-out $(AOM_NONDEPS),$^), --dep=$(vcp:.$(VCPROJ_SFX)=):aom)

aom.sln: $(wildcard *.$(VCPROJ_SFX))
	@echo "    [CREATE] $@"
	$(SRC_PATH_BARE)/build/make/gen_msvs_sln.sh \
            $(if $(filter aom.$(VCPROJ_SFX),$^),$(AOM_RDEPS)) \
            --dep=test_libaom:gtest \
            --ver=$(CONFIG_VS_VERSION)\
            --out=$@ $^
aom.sln.mk: aom.sln
	@true

PROJECTS-yes += aom.sln aom.sln.mk
-include aom.sln.mk

# Always install this file, as it is an unconditional post-build rule.
INSTALL_MAPS += src/%     $(SRC_PATH_BARE)/%
INSTALL-SRCS-yes            += $(target).mk
