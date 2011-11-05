##
##  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##


vpx.sln: $(wildcard *.vcproj)
	@echo "    [CREATE] $@"
	$(SRC_PATH_BARE)/build/make/gen_msvs_sln.sh \
            $(if $(filter %vpx.vcproj,$^),\
                $(foreach vcp,$(filter-out %vpx.vcproj %gtest.vcproj %obj_int_extract.vcproj,$^),\
                  --dep=$(vcp:.vcproj=):vpx) \
                $(foreach vcp,$(filter %_test.vcproj,$^),\
                  --dep=$(vcp:.vcproj=):gtest)) \
                  --dep=vpx:obj_int_extract \
                  --ver=$(CONFIG_VS_VERSION)\
                  --out=$@ $^
vpx.sln.mk: vpx.sln
	@true

PROJECTS-yes += vpx.sln vpx.sln.mk
-include vpx.sln.mk

# Always install this file, as it is an unconditional post-build rule.
INSTALL_MAPS += src/%     $(SRC_PATH_BARE)/%
INSTALL-SRCS-yes            += $(target).mk
