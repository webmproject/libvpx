# Copyright (c) 2013 The WebM project authors. All Rights Reserved.
#
# Use of this source code is governed by a BSD-style license
# that can be found in the LICENSE file in the root of the source
# tree. An additional intellectual property rights grant can be found
# in the file PATENTS.  All contributing project authors may
# be found in the AUTHORS file in the root of the source tree.
#
# This make file builds vpx_test app for android.
# The test app itself runs on the command line through adb shell
# The paths are really messed up as the libvpx make file
# expects to be made from a parent directory.
# TODO(joshualitt)
# Fix android make files so they can be built from anywhere, will require
# changing the libvpx make file and this one.
CUR_WD := $(call my-dir)
BINDINGS_DIR := $(CUR_WD)/../../..
LOCAL_PATH := $(CUR_WD)/../../..

#libvpx
include $(CLEAR_VARS)
include $(BINDINGS_DIR)/libvpx/build/make/Android.mk
# Restore path
# TODO joshualitt Fix makefiles so this is no longer needed
LOCAL_PATH := $(CUR_WD)/../..

#libgtest
include $(CLEAR_VARS)
LOCAL_CPP_EXTENSION := .cc
LOCAL_MODULE := gtest
LOCAL_C_INCLUDES := $(LOCAL_PATH)/third_party/googletest/src/
LOCAL_C_INCLUDES += $(LOCAL_PATH)/third_party/googletest/src/include/
LOCAL_SRC_FILES := ./third_party/googletest/src/src/gtest-all.cc
include $(BUILD_STATIC_LIBRARY)

#libnestegg
include $(CLEAR_VARS)
LOCAL_CPP_EXTENSION := .cc
LOCAL_MODULE := nestegg
NESTEGG_PATH := $(LOCAL_PATH)/nestegg
LOCAL_C_INCLUDES := $(NESTEGG_PATH)/include
LOCAL_C_INCLUDES += $(LOCAL_PATH)/
LOCAL_C_INCLUDES += $(NESTEGG_PATH)/halloc/
LOCAL_SRC_FILES := ./nestegg/halloc/src/halloc.c
LOCAL_SRC_FILES += ./nestegg/src/nestegg.c
include $(BUILD_STATIC_LIBRARY)

#libvpx_test
include $(CLEAR_VARS)
LOCAL_MODULE := libvpx_test
LOCAL_STATIC_LIBRARIES := gtest
LOCAL_STATIC_LIBRARIES += nestegg
LOCAL_STATIC_LIBRARIES += cpufeatures
LOCAL_SHARED_LIBRARIES := vpx
LOCAL_C_INCLUDES := $(LOCAL_PATH)/
LOCAL_C_INCLUDES += $(BINDINGS_DIR)/
LOCAL_C_INCLUDES += $(LOCAL_PATH)/third_party/googletest/src/include
LOCAL_SRC_FILES := ./args.c
LOCAL_SRC_FILES += ./md5_utils.c
LOCAL_SRC_FILES += ./test/decode_test_driver.cc
LOCAL_SRC_FILES += ./test/test_libvpx.cc
LOCAL_SRC_FILES += ./test/test_vector_test.cc
include $(BUILD_EXECUTABLE)
