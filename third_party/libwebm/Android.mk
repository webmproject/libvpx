LOCAL_PATH:= $(call my-dir)

include $(CLEAR_VARS)
LOCAL_MODULE:= libwebm
LOCAL_SRC_FILES:= common/hdr_util.cc \
                  mkvparser/mkvparser.cc \
                  mkvparser/mkvreader.cc \
                  mkvmuxer/mkvmuxer.cc \
                  mkvmuxer/mkvmuxerutil.cc \
                  mkvmuxer/mkvwriter.cc
include $(BUILD_STATIC_LIBRARY)
