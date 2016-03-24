// Copyright (c) 2016 The WebM project authors. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the LICENSE file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS.  All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
#include "hdr_util.h"

#include <cstddef>
#include <new>

#include "mkvparser/mkvparser.h"

namespace libwebm {
bool CopyPrimaryChromaticity(const mkvparser::PrimaryChromaticity& parser_pc,
                             PrimaryChromaticityPtr* muxer_pc) {
  muxer_pc->reset(new (std::nothrow)
                      mkvmuxer::PrimaryChromaticity(parser_pc.x, parser_pc.y));
  if (!muxer_pc->get())
    return false;
  return true;
}

bool MasteringMetadataValuePresent(double value) {
  return value != mkvparser::MasteringMetadata::kValueNotPresent;
}

bool CopyMasteringMetadata(const mkvparser::MasteringMetadata& parser_mm,
                           mkvmuxer::MasteringMetadata* muxer_mm) {
  if (MasteringMetadataValuePresent(parser_mm.luminance_max))
    muxer_mm->luminance_max = parser_mm.luminance_max;
  if (MasteringMetadataValuePresent(parser_mm.luminance_min))
    muxer_mm->luminance_min = parser_mm.luminance_min;

  PrimaryChromaticityPtr r_ptr(NULL);
  PrimaryChromaticityPtr g_ptr(NULL);
  PrimaryChromaticityPtr b_ptr(NULL);
  PrimaryChromaticityPtr wp_ptr(NULL);

  if (parser_mm.r) {
    if (!CopyPrimaryChromaticity(*parser_mm.r, &r_ptr))
      return false;
  }
  if (parser_mm.g) {
    if (!CopyPrimaryChromaticity(*parser_mm.g, &g_ptr))
      return false;
  }
  if (parser_mm.b) {
    if (!CopyPrimaryChromaticity(*parser_mm.b, &b_ptr))
      return false;
  }
  if (parser_mm.white_point) {
    if (!CopyPrimaryChromaticity(*parser_mm.white_point, &wp_ptr))
      return false;
  }

  if (!muxer_mm->SetChromaticity(r_ptr.get(), g_ptr.get(), b_ptr.get(),
                                 wp_ptr.get())) {
    return false;
  }

  return true;
}

bool ColourValuePresent(long long value) {
  return value != mkvparser::Colour::kValueNotPresent;
}

bool CopyColour(const mkvparser::Colour& parser_colour,
                mkvmuxer::Colour* muxer_colour) {
  if (!muxer_colour)
    return false;

  if (ColourValuePresent(parser_colour.matrix_coefficients))
    muxer_colour->matrix_coefficients = parser_colour.matrix_coefficients;
  if (ColourValuePresent(parser_colour.bits_per_channel))
    muxer_colour->bits_per_channel = parser_colour.bits_per_channel;
  if (ColourValuePresent(parser_colour.chroma_subsampling_horz))
    muxer_colour->chroma_subsampling_horz =
        parser_colour.chroma_subsampling_horz;
  if (ColourValuePresent(parser_colour.chroma_subsampling_vert))
    muxer_colour->chroma_subsampling_vert =
        parser_colour.chroma_subsampling_vert;
  if (ColourValuePresent(parser_colour.cb_subsampling_horz))
    muxer_colour->cb_subsampling_horz = parser_colour.cb_subsampling_horz;
  if (ColourValuePresent(parser_colour.cb_subsampling_vert))
    muxer_colour->cb_subsampling_vert = parser_colour.cb_subsampling_vert;
  if (ColourValuePresent(parser_colour.chroma_siting_horz))
    muxer_colour->chroma_siting_horz = parser_colour.chroma_siting_horz;
  if (ColourValuePresent(parser_colour.chroma_siting_vert))
    muxer_colour->chroma_siting_vert = parser_colour.chroma_siting_vert;
  if (ColourValuePresent(parser_colour.range))
    muxer_colour->range = parser_colour.range;
  if (ColourValuePresent(parser_colour.transfer_characteristics))
    muxer_colour->transfer_characteristics =
        parser_colour.transfer_characteristics;
  if (ColourValuePresent(parser_colour.primaries))
    muxer_colour->primaries = parser_colour.primaries;
  if (ColourValuePresent(parser_colour.max_cll))
    muxer_colour->max_cll = parser_colour.max_cll;
  if (ColourValuePresent(parser_colour.max_fall))
    muxer_colour->max_fall = parser_colour.max_fall;

  if (parser_colour.mastering_metadata) {
    mkvmuxer::MasteringMetadata muxer_mm;
    if (!CopyMasteringMetadata(*parser_colour.mastering_metadata, &muxer_mm))
      return false;
    if (!muxer_colour->SetMasteringMetadata(muxer_mm))
      return false;
  }
  return true;
}
}  // namespace libwebm
