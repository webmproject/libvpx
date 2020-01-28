/*
 *  Copyright (c) 2019 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VP9_SIMPLE_ENCODE_H_
#define VPX_VP9_SIMPLE_ENCODE_H_

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <vector>

namespace vp9 {

enum FrameType {
  kKeyFrame = 0,
  kInterFrame,
  kAlternateReference,
};

// The frame is split to 4x4 blocks.
// This structure contains the information of each 4x4 block.
struct PartitionInfo {
  int row;           // row pixel offset of current 4x4 block
  int column;        // column pixel offset of current 4x4 block
  int row_start;     // row pixel offset of the start of the prediction block
  int column_start;  // column pixel offset of the start of the prediction block
  int width;         // prediction block width
  int height;        // prediction block height
};

struct EncodeFrameInfo {
  int show_idx;
  FrameType frame_type;
};

// This structure is a copy of vp9 |nmv_component_counts|.
struct NewMotionvectorComponentCounts {
  std::vector<unsigned int> sign;
  std::vector<unsigned int> classes;
  std::vector<unsigned int> class0;
  std::vector<std::vector<unsigned int>> bits;
  std::vector<std::vector<unsigned int>> class0_fp;
  std::vector<unsigned int> fp;
  std::vector<unsigned int> class0_hp;
  std::vector<unsigned int> hp;
};

// This structure is a copy of vp9 |nmv_context_counts|.
struct NewMotionVectorContextCounts {
  std::vector<unsigned int> joints;
  std::vector<NewMotionvectorComponentCounts> comps;
};

#define UintArray2D std::vector<std::vector<unsigned int>>
#define UintArray3D std::vector<std::vector<std::vector<unsigned int>>>
#define UintArray5D \
  std::vector<std::vector<std::vector<std::vector<std::vector<unsigned int>>>>>
#define UintArray6D        \
  std::vector<std::vector< \
      std::vector<std::vector<std::vector<std::vector<unsigned int>>>>>>

// This structure is a copy of vp9 |tx_counts|.
struct TransformSizeCounts {
  // Transform size found in blocks of partition size 32x32.
  // First dimension: transform size contexts (2).
  // Second dimension: transform size type (3: 32x32, 16x16, 8x8)
  UintArray2D p32x32;
  // Transform size found in blocks of partition size 16x16.
  // First dimension: transform size contexts (2).
  // Second dimension: transform size type (2: 16x16, 8x8)
  UintArray2D p16x16;
  // Transform size found in blocks of partition size 8x8.
  // First dimension: transform size contexts (2).
  // Second dimension: transform size type (1: 8x8)
  UintArray2D p8x8;
  // Overall transform size count.
  std::vector<unsigned int> tx_totals;
};

// This structure is a copy of vp9 |FRAME_COUNTS|.
struct FrameCounts {
  // Intra prediction mode for luma plane. First dimension: block size (4).
  // Second dimension: intra prediction mode (10).
  UintArray2D y_mode;
  // Intra prediction mode for chroma plane. First and second dimension:
  // intra prediction mode (10).
  UintArray2D uv_mode;
  // Partition type. First dimension: partition contexts (16).
  // Second dimension: partition type (4).
  UintArray2D partition;
  // Transform coefficient.
  UintArray6D coef;
  // End of block (the position of the last non-zero transform coefficient)
  UintArray5D eob_branch;
  // Interpolation filter type. First dimension: switchable filter contexts (4).
  // Second dimension: filter types (3).
  UintArray2D switchable_interp;
  // Inter prediction mode (the motion vector type).
  // First dimension: inter mode contexts (7).
  // Second dimension: mode type (4).
  UintArray2D inter_mode;
  // Block is intra or inter predicted. First dimension: contexts (4).
  // Second dimension: type (0 for intra, 1 for inter).
  UintArray2D intra_inter;
  // Block is compound predicted (predicted from average of two blocks).
  // First dimension: contexts (5).
  // Second dimension: type (0 for single, 1 for compound prediction).
  UintArray2D comp_inter;
  // Type of the reference frame. Only one reference frame.
  // First dimension: context (5). Second dimension: context (2).
  // Third dimension: count (2).
  UintArray3D single_ref;
  // Type of the two reference frames.
  // First dimension: context (5). Second dimension: count (2).
  UintArray2D comp_ref;
  // Block skips transform and quantization, uses prediction as reconstruction.
  // First dimension: contexts (3). Second dimension: type (0 not skip, 1 skip).
  UintArray2D skip;
  // Transform size.
  TransformSizeCounts tx;
  // New motion vector.
  NewMotionVectorContextCounts mv;
};

struct EncodeFrameResult {
  int show_idx;
  FrameType frame_type;
  size_t coding_data_bit_size;
  size_t coding_data_byte_size;
  // The EncodeFrame will allocate a buffer, write the coding data into the
  // buffer and give the ownership of the buffer to coding_data.
  std::unique_ptr<unsigned char[]> coding_data;
  double psnr;
  uint64_t sse;
  int quantize_index;
  FrameCounts frame_counts;
  int num_rows_4x4;  // number of row units, in size of 4.
  int num_cols_4x4;  // number of column units, in size of 4.
  // A vector of the partition information of the frame.
  // The number of elements is |num_rows_4x4| * |num_cols_4x4|.
  // The frame is divided 4x4 blocks of |num_rows_4x4| rows and
  // |num_cols_4x4| columns.
  // Each 4x4 block contains the current pixel position (|row|, |column|),
  // the start pixel position of the partition (|row_start|, |column_start|),
  // and the |width|, |height| of the partition.
  // The current pixel position can be the same as the start pixel position
  // if the 4x4 block is the top-left block in the partition. Otherwise, they
  // are different.
  // Within the same partition, all 4x4 blocks have the same |row_start|,
  // |column_start|, |width| and |height|.
  // For example, if the frame is partitioned to a 32x32 block,
  // starting at (0, 0). Then, there're 64 4x4 blocks within this partition.
  // They all have the same |row_start|, |column_start|, |width|, |height|,
  // which can be used to figure out the start of the current partition and
  // the start of the next partition block.
  // Horizontal next: |column_start| + |width|,
  // Vertical next: |row_start| + |height|.
  std::vector<PartitionInfo> partition_info;
};

struct GroupOfPicture {
  // This list will be updated internally in StartEncode() and
  // EncodeFrame()/EncodeFrameWithQuantizeIndex().
  // In EncodeFrame()/EncodeFrameWithQuantizeIndex(), the update will only be
  // triggered when the coded frame is the last one in the previous group of
  // pictures.
  std::vector<EncodeFrameInfo> encode_frame_list;
  // Indicates the index of the next coding frame in encode_frame_list.
  // In other words, EncodeFrameInfo of the next coding frame can be
  // obtained with encode_frame_list[next_encode_frame_index].
  // Internally, next_encode_frame_index will be set to zero after the last
  // frame of the group of pictures is coded. Otherwise, next_encode_frame_index
  // will be increased after each EncodeFrame()/EncodeFrameWithQuantizeIndex()
  // call.
  int next_encode_frame_index;
  // Number of show frames in this group of pictures.
  int show_frame_count;
  // The show index/timestamp of the earliest show frame in the group of
  // pictures.
  int start_show_index;
};

class SimpleEncode {
 public:
  SimpleEncode(int frame_width, int frame_height, int frame_rate_num,
               int frame_rate_den, int target_bitrate, int num_frames,
               const char *infile_path);
  ~SimpleEncode();
  SimpleEncode(SimpleEncode &) = delete;
  SimpleEncode &operator=(const SimpleEncode &) = delete;

  // Makes encoder compute the first pass stats and store it internally for
  // future encode.
  void ComputeFirstPassStats();

  // Outputs the first pass stats represented by a 2-D vector.
  // One can use the frame index at first dimension to retrieve the stats for
  // each video frame. The stats of each video frame is a vector of 25 double
  // values. For details, please check FIRSTPASS_STATS in vp9_firstpass.h
  std::vector<std::vector<double>> ObserveFirstPassStats();

  // Initializes the encoder for actual encoding.
  // This function should be called after ComputeFirstPassStats().
  void StartEncode();

  // Frees the encoder.
  // This function should be called after StartEncode() or EncodeFrame().
  void EndEncode();

  // Given a key_frame_index, computes this key frame group's size.
  // The key frame group size includes one key frame plus the number of
  // following inter frames. Note that the key frame group size only counts the
  // show frames. The number of no show frames like alternate refereces are not
  // counted.
  int GetKeyFrameGroupSize(int key_frame_index) const;

  // Provides the group of pictures that the next coding frame is in.
  // Only call this function between StartEncode() and EndEncode()
  GroupOfPicture ObserveGroupOfPicture() const;

  // Gets encode_frame_info for the next coding frame.
  // Only call this function between StartEncode() and EndEncode()
  EncodeFrameInfo GetNextEncodeFrameInfo() const;

  // Encodes a frame
  // This function should be called after StartEncode() and before EndEncode().
  void EncodeFrame(EncodeFrameResult *encode_frame_result);

  // Encodes a frame with a specific quantize index.
  // This function should be called after StartEncode() and before EndEncode().
  void EncodeFrameWithQuantizeIndex(EncodeFrameResult *encode_frame_result,
                                    int quantize_index);

  // Gets the number of coding frames for the video. The coding frames include
  // show frame and no show frame.
  // This function should be called after ComputeFirstPassStats().
  int GetCodingFrameNum() const;

  // Gets the total number of pixels of YUV planes per frame.
  uint64_t GetFramePixelCount() const;

 private:
  class EncodeImpl;

  int frame_width_;   // frame width in pixels.
  int frame_height_;  // frame height in pixels.
  int num_rows_4x4_;  // number of row units, in size of 4.
  int num_cols_4x4_;  // number of column units, in size of 4.
  int frame_rate_num_;
  int frame_rate_den_;
  int target_bitrate_;
  int num_frames_;
  std::FILE *file_;
  std::unique_ptr<EncodeImpl> impl_ptr_;

  GroupOfPicture group_of_picture_;
};

}  // namespace vp9

#endif  // VPX_VP9_SIMPLE_ENCODE_H_
