/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


typedef struct core_codec *codec_ptr;
typedef struct interface_table *interface_ptr;

typedef struct
{
    void (*Initialize)();
    void (*Shutdown)();
    codec_ptr(*Create)();
    int (*compress_frame)(codec_ptr, unsigned int *frame_flags, YV12_BUFFER_CONFIG *sd, unsigned long *size, char *dest, INT64 time_stamp);
    int (*show_frame)(codec_ptr , YV12_BUFFER_CONFIG *dest, int deblock_level, int noise_level, int flags);
    void (*Remove)(codec_ptr *comp);
    interface_ptr(*get_interface)(unsigned int id);

} core_codec;

typedef struct
{
    int (*set_bitrate)(codec_ptr, END_USAGE usage, int Datarate);
    int (*get_bitrate)(codec_ptr, END_USAGE *usage, int *Datarate);
    int (*set_mode)(codec_ptr, MODE mode, int Speed, char *File);
    int (*get_mode)(codec_ptr, MODE *mode, int *Speed, char **File);
} codec_settings_basic;

typedef struct
{
    int (*set_bitrate)(codec_ptr, END_USAGE usage, int Datarate);
    int (*get_bitrate)(codec_ptr, END_USAGE *usage, int *Datarate);
    int (*set_mode)(codec_ptr, MODE  mode, int  Speed, char   *File);
    int (*get_mode)(codec_ptr, MODE *mode, int *Speed, char **File);
    int (*set_denoise)(codec_ptr, int  Level);
    int (*get_denoise)(codec_ptr, int *Level);
    int (*set_sharpness)(codec_ptr, int  sharpness);
    int (*get_sharpness)(codec_ptr, int *sharpness);
    int (*set_keyframing)(codec_ptr, int  Auto, int  max_distance);
    int (*get_keyframing)(codec_ptr, int *Auto, int *max_distance);
    int (*set_buffering)(codec_ptr, int  buffer_level, int  max_buffer_level);
    int (*get_buffering)(codec_ptr, int *buffer_level, int *max_buffer_level);
    int (*set_adjust_frame_rate)(codec_ptr, int Allowed, int at_buffer_level_pct);
    int (*get_adjust_frame_rate)(codec_ptr, int *Allowed, int *at_buffer_level_pct);
    int (*set_adjust_frame_size)(codec_ptr, int Allowed, int down_at_buffer_level_pct, int up_at_buffer_level_pct);
    int (*get_adjust_frame_size)(codec_ptr, int *Allowed, int *down_at_buffer_level_pct, int *up_at_buffer_level_pct);
    int (*set_adjust_quality)(codec_ptr, int Allowed, int min_quantizer, int max_quantizer);
    int (*get_adjust_quality)(codec_ptr, int *Allowed, int *min_quantizer, int *max_quantizer);
    int (*set_vbrparms)(codec_ptr, int Bias, int Min, int Max);
    int (*get_vbrparms)(codec_ptr, int *Bias, int *Min, int *Max);

} codec_settings_v1;

typedef struct
{
    int (*request_recovery)(codec_ptr);
    int (*request_droppable)(codec_ptr);
    int (*internal_size)(codec_ptr, VPX_SCALING Vertical, VPX_SCALING Horizontal);
    int (*update_last)(codec_ptr);
    int (*update_gold)(codec_ptr);
    int (*use_only_last)(codec_ptr);
    int (*use_only_gold)(codec_ptr);
    int (*update_entropy)(codec_ptr);

} codec_realtime_requests;
