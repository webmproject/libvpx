/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VPX_FRAME_BUFFER_H_
#define VPX_VPX_FRAME_BUFFER_H_

/*!\file
 * \brief Describes the decoder external frame buffer interface.
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "./vpx_image.h"
#include "./vpx_integer.h"

/*!\brief The maximum number of work buffers used by libvpx.
 *  Support maximum 4 threads to decode video in parallel.
 *  Each thread will use one work buffer.
 * TODO(hkuang): Add support to set number of worker threads dynamically.
 */
#define VPX_MAXIMUM_WORK_BUFFERS 8

/*!\brief The maximum number of reference buffers that a VP9 encoder may use.
 */
#define VP9_MAXIMUM_REF_BUFFERS 8

/*!\brief The maximum number of planes that can be allocated by an external
 * frame buffer.
 */
#define VP9_MAXIMUM_EXTERNAL_PLANES 3

/*!\brief External frame buffer allocation type.
 * This enum defines if the external frame buffer contains one allocation for
 * all the planes contained in the frame buffer, or one allocation per plane.
 */
typedef enum vpx_codec_frame_buffer_type {
  VPX_CODEC_FRAME_BUFFER_TYPE_SIZE = 1,
  VPX_CODEC_FRAME_BUFFER_TYPE_PLANES = 2
} vpx_codec_frame_buffer_type_t;

/*!\brief External frame buffer
 *
 * This structure holds allocated frame buffers used by the decoder.
 */
typedef struct vpx_codec_frame_buffer {
  vpx_img_fmt_t fmt;  /**< Requested framebuffer format. */
  size_t width;  /**< Logical width of the frame. */
  size_t height;  /**< Logical height of the frame. */

  /*! Whether data and size, or plane and stride should be used. */
  vpx_codec_frame_buffer_type_t type;

  uint8_t *data;  /**< Pointer to the data buffer */
  size_t size;  /**< Size of data in bytes */

  /*! Pointers for each plane buffer */
  uint8_t *plane[VP9_MAXIMUM_EXTERNAL_PLANES];
  int stride[VP9_MAXIMUM_EXTERNAL_PLANES]; /**< Strides for each plane */

  void *priv;  /**< Frame's private data */
} vpx_codec_frame_buffer_t;

/*!\brief get frame buffer callback prototype
 *
 * This callback is invoked by the decoder to allocate storage for the frame
 * buffer in order for the decode call to complete.
 * The callback can allocate one buffer for all the frame buffers, or one buffer
 * per plane with the number and size of planes determined by fb->fmt.
 * If the application decides to allocate only one buffer, fb->type must be set
 * to #VPX_CODEC_FRAME_BUFFER_SIZE, the buffer size must be at least |min_size|
 * in bytes and the pointer to the buffer must be assigned to fb->data.
 * Then the callback must set fb->size to the allocated size.
 * The application does not need to align the allocated data.
 * If the application decides to allocate one buffer per plane, fb->type must be
 * set to #VPX_CODEC_FRAME_BUFFER_PLANES, the buffer pointers must be assigned
 * to fb->plane, and their strides to fb->stride.
 * The application provided buffers and strides must be aligned to 32 bytes.
 * The callback must zero out all the buffers allocated.
 *
 * The callback is triggered when the decoder needs a frame buffer to
 * decode a compressed image into. This function may be called more than once
 * for every call to vpx_codec_decode. The application may set fb->priv to
 * some data which will be passed back in the ximage and the release function
 * call. |fb| is guaranteed to not be NULL. On success the callback must
 * return 0. Any failure the callback must return a value less than 0.
 *
 * \param[in] priv         Callback's private data
 * \param[in] new_size     Size in bytes needed by the buffer
 * \param[in,out] fb       Pointer to vpx_codec_frame_buffer_t
 */
typedef int (*vpx_get_frame_buffer_cb_fn_t)(
    void *priv, size_t min_size, vpx_codec_frame_buffer_t *fb);

/*!\brief release frame buffer callback prototype
 *
 * This callback is invoked by the decoder when the frame buffer is not
 * referenced by any other buffers. |fb| is guaranteed to not be NULL. On
 * success the callback must return 0. Any failure the callback must return
 * a value less than 0.
 *
 * \param[in] priv         Callback's private data
 * \param[in] fb           Pointer to vpx_codec_frame_buffer_t
 */
typedef int (*vpx_release_frame_buffer_cb_fn_t)(
    void *priv, vpx_codec_frame_buffer_t *fb);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VPX_VPX_FRAME_BUFFER_H_
