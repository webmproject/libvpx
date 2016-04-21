/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be
 *  found  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>

#include "vp9/common/vp9_motion_model.h"
#include "vp9/encoder/vp9_opticalflow.h"
#include "vp9/encoder/vp9_ransac.h"
#include "vp9/encoder/vp9_resize.h"

// ITERATIONS is the number of iterative refinements to make on an initial
// flow estimate. 3 seems to work well but feel free to play with this.
#define ITERATIONS 3
#define MAX_MEDIAN_LENGTH  121
#define MAX_LEVELS 8
#define MIN_LEVEL_SIZE 30
#define MAX_ERROR  0.00001
#define BLUR_SIZE 11
#define BLUR_SIGMA 1.5
#define BLUR_AMP 1
#define MEDIAN_FILTER_SIZE 5
#define RELATIVE_ERROR(a, b) (fabs(((a) - (b))))
#define WRAP_PIXEL(p, size) ((p) < 0 ? abs((p)) - 1 : \
                            ((p) >= (size) ? (2 * (size) - 1 - (p)) : (p)))

// Struct for an image pyramid
typedef struct {
  int n_levels;
  int widths[MAX_LEVELS];
  int heights[MAX_LEVELS];
  int strides[MAX_LEVELS];
  unsigned char **levels;
} imPyramid;

typedef struct {
  int stride;
  double *confidence;
  double *u;
  double *v;
} flowMV;

static void * safe_malloc(size_t size) {
  void *v = malloc(size);
  if (!v) {
    fprintf(stderr, "Not enough memory\n");
    exit(EXIT_FAILURE);
  }
  return v;
}

// Produces an image pyramid where each level is resized by a
// specified resize factor
static void image_pyramid(unsigned char *img, imPyramid *pyr, int width,
                          int height, int stride, int n_levels) {
  int i;
  int max_levels = 1;
  int scaled = width < height ? (width/MIN_LEVEL_SIZE) :
               (height/MIN_LEVEL_SIZE);
  pyr->widths[0] = width;
  pyr->heights[0] = height;
  pyr->strides[0] = stride;

  if (n_levels == 1) {
    pyr->n_levels = 1;
    pyr->levels = (unsigned char**)safe_malloc(sizeof(*pyr->levels));
    pyr->levels[0] = img;
    return;
  }
  // compute number of possible levels of the pyramid. The smallest dim of the
  // smallest level should be no less than 30 px.
  while (scaled > 1) {
    max_levels++;
    scaled >>= 1;
  }

  if (n_levels > max_levels)
    n_levels = max_levels;

  pyr->n_levels = n_levels;
  pyr->levels = (unsigned char**)safe_malloc(n_levels * sizeof(*pyr->levels));
  pyr->levels[0] = img;

  // compute each level
  for (i = 1; i < n_levels; ++i) {
    pyr->widths[i] = pyr->widths[i-1] >> 1;
    pyr->heights[i] = pyr->heights[i-1] >> 1;
    pyr->strides[i] = pyr->widths[i];
    pyr->levels[i] = (unsigned char*)safe_malloc(pyr->widths[i] *
                                                 pyr->heights[i] *
                                                 sizeof(*pyr->levels[i]));

    vp9_resize_plane(pyr->levels[i-1], pyr->heights[i-1], pyr->widths[i-1],
                     pyr->strides[i-1], pyr->levels[i], pyr->heights[i],
                     pyr->widths[i], pyr->strides[i]);
  }
}

static void destruct_pyramid(imPyramid *pyr) {
  int i;

  // start at level 1 because level 0 is a reference to the original image
  for (i = 1; i < pyr->n_levels; ++i)
    free(pyr->levels[i]);
  free(pyr->levels);
}

// Convolution with reflected content padding to minimize edge artifacts.
static void convolve_char(double *filter, unsigned char *image,
                          double *convolved, int filt_width, int filt_height,
                          int image_width, int image_height, int image_stride,
                          int convolved_stride) {
  int i, j, x, y, imx, imy;
  double pixel_val;
  double filter_val;
  double image_val;
  int x_pad_size = filt_width >> 1;
  int y_pad_size = filt_height >> 1;
  for (j = 0; j < image_height; ++j)
    for (i = 0; i < image_width; ++i) {
      pixel_val = 0;
      for (y = (-y_pad_size); y<= y_pad_size; y++)
        for (x = (-x_pad_size); x <= x_pad_size; x++) {
          filter_val = filter[(x + x_pad_size) + (y + y_pad_size) * filt_width];
          imx = WRAP_PIXEL(i + x, image_width);
          imy = WRAP_PIXEL(j + y, image_height);
          image_val = image[imx + imy * image_stride];
          pixel_val += (filter_val * image_val);
        }
      convolved[i + j * convolved_stride] =  pixel_val;
    }
}

// Convolution with reflected content padding to minimize edge artifacts.
static void convolve_double(const double *filter, double *image,
                            double *convolved, int filt_width,
                            int filt_height, int image_width,
                            int image_height, int image_stride,
                            int convolved_stride) {
  int i, j, x, y, imx, imy;
  double pixel_val;
  double filter_val;
  double image_val;
  int x_pad_size = filt_width >> 1;
  int y_pad_size = filt_height >> 1;
  for (j = 0; j < image_height; ++j)
    for (i = 0; i < image_width; ++i) {
      pixel_val = 0;
      for (y = (-y_pad_size); y<= y_pad_size; y++)
        for (x = (-x_pad_size); x <= x_pad_size; x++) {
          filter_val = filter[(x + x_pad_size) + (y + y_pad_size) * filt_width];
          imx = WRAP_PIXEL(i + x, image_width);
          imy = WRAP_PIXEL(j + y, image_height);
          image_val = image[imx + imy * image_stride];
          pixel_val += (filter_val * image_val);
        }
      convolved[i + j * convolved_stride] = pixel_val;
    }
}

// computes x and y spatial derivatives of an image.
static void differentiate(double *img, int width,  int height,
                          int stride, double *dx_img,  double *dy_img,
                          int deriv_stride) {
  const double spatial_deriv_kernel[] = {-0.0833, 0.6667, 0.0,
                                         -0.6667, 0.0833};
  convolve_double(spatial_deriv_kernel, img, dx_img, 5, 1, width, height,
                  stride, deriv_stride);
  convolve_double(spatial_deriv_kernel, img, dy_img, 1, 5, width, height,
                  stride, deriv_stride);
}

// creates a 2D gaussian kernel
static void gaussian_kernel(double *mat, int size, double sig, double amp) {
  int x, y;
  double filter_val;
  int half;
  double denominator = 1 / (2.0 * sig * sig);
  double sum = 0.0f;

  size |= 1;
  half = size >> 1;
  for (y = -half; y<= half; y++)
    for (x = -half; x<= half; x++) {
      filter_val = amp * exp((double)(-1.0 * ((x * x + y * y) * denominator)));
      mat[(x + half) + ((y + half) * size)] = filter_val;
      sum += filter_val;
    }

  // normalize the filter
  if (sum > MAX_ERROR) {
    sum = 1 / sum;
    for (x = 0; x < size * size; x++) {
      mat[x] = mat[x] * sum;
    }
  }
}

// blurs image with a gaussian kernel
static void blur_img(unsigned char *img, int width, int height,
                     int stride, double *smooth_img, int smooth_stride,
                     int smoothing_size, double smoothing_sig,
                     double smoothing_amp) {
  double *smoothing_kernel = (double *)safe_malloc(smoothing_size *
                                                   smoothing_size *
                                                   sizeof(*smoothing_kernel));
  gaussian_kernel(smoothing_kernel, smoothing_size, smoothing_sig,
                  smoothing_amp);
  convolve_char(smoothing_kernel, img, smooth_img, smoothing_size,
                smoothing_size, width, height, stride, smooth_stride);
  free(smoothing_kernel);
}

// Implementation of Hoare's select for linear time median selection.
// Takes in a pointer to the beginning of a list of values, the length
// of the list, an integer representing the index of the number to be
// selected, and an unsigned integer used to seed a random number generator.
static double selection(double *vals, int length, int k, unsigned int seed) {
  int pivot_ind, x;
  double pivot;
  double L[MAX_MEDIAN_LENGTH], G[MAX_MEDIAN_LENGTH];
  int l = 0, e = 0, g = 0;
  pivot_ind = (int) ((length / ((double)RAND_MAX)) * rand_r(&seed));
  pivot = vals[pivot_ind];
  for (x = 0; x < length; ++x) {
    if (RELATIVE_ERROR(vals[x], pivot) <= MAX_ERROR) {
      e++;
    } else if (vals[x] < pivot) {
      L[l] = vals[x];
      l++;
    } else if (vals[x] > pivot) {
      G[g] = vals[x];
      g++;
    }
  }

  if (k <= l)
    return selection(L, l, k, seed);
  else if (k <= l + e)
    return pivot;
  else
    return selection(G, g, k - l - e, seed);
}

// Performs median filtering for denoising. This implementation uses hoare's
// select to find the median in a block. block_x, block_y are
// the x and y coordinates of the center of the block in the source.
static void median_filter(double *source, double *filtered, int block_size,
                          int width, int height, int source_stride,
                          int filtered_stride) {
  int i, j, block_x, block_y, imx, imy;
  double pivot, val;
  int length = block_size * block_size;
  double L[MAX_MEDIAN_LENGTH], G[MAX_MEDIAN_LENGTH];
  int k = length >> 1 | 1;
  int l, e, g;
  int pad_size = block_size >> 1;
  unsigned int seed = (unsigned int)*source;

  // find the median within each block. Reflected content is used for padding.
  for (block_y = 0; block_y < height; ++block_y)
    for (block_x = 0; block_x < width; ++block_x) {
      l = 0, e = 0, g = 0;
      memset(L, 0, length * sizeof(*L));
      memset(G, 0, length * sizeof(*G));
      pivot = source[block_x + block_y * source_stride];
      for (j = -pad_size; j <= pad_size; ++j)
        for (i = -pad_size; i <= pad_size; ++i) {
          imx = WRAP_PIXEL(i + block_x, width);
          imy = WRAP_PIXEL(j + block_y, height);
          val = source[imx + imy * source_stride];

          // pulling out the first iteration of selection so we don't have
          // to iterate through the block to flatten it and then
          // iterate through those same values to put them into
          // the L E G bins in selection. Ugly but more efficent.
          if (RELATIVE_ERROR(val, pivot) <= MAX_ERROR) {
            e++;
          } else if (val < pivot) {
            L[l] = val;
            l++;
          } else if (val > pivot) {
            G[g] = val;
            g++;
          }
        }
      if (k <= l)
        filtered[block_x + block_y * filtered_stride] =
            selection(L, l, k, seed);
      else if (k <= l + e)
        filtered[block_x + block_y * filtered_stride] = pivot;
      else
        filtered[block_x + block_y * filtered_stride] = selection(G, g, k -
                                                                  l - e, seed);
  }
}

static inline void pointwise_matrix_sub(double *mat1, double *mat2,
                                        double *diff, int width, int height,
                                        int stride1, int stride2,
                                        int diffstride) {
  int i, j;
  for (j = 0; j < height; ++j)
    for (i = 0; i < width; ++i) {
      diff[i + j * diffstride] = mat1[i + j * stride1] - mat2[i + j * stride2];
    }
}

static inline void pointwise_matrix_add(double *mat1, double *mat2,
                                        double *sum, int width, int height,
                                        int stride1, int stride2,
                                        int sumstride) {
  int i, j;
  for (j = 0; j < height; ++j)
    for (i = 0; i < width; ++i) {
      sum[i + j * sumstride] = mat1[i + j * stride1] + mat2[i + j * stride2];
    }
}

// Solves lucas kanade equation at any give pixel in the image using a local
// neighborhood for support. loc_x and loc_y are the x and y components
// of the center of the local neighborhood. Assumes dx, dy and dt have same
// stride. window is a wind_size x wind_size set of weights for the
// window_size x window_size neighborhood surrounding loc_x and loc_y.
static void optical_flow_per_pixel(double *dx, double *dy, double *dt,
                                   const double *window, int wind_size,
                                   flowMV *flow, int width, int height,
                                   int stride, int loc_x, int loc_y) {
  int i, j, iw, jw, im_ind;
  double g = 0;
  // M and b are matrices used to solve the equation a = M^-1 * b where a
  // are the desired optical flow parameters
  double M[4] = {0, 0, 0, 0};
  double b[2] = {0, 0};
  double det = 0;
  double trace2 = 0;
  double corner_score = 0;
  int step = wind_size >> 1;
  double gdx = 0, gdy = 0;

  for (j = loc_y - step; j <= loc_y + step; ++j)
    for (i = loc_x - step; i <= loc_x + step; ++i) {
      // if pixel is out of bounds, use reflected image content
      iw = WRAP_PIXEL(i, width);
      jw = WRAP_PIXEL(j, height);
      im_ind = iw + jw * stride;
      g = window[(i - loc_x + step) + (j - loc_y + step) * wind_size];
      gdx = g * dx[im_ind];
      gdy = g * dy[im_ind];
      M[0] += gdx * dx[im_ind];
      M[1] += gdx * dy[im_ind];
      M[3] += gdy * dy[im_ind];
      b[0] += -gdx * dt[im_ind];
      b[1] += -gdy * dt[im_ind];
    }
  M[2] = M[1];
  det = (M[0] * M[3]) - (M[1] * M[1]);
  trace2 = (M[0] + M[3]) * (M[0] + M[3]);
  if (RELATIVE_ERROR(det, 0) > MAX_ERROR) {
    const double det_inv = 1 / det;
    const double mult_b0 = det_inv * b[0];
    const double mult_b1 = det_inv * b[1];
    flow->u[loc_x + loc_y * flow->stride] = M[3] * mult_b0 - M[1] * mult_b1;
    flow->v[loc_x + loc_y * flow->stride] = -M[2] * mult_b0 + M[0] * mult_b1;
  } else {
    if (M[0] == 0 && M[3] == 0) {
      flow->u[loc_x + loc_y * flow->stride] = 0;
      flow->v[loc_x + loc_y * flow->stride] = 0;
    } else {
      const double trace2_inv = 1 / trace2;
      const double mult_b0 = trace2_inv * b[0];
      const double mult_b1 = trace2_inv * b[1];
      flow->u[loc_x + loc_y * flow->stride] = M[0] * mult_b0 + M[1] * mult_b1;
      flow->v[loc_x + loc_y * flow->stride] = M[2] * mult_b0 + M[3] * mult_b1;
    }
  }
  // compute inverse harris corner score as confidence metric
  corner_score = 1 / (det - (0.01 * trace2));
  flow->confidence[loc_x + loc_y * flow->stride] = corner_score > 2 ?
                                                   2 : corner_score;
}

// computes lucas kanade optical flow. Note that this assumes images that
// already denoised and differentiated.
static void lucas_kanade_base(double *smooth_frm, double *smooth_ref,
                              double *dx, double *dy, double *dt, flowMV *flow,
                              int width, int height, int smooth_frm_stride,
                              int smooth_ref_stride, int deriv_stride) {
  int i, j;
  const int local_neighborhood_sz = 5;
  // windowing function for local neighborhood weights
  const double window[] = {0.0039, 0.0156, 0.0234, 0.0156, 0.0039,
                           0.0156, 0.0625, 0.0938, 0.0625, 0.0156,
                           0.0234, 0.0938, 0.1406, 0.0938, 0.0234,
                           0.0156, 0.0625, 0.0938, 0.0625, 0.0156,
                           0.0039, 0.0156, 0.0234, 0.0156, 0.0039};
  // compute temporal derivative
  pointwise_matrix_sub(smooth_ref, smooth_frm, dt, width, height,
                       smooth_ref_stride, smooth_frm_stride, deriv_stride);

  for (j = 0; j < height; ++j)
    for (i = 0; i < width; ++i) {
      optical_flow_per_pixel(dx, dy, dt, window, local_neighborhood_sz,
                             flow, width, height, deriv_stride, i, j);
    }
}

// Improves an initial approximation for the vector field by iteratively
// warping one image towards the other.
static void iterative_refinement(unsigned char *ref, double *smooth_frm,
                                 double *dx, double *dy, double *dt,
                                 flowMV *flow, int width, int height,
                                 int ref_stride, int smooth_frm_stride,
                                 int deriv_stride, int n_refinements) {
  int i, j, n;
  double x, y;
  unsigned char *estimate = (unsigned char*)safe_malloc(width * height *
                                                        sizeof(*estimate));
  double *smooth_estimate = (double*)safe_malloc(width * height *
                                                 sizeof(*smooth_estimate));
  flowMV new_flow;
  new_flow.u = (double*)safe_malloc(width * height * sizeof(*new_flow.u));
  new_flow.v = (double*)safe_malloc(width * height * sizeof(*new_flow.v));
  new_flow.confidence = (double*)safe_malloc(width * height *
                                             sizeof(*new_flow.confidence));
  new_flow.stride = width;
  // warp one image toward the other
  for (n = 0; n < n_refinements; ++n) {
    for (j = 0; j < height; ++j)
      for (i = 0; i < width; ++i) {
        x = i - flow->u[i + j * flow->stride];
        y = j - flow->v[i + j * flow->stride];
        estimate[i + j * width] = interpolate(ref, x, y, width,
                                              height, ref_stride);
      }

    // compute flow between frame and warped estimate
    blur_img(estimate, width, height, width, smooth_estimate, width, BLUR_SIZE,
             BLUR_SIGMA, BLUR_AMP);
    lucas_kanade_base(smooth_frm, smooth_estimate, dx, dy, dt, &new_flow,
                      width, height, smooth_frm_stride, width, deriv_stride);

    // add residual and apply median filter for denoising
    pointwise_matrix_add(flow->u, new_flow.u, new_flow.u, width, height,
                         flow->stride, new_flow.stride, new_flow.stride);
    pointwise_matrix_add(flow->v, new_flow.v, new_flow.v, width, height,
                         flow->stride, new_flow.stride, new_flow.stride);
    median_filter(new_flow.u, flow->u, MEDIAN_FILTER_SIZE, width, height,
                  new_flow.stride, flow->stride);
    median_filter(new_flow.v, flow->v, MEDIAN_FILTER_SIZE, width, height,
                  new_flow.stride, flow->stride);
    median_filter(new_flow.confidence, flow->confidence, MEDIAN_FILTER_SIZE,
                  width, height, new_flow.stride, flow->stride);
  }

  free(smooth_estimate);
  free(estimate);
  free(new_flow.u);
  free(new_flow.v);
  free(new_flow.confidence);
}

// interface for computing optical flow.
void compute_flow(unsigned char *frm, unsigned char *ref,
                  double *u, double *v, double *confidence,
                  int width, int height, int frm_stride,
                  int ref_stride, int n_levels) {
  double *smooth_ref = (double*)safe_malloc(width * height *
                                            sizeof(*smooth_ref));
  double *smooth_frm = (double*)safe_malloc(width * height *
                                            sizeof(*smooth_frm));
  double *dx_frm = (double*)safe_malloc(width * height *
                                        sizeof(*dx_frm));
  double *dy_frm = (double*)safe_malloc(width * height *
                                       sizeof(*dx_frm));
  double *dt = (double *)safe_malloc(width * height * sizeof(*dt));
  flowMV flow;
  flow.u = u;
  flow.v = v;
  flow.confidence = confidence;
  flow.stride = width;
  if (n_levels == 1) {
    // Uses lucas kanade and iterative estimation without coarse to fine.
    // This is more efficient for small motion but becomes less accurate
    // as motion gets larger.
    blur_img(ref, width, height, ref_stride, smooth_ref, width, BLUR_SIZE,
             BLUR_SIGMA, BLUR_AMP);
    blur_img(frm, width, height, frm_stride, smooth_frm, width, BLUR_SIZE,
             BLUR_SIGMA, BLUR_AMP);
    differentiate(smooth_frm, width, height, width, dx_frm, dy_frm, width);
    lucas_kanade_base(smooth_frm, smooth_ref, dx_frm, dy_frm, dt, &flow,
                      width, height, width, width, width);
    iterative_refinement(ref, smooth_frm, dx_frm, dy_frm, dt, &flow, width,
                         height, ref_stride, width, width, ITERATIONS);
  } else {
    int i, j, k;
    // w, h, s_r, s_f are intermediate width, height, ref stride, frame stride
    // as the resolution changes up the pyramid
    int w, h, s_r, s_f, w_upscale, h_upscale, s_upscale;
    double x, y;
    double *smooth_frm_approx = (double*)safe_malloc(width * height *
                                             sizeof(*smooth_frm_approx));
    double *dx_frm_approx = (double*)safe_malloc(width * height *
                                             sizeof(*dx_frm_approx));
    double *dy_frm_approx = (double*)safe_malloc(width * height *
                                             sizeof(*dy_frm_approx));
    double *u_upscale = (double*)safe_malloc(width * height *
                                             sizeof(*u_upscale));
    double *v_upscale = (double*)safe_malloc(width * height *
                                             sizeof(*v_upscale));
    unsigned char *frm_approx = (unsigned char*)safe_malloc(width * height *
                                                           sizeof(*frm_approx));
    imPyramid *frm_pyramid = (imPyramid*)safe_malloc(sizeof(imPyramid));
    imPyramid *ref_pyramid = (imPyramid*)safe_malloc(sizeof(imPyramid));
    // create image pyramids
    image_pyramid(frm, frm_pyramid, width, height, frm_stride, n_levels);
    image_pyramid(ref, ref_pyramid, width, height, ref_stride, n_levels);
    n_levels = frm_pyramid->n_levels;
    w = frm_pyramid->widths[n_levels - 1];
    h = frm_pyramid->heights[n_levels - 1];
    s_r = ref_pyramid->strides[n_levels - 1];
    s_f = frm_pyramid->strides[n_levels - 1];

    // for each level in the pyramid starting with the coarsest
    for (i = n_levels - 1; i >= 0; --i) {
      assert(frm_pyramid->widths[i] == ref_pyramid->widths[i]);
      assert(frm_pyramid->heights[i] == ref_pyramid->heights[i]);
      blur_img(ref_pyramid->levels[i], w, h, s_r, smooth_ref, width, BLUR_SIZE,
               BLUR_SIGMA, BLUR_AMP);
      blur_img(frm_pyramid->levels[i], w, h, s_f, smooth_frm, width, BLUR_SIZE,
               BLUR_SIGMA, BLUR_AMP);
      differentiate(smooth_frm, w, h, width, dx_frm, dy_frm, width);


      //  Compute optical flow at this level between the reference frame and
      //  the estimate produced from warping. If this is the first iteration
      //  (meaning no warping has happened yet) then we have no approximation
      //  for the frame and only have to worry about flow between the original
      //  frame and reference. Every subsequent iteration requires computing
      //  and estimate for the frame based on previously computed flow vectors.

      if (i < n_levels - 1) {
        blur_img(frm_approx, w, h, width, smooth_frm_approx, width, BLUR_SIZE,
                 BLUR_SIGMA, BLUR_AMP);
        differentiate(smooth_frm_approx, w, h, width, dx_frm_approx,
                      dy_frm_approx, width);
        lucas_kanade_base(smooth_frm_approx, smooth_ref, dx_frm_approx,
                          dy_frm_approx, dt, &flow, w, h, width, width, width);
        pointwise_matrix_add(flow.u, u_upscale, u_upscale, w, h, flow.stride,
                             width, width);
        pointwise_matrix_add(flow.v, v_upscale, v_upscale, w, h, flow.stride,
                             width, width);
        median_filter(u_upscale, flow.u, MEDIAN_FILTER_SIZE, w, h, width,
                      flow.stride);
        median_filter(v_upscale, flow.v, MEDIAN_FILTER_SIZE, w, h, width,
                      flow.stride);
      } else {
        lucas_kanade_base(smooth_frm, smooth_ref, dx_frm,
                          dy_frm, dt, &flow, w, h, width, width, width);
      }
      iterative_refinement(ref_pyramid->levels[i], smooth_frm, dx_frm, dy_frm,
                           dt, &flow, w, h, s_r, width, width, ITERATIONS);


      // if we're at the finest level, we're ready to return u and v
      if (i == 0) {
        assert(w == width);
        assert(h == height);
        destruct_pyramid(frm_pyramid);
        destruct_pyramid(ref_pyramid);
        free(frm_pyramid);
        free(ref_pyramid);
        free(frm_approx);
        free(smooth_frm_approx);
        free(dx_frm_approx);
        free(dy_frm_approx);
        free(u_upscale);
        free(v_upscale);
        break;
      }

      w_upscale = ref_pyramid->widths[i - 1];
      h_upscale = ref_pyramid->heights[i - 1];
      s_upscale = ref_pyramid->strides[i - 1];

      // warp image according to optical flow estimate
      for (j = 0; j < h_upscale; ++j)
        for (k = 0; k < w_upscale; ++k) {
          u_upscale[k + j * w_upscale] = flow.u[(int)(k >> 1) +
                                                 (int)(j >> 1) * flow.stride];
          v_upscale[k + j * w_upscale] = flow.v[(int)(k >> 1) +
                                                 (int)(j >> 1) * flow.stride];
          x = k - u_upscale[k + j * w_upscale];
          y = j - v_upscale[k + j * w_upscale];
          frm_approx[k + j * width] = interpolate(ref_pyramid->levels[i - 1],
                                                  x, y, w_upscale, h_upscale,
                                                  s_upscale);
        }

      // assign dimensions for next level
      w = frm_pyramid->widths[i - 1];
      h = frm_pyramid->heights[i - 1];
      s_r = ref_pyramid->strides[i - 1];
      s_f = frm_pyramid->strides[i - 1];
    }
  }
  free(smooth_ref);
  free(smooth_frm);
  free(dx_frm);
  free(dy_frm);
  free(dt);
}
