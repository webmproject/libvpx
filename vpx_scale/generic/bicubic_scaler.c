/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <float.h>
#include <math.h>
#include <stdio.h>
#include "vpx_mem/vpx_mem.h"
#include "vpxscale_arbitrary.h"

#define FIXED_POINT

#define MAX_IN_WIDTH        800
#define MAX_IN_HEIGHT       600
#define MAX_OUT_WIDTH       800
#define MAX_OUT_HEIGHT      600
#define MAX_OUT_DIMENSION   ((MAX_OUT_WIDTH > MAX_OUT_HEIGHT) ? \
                             MAX_OUT_WIDTH : MAX_OUT_HEIGHT)

BICUBIC_SCALER_STRUCT g_b_scaler;
static int g_first_time = 1;

#pragma DATA_SECTION(g_hbuf, "VP6_HEAP")
#pragma DATA_ALIGN (g_hbuf, 32);
unsigned char g_hbuf[MAX_OUT_DIMENSION];

#pragma DATA_SECTION(g_hbuf_uv, "VP6_HEAP")
#pragma DATA_ALIGN (g_hbuf_uv, 32);
unsigned char g_hbuf_uv[MAX_OUT_DIMENSION];


#ifdef FIXED_POINT
static int a_i = 0.6 * 65536;
#else
static float a = -0.6;
#endif

#ifdef FIXED_POINT
//         3     2
// C0 = a*t - a*t
//
static short c0_fixed(unsigned int t) {
  // put t in Q16 notation
  unsigned short v1, v2;

  // Q16
  v1 = (a_i * t) >> 16;
  v1 = (v1 * t) >> 16;

  // Q16
  v2 = (a_i * t) >> 16;
  v2 = (v2 * t) >> 16;
  v2 = (v2 * t) >> 16;

  // Q12
  return -((v1 - v2) >> 4);
}

//                     2          3
// C1 = a*t + (3-2*a)*t  - (2-a)*t
//
static short c1_fixed(unsigned int t) {
  unsigned short v1, v2, v3;
  unsigned short two, three;

  // Q16
  v1 = (a_i * t) >> 16;

  // Q13
  two = 2 << 13;
  v2 = two - (a_i >> 3);
  v2 = (v2 * t) >> 16;
  v2 = (v2 * t) >> 16;
  v2 = (v2 * t) >> 16;

  // Q13
  three = 3 << 13;
  v3 = three - (2 * (a_i >> 3));
  v3 = (v3 * t) >> 16;
  v3 = (v3 * t) >> 16;

  // Q12
  return (((v1 >> 3) - v2 + v3) >> 1);

}

//                 2          3
// C2 = 1 - (3-a)*t  + (2-a)*t
//
static short c2_fixed(unsigned int t) {
  unsigned short v1, v2, v3;
  unsigned short two, three;

  // Q13
  v1 = 1 << 13;

  // Q13
  three = 3 << 13;
  v2 = three - (a_i >> 3);
  v2 = (v2 * t) >> 16;
  v2 = (v2 * t) >> 16;

  // Q13
  two = 2 << 13;
  v3 = two - (a_i >> 3);
  v3 = (v3 * t) >> 16;
  v3 = (v3 * t) >> 16;
  v3 = (v3 * t) >> 16;

  // Q12
  return (v1 - v2 + v3) >> 1;
}

//                 2      3
// C3 = a*t - 2*a*t  + a*t
//
static short c3_fixed(unsigned int t) {
  int v1, v2, v3;

  // Q16
  v1 = (a_i * t) >> 16;

  // Q15
  v2 = 2 * (a_i >> 1);
  v2 = (v2 * t) >> 16;
  v2 = (v2 * t) >> 16;

  // Q16
  v3 = (a_i * t) >> 16;
  v3 = (v3 * t) >> 16;
  v3 = (v3 * t) >> 16;

  // Q12
  return ((v2 - (v1 >> 1) - (v3 >> 1)) >> 3);
}
#else
//          3     2
// C0 = -a*t + a*t
//
float C0(float t) {
  return -a * t * t * t + a * t * t;
}

//                      2          3
// C1 = -a*t + (2*a+3)*t  - (a+2)*t
//
float C1(float t) {
  return -(a + 2.0f) * t * t * t + (2.0f * a + 3.0f) * t * t - a * t;
}

//                 2          3
// C2 = 1 - (a+3)*t  + (a+2)*t
//
float C2(float t) {
  return (a + 2.0f) * t * t * t - (a + 3.0f) * t * t + 1.0f;
}

//                 2      3
// C3 = a*t - 2*a*t  + a*t
//
float C3(float t) {
  return a * t * t * t - 2.0f * a * t * t + a * t;
}
#endif

#if 0
int compare_real_fixed() {
  int i, errors = 0;
  float mult = 1.0 / 10000.0;
  unsigned int fixed_mult = mult * 4294967296;// 65536;
  unsigned int phase_offset_int;
  float phase_offset_real;

  for (i = 0; i < 10000; i++) {
    int fixed0, fixed1, fixed2, fixed3, fixed_total;
    int real0, real1, real2, real3, real_total;

    phase_offset_real = (float)i * mult;
    phase_offset_int = (fixed_mult * i) >> 16;
//      phase_offset_int = phase_offset_real * 65536;

    fixed0 = c0_fixed(phase_offset_int);
    real0 = C0(phase_offset_real) * 4096.0;

    if ((abs(fixed0) > (abs(real0) + 1)) || (abs(fixed0) < (abs(real0) - 1)))
      errors++;

    fixed1 = c1_fixed(phase_offset_int);
    real1 = C1(phase_offset_real) * 4096.0;

    if ((abs(fixed1) > (abs(real1) + 1)) || (abs(fixed1) < (abs(real1) - 1)))
      errors++;

    fixed2 = c2_fixed(phase_offset_int);
    real2 = C2(phase_offset_real) * 4096.0;

    if ((abs(fixed2) > (abs(real2) + 1)) || (abs(fixed2) < (abs(real2) - 1)))
      errors++;

    fixed3 = c3_fixed(phase_offset_int);
    real3 = C3(phase_offset_real) * 4096.0;

    if ((abs(fixed3) > (abs(real3) + 1)) || (abs(fixed3) < (abs(real3) - 1)))
      errors++;

    fixed_total = fixed0 + fixed1 + fixed2 + fixed3;
    real_total = real0 + real1 + real2 + real3;

    if ((fixed_total > 4097) || (fixed_total < 4094))
      errors++;

    if ((real_total > 4097) || (real_total < 4095))
      errors++;
  }

  return errors;
}
#endif

// Find greatest common denominator between two integers.  Method used here is
//  slow compared to Euclid's algorithm, but does not require any division.
int gcd(int a, int b) {
  // Problem with this algorithm is that if a or b = 0 this function
  //  will never exit.  Don't want to return 0 because any computation
  //  that was based on a common denoninator and tried to reduce by
  //  dividing by 0 would fail.  Best solution that could be thought of
  //  would to be fail by returing a 1;
  if (a <= 0 || b <= 0)
    return 1;

  while (a != b) {
    if (b > a)
      b = b - a;
    else {
      int tmp = a;// swap large and
      a = b; // small
      b = tmp;
    }
  }

  return b;
}

void bicubic_coefficient_init() {
  vpx_memset(&g_b_scaler, 0, sizeof(BICUBIC_SCALER_STRUCT));
  g_first_time = 0;
}

void bicubic_coefficient_destroy() {
  if (!g_first_time) {
    vpx_free(g_b_scaler.l_w);

    vpx_free(g_b_scaler.l_h);

    vpx_free(g_b_scaler.l_h_uv);

    vpx_free(g_b_scaler.c_w);

    vpx_free(g_b_scaler.c_h);

    vpx_free(g_b_scaler.c_h_uv);

    vpx_memset(&g_b_scaler, 0, sizeof(BICUBIC_SCALER_STRUCT));
  }
}

// Create the coeffients that will be used for the cubic interpolation.
//  Because scaling does not have to be equal in the vertical and horizontal
//  regimes the phase offsets will be different.  There are 4 coefficents
//  for each point, two on each side.  The layout is that there are the
//  4 coefficents for each phase in the array and then the next phase.
int bicubic_coefficient_setup(int in_width, int in_height, int out_width, int out_height) {
  int i;
#ifdef FIXED_POINT
  int phase_offset_int;
  unsigned int fixed_mult;
  int product_val = 0;
#else
  float phase_offset;
#endif
  int gcd_w, gcd_h, gcd_h_uv, d_w, d_h, d_h_uv;

  if (g_first_time)
    bicubic_coefficient_init();


  // check to see if the coefficents have already been set up correctly
  if ((in_width == g_b_scaler.in_width) && (in_height == g_b_scaler.in_height)
      && (out_width == g_b_scaler.out_width) && (out_height == g_b_scaler.out_height))
    return 0;

  g_b_scaler.in_width = in_width;
  g_b_scaler.in_height = in_height;
  g_b_scaler.out_width = out_width;
  g_b_scaler.out_height = out_height;

  // Don't want to allow crazy scaling, just try and prevent a catastrophic
  //  failure here.  Want to fail after setting the member functions so if
  //  if the scaler is called the member functions will not scale.
  if (out_width <= 0 || out_height <= 0)
    return -1;

  // reduce in/out width and height ratios using the gcd
  gcd_w = gcd(out_width, in_width);
  gcd_h = gcd(out_height, in_height);
  gcd_h_uv = gcd(out_height, in_height / 2);

  // the numerator width and height are to be saved in
  //  globals so they can be used during the scaling process
  //  without having to be recalculated.
  g_b_scaler.nw = out_width / gcd_w;
  d_w = in_width / gcd_w;

  g_b_scaler.nh = out_height / gcd_h;
  d_h = in_height / gcd_h;

  g_b_scaler.nh_uv = out_height / gcd_h_uv;
  d_h_uv = (in_height / 2) / gcd_h_uv;

  // allocate memory for the coefficents
  vpx_free(g_b_scaler.l_w);

  vpx_free(g_b_scaler.l_h);

  vpx_free(g_b_scaler.l_h_uv);

  g_b_scaler.l_w = (short *)vpx_memalign(32, out_width * 2);
  g_b_scaler.l_h = (short *)vpx_memalign(32, out_height * 2);
  g_b_scaler.l_h_uv = (short *)vpx_memalign(32, out_height * 2);

  vpx_free(g_b_scaler.c_w);

  vpx_free(g_b_scaler.c_h);

  vpx_free(g_b_scaler.c_h_uv);

  g_b_scaler.c_w = (short *)vpx_memalign(32, g_b_scaler.nw * 4 * 2);
  g_b_scaler.c_h = (short *)vpx_memalign(32, g_b_scaler.nh * 4 * 2);
  g_b_scaler.c_h_uv = (short *)vpx_memalign(32, g_b_scaler.nh_uv * 4 * 2);

  g_b_scaler.hbuf = g_hbuf;
  g_b_scaler.hbuf_uv = g_hbuf_uv;

  // Set up polyphase filter taps.  This needs to be done before
  //  the scaling because of the floating point math required.  The
  //  coefficients are multiplied by 2^12 so that fixed point math
  //  can be used in the main scaling loop.
#ifdef FIXED_POINT
  fixed_mult = (1.0 / (float)g_b_scaler.nw) * 4294967296;

  product_val = 0;

  for (i = 0; i < g_b_scaler.nw; i++) {
    if (product_val > g_b_scaler.nw)
      product_val -= g_b_scaler.nw;

    phase_offset_int = (fixed_mult * product_val) >> 16;

    g_b_scaler.c_w[i * 4]   = c3_fixed(phase_offset_int);
    g_b_scaler.c_w[i * 4 + 1] = c2_fixed(phase_offset_int);
    g_b_scaler.c_w[i * 4 + 2] = c1_fixed(phase_offset_int);
    g_b_scaler.c_w[i * 4 + 3] = c0_fixed(phase_offset_int);

    product_val += d_w;
  }


  fixed_mult = (1.0 / (float)g_b_scaler.nh) * 4294967296;

  product_val = 0;

  for (i = 0; i < g_b_scaler.nh; i++) {
    if (product_val > g_b_scaler.nh)
      product_val -= g_b_scaler.nh;

    phase_offset_int = (fixed_mult * product_val) >> 16;

    g_b_scaler.c_h[i * 4]   = c0_fixed(phase_offset_int);
    g_b_scaler.c_h[i * 4 + 1] = c1_fixed(phase_offset_int);
    g_b_scaler.c_h[i * 4 + 2] = c2_fixed(phase_offset_int);
    g_b_scaler.c_h[i * 4 + 3] = c3_fixed(phase_offset_int);

    product_val += d_h;
  }

  fixed_mult = (1.0 / (float)g_b_scaler.nh_uv) * 4294967296;

  product_val = 0;

  for (i = 0; i < g_b_scaler.nh_uv; i++) {
    if (product_val > g_b_scaler.nh_uv)
      product_val -= g_b_scaler.nh_uv;

    phase_offset_int = (fixed_mult * product_val) >> 16;

    g_b_scaler.c_h_uv[i * 4]   = c0_fixed(phase_offset_int);
    g_b_scaler.c_h_uv[i * 4 + 1] = c1_fixed(phase_offset_int);
    g_b_scaler.c_h_uv[i * 4 + 2] = c2_fixed(phase_offset_int);
    g_b_scaler.c_h_uv[i * 4 + 3] = c3_fixed(phase_offset_int);

    product_val += d_h_uv;
  }

#else

  for (i = 0; i < g_nw; i++) {
    phase_offset = (float)((i * d_w) % g_nw) / (float)g_nw;
    g_c_w[i * 4]   = (C3(phase_offset) * 4096.0);
    g_c_w[i * 4 + 1] = (C2(phase_offset) * 4096.0);
    g_c_w[i * 4 + 2] = (C1(phase_offset) * 4096.0);
    g_c_w[i * 4 + 3] = (C0(phase_offset) * 4096.0);
  }

  for (i = 0; i < g_nh; i++) {
    phase_offset = (float)((i * d_h) % g_nh) / (float)g_nh;
    g_c_h[i * 4]   = (C0(phase_offset) * 4096.0);
    g_c_h[i * 4 + 1] = (C1(phase_offset) * 4096.0);
    g_c_h[i * 4 + 2] = (C2(phase_offset) * 4096.0);
    g_c_h[i * 4 + 3] = (C3(phase_offset) * 4096.0);
  }

  for (i = 0; i < g_nh_uv; i++) {
    phase_offset = (float)((i * d_h_uv) % g_nh_uv) / (float)g_nh_uv;
    g_c_h_uv[i * 4]   = (C0(phase_offset) * 4096.0);
    g_c_h_uv[i * 4 + 1] = (C1(phase_offset) * 4096.0);
    g_c_h_uv[i * 4 + 2] = (C2(phase_offset) * 4096.0);
    g_c_h_uv[i * 4 + 3] = (C3(phase_offset) * 4096.0);
  }

#endif

  // Create an array that corresponds input lines to output lines.
  //  This doesn't require floating point math, but it does require
  //  a division and because hardware division is not present that
  //  is a call.
  for (i = 0; i < out_width; i++) {
    g_b_scaler.l_w[i] = (i * d_w) / g_b_scaler.nw;

    if ((g_b_scaler.l_w[i] + 2) <= in_width)
      g_b_scaler.max_usable_out_width = i;

  }

  for (i = 0; i < out_height + 1; i++) {
    g_b_scaler.l_h[i] = (i * d_h) / g_b_scaler.nh;
    g_b_scaler.l_h_uv[i] = (i * d_h_uv) / g_b_scaler.nh_uv;
  }

  return 0;
}

int bicubic_scale(int in_width, int in_height, int in_stride,
                  int out_width, int out_height, int out_stride,
                  unsigned char *input_image, unsigned char *output_image) {
  short *RESTRICT l_w, * RESTRICT l_h;
  short *RESTRICT c_w, * RESTRICT c_h;
  unsigned char *RESTRICT ip, * RESTRICT op;
  unsigned char *RESTRICT hbuf;
  int h, w, lw, lh;
  int temp_sum;
  int phase_offset_w, phase_offset_h;

  c_w = g_b_scaler.c_w;
  c_h = g_b_scaler.c_h;

  op = output_image;

  l_w = g_b_scaler.l_w;
  l_h = g_b_scaler.l_h;

  phase_offset_h = 0;

  for (h = 0; h < out_height; h++) {
    // select the row to work on
    lh = l_h[h];
    ip = input_image + (in_stride * lh);

    // vp8_filter the row vertically into an temporary buffer.
    //  If the phase offset == 0 then all the multiplication
    //  is going to result in the output equalling the input.
    //  So instead point the temporary buffer to the input.
    //  Also handle the boundry condition of not being able to
    //  filter that last lines.
    if (phase_offset_h && (lh < in_height - 2)) {
      hbuf = g_b_scaler.hbuf;

      for (w = 0; w < in_width; w++) {
        temp_sum =  c_h[phase_offset_h * 4 + 3] * ip[w - in_stride];
        temp_sum += c_h[phase_offset_h * 4 + 2] * ip[w];
        temp_sum += c_h[phase_offset_h * 4 + 1] * ip[w + in_stride];
        temp_sum += c_h[phase_offset_h * 4]   * ip[w + 2 * in_stride];

        hbuf[w] = temp_sum >> 12;
      }
    } else
      hbuf = ip;

    // increase the phase offset for the next time around.
    if (++phase_offset_h >= g_b_scaler.nh)
      phase_offset_h = 0;

    // now filter and expand it horizontally into the final
    //  output buffer
    phase_offset_w = 0;

    for (w = 0; w < out_width; w++) {
      // get the index to use to expand the image
      lw = l_w[w];

      temp_sum =  c_w[phase_offset_w * 4]   * hbuf[lw - 1];
      temp_sum += c_w[phase_offset_w * 4 + 1] * hbuf[lw];
      temp_sum += c_w[phase_offset_w * 4 + 2] * hbuf[lw + 1];
      temp_sum += c_w[phase_offset_w * 4 + 3] * hbuf[lw + 2];
      temp_sum = temp_sum >> 12;

      if (++phase_offset_w >= g_b_scaler.nw)
        phase_offset_w = 0;

      // boundry conditions
      if ((lw + 2) >= in_width)
        temp_sum = hbuf[lw];

      if (lw == 0)
        temp_sum = hbuf[0];

      op[w] = temp_sum;
    }

    op += out_stride;
  }

  return 0;
}

void bicubic_scale_frame_reset() {
  g_b_scaler.out_width = 0;
  g_b_scaler.out_height = 0;
}

void bicubic_scale_frame(YV12_BUFFER_CONFIG *src, YV12_BUFFER_CONFIG *dst,
                         int new_width, int new_height) {

  dst->y_width = new_width;
  dst->y_height = new_height;
  dst->uv_width = new_width / 2;
  dst->uv_height = new_height / 2;

  dst->y_stride = dst->y_width;
  dst->uv_stride = dst->uv_width;

  bicubic_scale(src->y_width, src->y_height, src->y_stride,
                new_width, new_height, dst->y_stride,
                src->y_buffer, dst->y_buffer);

  bicubic_scale(src->uv_width, src->uv_height, src->uv_stride,
                new_width / 2, new_height / 2, dst->uv_stride,
                src->u_buffer, dst->u_buffer);

  bicubic_scale(src->uv_width, src->uv_height, src->uv_stride,
                new_width / 2, new_height / 2, dst->uv_stride,
                src->v_buffer, dst->v_buffer);
}
