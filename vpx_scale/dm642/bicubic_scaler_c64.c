/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#include <float.h>
#include <math.h>
#include <stdio.h>
#include "vpx_mem/vpx_mem.h"
#include "vpxscale_arbitrary.h"

extern BICUBIC_SCALER_STRUCT g_b_scaler;

int bicubic_scale_c64(int in_width, int in_height, int in_stride,
                      int out_width, int out_height, int out_stride,
                      unsigned char *input_image, unsigned char *output_image)
{
    short *restrict l_w, * restrict l_h;
    short *restrict c_w, * restrict c_h;
    unsigned char *restrict ip, * restrict op, *restrict op_w;
    unsigned char *restrict hbuf;
    int h, w, lw, lh;
    int phase_offset_w, phase_offset_h;
    double coeff;
    int max_phase;

    c_w = g_b_scaler.c_w;
    c_h = g_b_scaler.c_h;

    op = output_image;

    l_w = g_b_scaler.l_w;
    l_h = g_b_scaler.l_h;

    phase_offset_h = 0;

    for (h = 0; h < out_height; h++)
    {
        // select the row to work on
        lh = l_h[h];
        ip = input_image + (in_stride * lh);

        coeff = _memd8_const(&c_h[phase_offset_h*4]);

        // vp8_filter the row vertically into an temporary buffer.
        //  If the phase offset == 0 then all the multiplication
        //  is going to result in the output equalling the input.
        //  So instead point the temporary buffer to the input.
        //  Also handle the boundry condition of not being able to
        //  filter that last lines.
        if (phase_offset_h && (lh < in_height - 2))
        {
            hbuf = g_b_scaler.hbuf;

            for (w = 0; w < in_width; w += 4)
            {
                int ip1, ip2, ip3, ip4;
                int y13_12, y11_10, y23_22, y21_20, y33_32, y31_30, y43_42, y41_40;
                int y10_20, y11_21, y12_22, y13_23, y30_40, y31_41, y32_42, y33_43;
                int s1, s2, s3, s4;

                ip1 = _mem4_const(&ip[w - in_stride]);
                ip2 = _mem4_const(&ip[w]);
                ip3 = _mem4_const(&ip[w + in_stride]);
                ip4 = _mem4_const(&ip[w + 2*in_stride]);

                // realignment of data.  Unpack the data so that it is in short
                //  format instead of bytes.
                y13_12 = _unpkhu4(ip1);
                y11_10 = _unpklu4(ip1);
                y23_22 = _unpkhu4(ip2);
                y21_20 = _unpklu4(ip2);
                y33_32 = _unpkhu4(ip3);
                y31_30 = _unpklu4(ip3);
                y43_42 = _unpkhu4(ip4);
                y41_40 = _unpklu4(ip4);

                // repack the data so that elements 1 and 2 are together.  this
                //  lines up so that a dot product with the coefficients can be
                //  done.
                y10_20 = _pack2(y11_10, y21_20);
                y11_21 = _packh2(y11_10, y21_20);
                y12_22 = _pack2(y13_12, y23_22);
                y13_23 = _packh2(y13_12, y23_22);

                s1 = _dotp2(_hi(coeff), y10_20);
                s2 = _dotp2(_hi(coeff), y11_21);
                s3 = _dotp2(_hi(coeff), y12_22);
                s4 = _dotp2(_hi(coeff), y13_23);

                y30_40 = _pack2(y31_30, y41_40);
                y31_41 = _packh2(y31_30, y41_40);
                y32_42 = _pack2(y33_32, y43_42);
                y33_43 = _packh2(y33_32, y43_42);

                // now repack elements 3 and 4 together.
                s1 += _dotp2(_lo(coeff), y30_40);
                s2 += _dotp2(_lo(coeff), y31_41);
                s3 += _dotp2(_lo(coeff), y32_42);
                s4 += _dotp2(_lo(coeff), y33_43);

                s1 = s1 >> 12;
                s2 = s2 >> 12;
                s3 = s3 >> 12;
                s4 = s4 >> 12;

                s1 = _pack2(s2, s1);
                s2 = _pack2(s4, s3);

                _amem4(&hbuf[w])  = _spacku4(s2, s1);
            }
        }
        else
            hbuf = ip;

        // increase the phase offset for the next time around.
        if (++phase_offset_h >= g_b_scaler.nh)
            phase_offset_h = 0;

        op_w = op;

        // will never be able to interpolate first pixel, so just copy it
        // over here.
        phase_offset_w = 1;
        *op_w++ = hbuf[0];

        if (1 >= g_b_scaler.nw) phase_offset_w = 0;

        max_phase = g_b_scaler.nw;

        for (w = 1; w < out_width; w++)
        {
            double coefficients;
            int hbuf_high, hbuf_low, hbuf_both;
            int sum_high, sum_low, sum;

            // get the index to use to expand the image
            lw = l_w[w];
            coefficients = _amemd8_const(&c_w[phase_offset_w*4]);
            hbuf_both = _mem4_const(&hbuf[lw-1]);

            hbuf_high = _unpkhu4(hbuf_both);
            hbuf_low  = _unpklu4(hbuf_both);

            sum_high = _dotp2(_hi(coefficients), hbuf_high);
            sum_low  = _dotp2(_lo(coefficients), hbuf_low);

            sum = (sum_high + sum_low) >> 12;

            if (++phase_offset_w >= max_phase)
                phase_offset_w = 0;

            if ((lw + 2) >= in_width)
                sum = hbuf[lw];

            *op_w++ = sum;
        }

        op += out_stride;
    }

    return 0;
}

void bicubic_scale_frame_c64(YV12_BUFFER_CONFIG *src, YV12_BUFFER_CONFIG *dst,
                             int new_width, int new_height)
{

    dst->y_width = new_width;
    dst->y_height = new_height;
    dst->uv_width = new_width / 2;
    dst->uv_height = new_height / 2;

    dst->y_stride = dst->y_width;
    dst->uv_stride = dst->uv_width;

    bicubic_scale_c64(src->y_width, src->y_height, src->y_stride,
                      new_width, new_height, dst->y_stride,
                      src->y_buffer, dst->y_buffer);

    bicubic_scale_c64(src->uv_width, src->uv_height, src->uv_stride,
                      new_width / 2, new_height / 2, dst->uv_stride,
                      src->u_buffer, dst->u_buffer);

    bicubic_scale_c64(src->uv_width, src->uv_height, src->uv_stride,
                      new_width / 2, new_height / 2, dst->uv_stride,
                      src->v_buffer, dst->v_buffer);
}
