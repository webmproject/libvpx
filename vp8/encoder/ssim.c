/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_scale/yv12config.h"
#include "math.h"
#include "onyx_int.h"

#if CONFIG_RUNTIME_CPU_DETECT
#define IF_RTCD(x)  (x)
#else
#define IF_RTCD(x)  NULL
#endif
// Google version of SSIM
// SSIM
#define KERNEL 3
#define KERNEL_SIZE  (2 * KERNEL + 1)

typedef unsigned char uint8;
typedef unsigned int uint32;

static const int K[KERNEL_SIZE] =
{
    1, 4, 11, 16, 11, 4, 1    // 16 * exp(-0.3 * i * i)
};
static const double ki_w = 1. / 2304.;  // 1 / sum(i:0..6, j..6) K[i]*K[j]
double get_ssimg(const uint8 *org, const uint8 *rec,
                 int xo, int yo, int W, int H,
                 const int stride1, const int stride2
                )
{
    // TODO(skal): use summed tables
    int y, x;

    const int ymin = (yo - KERNEL < 0) ? 0 : yo - KERNEL;
    const int ymax = (yo + KERNEL > H - 1) ? H - 1 : yo + KERNEL;
    const int xmin = (xo - KERNEL < 0) ? 0 : xo - KERNEL;
    const int xmax = (xo + KERNEL > W - 1) ? W - 1 : xo + KERNEL;
    // worst case of accumulation is a weight of 48 = 16 + 2 * (11 + 4 + 1)
    // with a diff of 255, squares. That would a max error of 0x8ee0900,
    // which fits into 32 bits integers.
    uint32 w = 0, xm = 0, ym = 0, xxm = 0, xym = 0, yym = 0;
    org += ymin * stride1;
    rec += ymin * stride2;

    for (y = ymin; y <= ymax; ++y, org += stride1, rec += stride2)
    {
        const int Wy = K[KERNEL + y - yo];

        for (x = xmin; x <= xmax; ++x)
        {
            const  int Wxy = Wy * K[KERNEL + x - xo];
            // TODO(skal): inlined assembly
            w   += Wxy;
            xm  += Wxy * org[x];
            ym  += Wxy * rec[x];
            xxm += Wxy * org[x] * org[x];
            xym += Wxy * org[x] * rec[x];
            yym += Wxy * rec[x] * rec[x];
        }
    }

    {
        const double iw = 1. / w;
        const double iwx = xm * iw;
        const double iwy = ym * iw;
        double sxx = xxm * iw - iwx * iwx;
        double syy = yym * iw - iwy * iwy;

        // small errors are possible, due to rounding. Clamp to zero.
        if (sxx < 0.) sxx = 0.;

        if (syy < 0.) syy = 0.;

        {
            const double sxsy = sqrt(sxx * syy);
            const double sxy = xym * iw - iwx * iwy;
            static const double C11 = (0.01 * 0.01) * (255 * 255);
            static const double C22 = (0.03 * 0.03) * (255 * 255);
            static const double C33 = (0.015 * 0.015) * (255 * 255);
            const double l = (2. * iwx * iwy + C11) / (iwx * iwx + iwy * iwy + C11);
            const double c = (2. * sxsy      + C22) / (sxx + syy + C22);

            const double s = (sxy + C33) / (sxsy + C33);
            return l * c * s;

        }
    }

}

double get_ssimfull_kernelg(const uint8 *org, const uint8 *rec,
                            int xo, int yo, int W, int H,
                            const int stride1, const int stride2)
{
    // TODO(skal): use summed tables
    // worst case of accumulation is a weight of 48 = 16 + 2 * (11 + 4 + 1)
    // with a diff of 255, squares. That would a max error of 0x8ee0900,
    // which fits into 32 bits integers.
    int y_, x_;
    uint32 xm = 0, ym = 0, xxm = 0, xym = 0, yym = 0;
    org += (yo - KERNEL) * stride1;
    org += (xo - KERNEL);
    rec += (yo - KERNEL) * stride2;
    rec += (xo - KERNEL);

    for (y_ = 0; y_ < KERNEL_SIZE; ++y_, org += stride1, rec += stride2)
    {
        const int Wy = K[y_];

        for (x_ = 0; x_ < KERNEL_SIZE; ++x_)
        {
            const int Wxy = Wy * K[x_];
            // TODO(skal): inlined assembly
            const int org_x = org[x_];
            const int rec_x = rec[x_];
            xm  += Wxy * org_x;
            ym  += Wxy * rec_x;
            xxm += Wxy * org_x * org_x;
            xym += Wxy * org_x * rec_x;
            yym += Wxy * rec_x * rec_x;
        }
    }

    {
        const double iw = ki_w;
        const double iwx = xm * iw;
        const double iwy = ym * iw;
        double sxx = xxm * iw - iwx * iwx;
        double syy = yym * iw - iwy * iwy;

        // small errors are possible, due to rounding. Clamp to zero.
        if (sxx < 0.) sxx = 0.;

        if (syy < 0.) syy = 0.;

        {
            const double sxsy = sqrt(sxx * syy);
            const double sxy = xym * iw - iwx * iwy;
            static const double C11 = (0.01 * 0.01) * (255 * 255);
            static const double C22 = (0.03 * 0.03) * (255 * 255);
            static const double C33 = (0.015 * 0.015) * (255 * 255);
            const double l = (2. * iwx * iwy + C11) / (iwx * iwx + iwy * iwy + C11);
            const double c = (2. * sxsy      + C22) / (sxx + syy + C22);
            const double s = (sxy + C33) / (sxsy + C33);
            return l * c * s;
        }
    }
}

double calc_ssimg(const uint8 *org, const uint8 *rec,
                  const int image_width, const int image_height,
                  const int stride1, const int stride2
                 )
{
    int j, i;
    double SSIM = 0.;

    for (j = 0; j < KERNEL; ++j)
    {
        for (i = 0; i < image_width; ++i)
        {
            SSIM += get_ssimg(org, rec, i, j, image_width, image_height, stride1, stride2);
        }
    }

    for (j = KERNEL; j < image_height - KERNEL; ++j)
    {
        for (i = 0; i < KERNEL; ++i)
        {
            SSIM += get_ssimg(org, rec, i, j, image_width, image_height, stride1, stride2);
        }

        for (i = KERNEL; i < image_width - KERNEL; ++i)
        {
            SSIM += get_ssimfull_kernelg(org, rec, i, j,
                                         image_width, image_height, stride1, stride2);
        }

        for (i = image_width - KERNEL; i < image_width; ++i)
        {
            SSIM += get_ssimg(org, rec, i, j, image_width, image_height, stride1, stride2);
        }
    }

    for (j = image_height - KERNEL; j < image_height; ++j)
    {
        for (i = 0; i < image_width; ++i)
        {
            SSIM += get_ssimg(org, rec, i, j, image_width, image_height, stride1, stride2);
        }
    }

    return SSIM;
}


double vp8_calc_ssimg
(
    YV12_BUFFER_CONFIG *source,
    YV12_BUFFER_CONFIG *dest,
    double *ssim_y,
    double *ssim_u,
    double *ssim_v
)
{
    double ssim_all = 0;
    int ysize  = source->y_width * source->y_height;
    int uvsize = ysize / 4;

    *ssim_y = calc_ssimg(source->y_buffer, dest->y_buffer,
                         source->y_width, source->y_height,
                         source->y_stride, dest->y_stride);


    *ssim_u = calc_ssimg(source->u_buffer, dest->u_buffer,
                         source->uv_width, source->uv_height,
                         source->uv_stride, dest->uv_stride);


    *ssim_v = calc_ssimg(source->v_buffer, dest->v_buffer,
                         source->uv_width, source->uv_height,
                         source->uv_stride, dest->uv_stride);

    ssim_all = (*ssim_y + *ssim_u + *ssim_v) / (ysize + uvsize + uvsize);
    *ssim_y /= ysize;
    *ssim_u /= uvsize;
    *ssim_v /= uvsize;
    return ssim_all;
}


void ssim_parms_c
(
    unsigned char *s,
    int sp,
    unsigned char *r,
    int rp,
    unsigned long *sum_s,
    unsigned long *sum_r,
    unsigned long *sum_sq_s,
    unsigned long *sum_sq_r,
    unsigned long *sum_sxr
)
{
    int i,j;
    for(i=0;i<16;i++,s+=sp,r+=rp)
     {
         for(j=0;j<16;j++)
         {
             *sum_s += s[j];
             *sum_r += r[j];
             *sum_sq_s += s[j] * s[j];
             *sum_sq_r += r[j] * r[j];
             *sum_sxr += s[j] * r[j];
         }
     }
}
void ssim_parms_8x8_c
(
    unsigned char *s,
    int sp,
    unsigned char *r,
    int rp,
    unsigned long *sum_s,
    unsigned long *sum_r,
    unsigned long *sum_sq_s,
    unsigned long *sum_sq_r,
    unsigned long *sum_sxr
)
{
    int i,j;
    for(i=0;i<8;i++,s+=sp,r+=rp)
     {
         for(j=0;j<8;j++)
         {
             *sum_s += s[j];
             *sum_r += r[j];
             *sum_sq_s += s[j] * s[j];
             *sum_sq_r += r[j] * r[j];
             *sum_sxr += s[j] * r[j];
         }
     }
}

const static long long c1 =  426148; // (256^2*(.01*255)^2
const static long long c2 = 3835331; //(256^2*(.03*255)^2

static double similarity
(
    unsigned long sum_s,
    unsigned long sum_r,
    unsigned long sum_sq_s,
    unsigned long sum_sq_r,
    unsigned long sum_sxr,
    int count
)
{
    long long ssim_n = (2*sum_s*sum_r+ c1)*(2*count*sum_sxr-2*sum_s*sum_r+c2);

    long long ssim_d = (sum_s*sum_s +sum_r*sum_r+c1)*
            (count*sum_sq_s-sum_s*sum_s + count*sum_sq_r-sum_r*sum_r +c2) ;

    return ssim_n * 1.0 / ssim_d;
}

static double ssim_16x16(unsigned char *s,int sp, unsigned char *r,int rp,
            const vp8_variance_rtcd_vtable_t *rtcd)
{
    unsigned long sum_s=0,sum_r=0,sum_sq_s=0,sum_sq_r=0,sum_sxr=0;
    rtcd->ssimpf(s, sp, r, rp, &sum_s, &sum_r, &sum_sq_s, &sum_sq_r, &sum_sxr);
    return similarity(sum_s, sum_r, sum_sq_s, sum_sq_r, sum_sxr, 256);
}
static double ssim_8x8(unsigned char *s,int sp, unsigned char *r,int rp,
                const vp8_variance_rtcd_vtable_t *rtcd)
{
    unsigned long sum_s=0,sum_r=0,sum_sq_s=0,sum_sq_r=0,sum_sxr=0;
    rtcd->ssimpf_8x8(s, sp, r, rp, &sum_s, &sum_r, &sum_sq_s, &sum_sq_r, &sum_sxr);
    return similarity(sum_s, sum_r, sum_sq_s, sum_sq_r, sum_sxr, 64);
}

// TODO: (jbb) tried to scale this function such that we may be able to use it
// for distortion metric in mode selection code ( provided we do a reconstruction)
long dssim(unsigned char *s,int sp, unsigned char *r,int rp,
           const vp8_variance_rtcd_vtable_t *rtcd)
{
    unsigned long sum_s=0,sum_r=0,sum_sq_s=0,sum_sq_r=0,sum_sxr=0;
    double ssim3;
    long long ssim_n;
    long long ssim_d;

    rtcd->ssimpf(s, sp, r, rp, &sum_s, &sum_r, &sum_sq_s, &sum_sq_r, &sum_sxr);
    ssim_n = (2*sum_s*sum_r+ c1)*(2*256*sum_sxr-2*sum_s*sum_r+c2);

    ssim_d = (sum_s*sum_s +sum_r*sum_r+c1)*
            (256*sum_sq_s-sum_s*sum_s + 256*sum_sq_r-sum_r*sum_r +c2) ;

    ssim3 = 256 * (ssim_d-ssim_n) / ssim_d;
    return (long)( 256*ssim3 * ssim3 );
}
// TODO: (jbb) this 8x8 window might be too big + we may want to pick pixels
// such that the window regions overlap block boundaries to penalize blocking
// artifacts.

double vp8_ssim2
(
    unsigned char *img1,
    unsigned char *img2,
    int stride_img1,
    int stride_img2,
    int width,
    int height,
    const vp8_variance_rtcd_vtable_t *rtcd
)
{
    int i,j;

    double ssim_total=0;

    // we can sample points as frequently as we like start with 1 per 8x8
    for(i=0; i < height; i+=8, img1 += stride_img1*8, img2 += stride_img2*8)
    {
        for(j=0; j < width; j+=8 )
        {
            ssim_total += ssim_8x8(img1, stride_img1, img2, stride_img2, rtcd);
        }
    }
    ssim_total /= (width/8 * height /8);
    return ssim_total;

}
double vp8_calc_ssim
(
    YV12_BUFFER_CONFIG *source,
    YV12_BUFFER_CONFIG *dest,
    int lumamask,
    double *weight,
    const vp8_variance_rtcd_vtable_t *rtcd
)
{
    double a, b, c;
    double ssimv;

    a = vp8_ssim2(source->y_buffer, dest->y_buffer,
                 source->y_stride, dest->y_stride, source->y_width,
                 source->y_height, rtcd);

    b = vp8_ssim2(source->u_buffer, dest->u_buffer,
                 source->uv_stride, dest->uv_stride, source->uv_width,
                 source->uv_height, rtcd);

    c = vp8_ssim2(source->v_buffer, dest->v_buffer,
                 source->uv_stride, dest->uv_stride, source->uv_width,
                 source->uv_height, rtcd);

    ssimv = a * .8 + .1 * (b + c);

    *weight = 1;

    return ssimv;
}
