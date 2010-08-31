/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_scale/yv12config.h"
#include "math.h"

#define C1 (float)(64 * 64 * 0.01*255*0.01*255)
#define C2 (float)(64 * 64 * 0.03*255*0.03*255)

static int width_y;
static int height_y;
static int height_uv;
static int width_uv;
static int stride_uv;
static int stride;
static int lumimask;
static int luminance;
static double plane_summed_weights = 0;

static short img12_sum_block[8*4096*4096*2] ;

static short img1_sum[8*4096*2];
static short img2_sum[8*4096*2];
static int   img1_sq_sum[8*4096*2];
static int   img2_sq_sum[8*4096*2];
static int   img12_mul_sum[8*4096*2];


double vp8_similarity
(
    int mu_x,
    int mu_y,
    int pre_mu_x2,
    int pre_mu_y2,
    int pre_mu_xy2
)
{
    int mu_x2, mu_y2, mu_xy, theta_x2, theta_y2, theta_xy;

    mu_x2 = mu_x * mu_x;
    mu_y2 = mu_y * mu_y;
    mu_xy = mu_x * mu_y;

    theta_x2 = 64 * pre_mu_x2 - mu_x2;
    theta_y2 = 64 * pre_mu_y2 - mu_y2;
    theta_xy = 64 * pre_mu_xy2 - mu_xy;

    return (2 * mu_xy + C1) * (2 * theta_xy + C2) / ((mu_x2 + mu_y2 + C1) * (theta_x2 + theta_y2 + C2));
}

double vp8_ssim
(
    const unsigned char *img1,
    const unsigned char *img2,
    int stride_img1,
    int stride_img2,
    int width,
    int height
)
{
    int x, y, x2, y2, img1_block, img2_block, img1_sq_block, img2_sq_block, img12_mul_block, temp;

    double plane_quality, weight, mean;

    short *img1_sum_ptr1, *img1_sum_ptr2;
    short *img2_sum_ptr1, *img2_sum_ptr2;
    int *img1_sq_sum_ptr1, *img1_sq_sum_ptr2;
    int *img2_sq_sum_ptr1, *img2_sq_sum_ptr2;
    int *img12_mul_sum_ptr1, *img12_mul_sum_ptr2;

    plane_quality = 0;

    if (lumimask)
        plane_summed_weights = 0.0f;
    else
        plane_summed_weights = (height - 7) * (width - 7);

    //some prologue for the main loop
    temp = 8 * width;

    img1_sum_ptr1      = img1_sum + temp;
    img2_sum_ptr1      = img2_sum + temp;
    img1_sq_sum_ptr1   = img1_sq_sum + temp;
    img2_sq_sum_ptr1   = img2_sq_sum + temp;
    img12_mul_sum_ptr1 = img12_mul_sum + temp;

    for (x = 0; x < width; x++)
    {
        img1_sum[x]      = img1[x];
        img2_sum[x]      = img2[x];
        img1_sq_sum[x]   = img1[x] * img1[x];
        img2_sq_sum[x]   = img2[x] * img2[x];
        img12_mul_sum[x] = img1[x] * img2[x];

        img1_sum_ptr1[x]      = 0;
        img2_sum_ptr1[x]      = 0;
        img1_sq_sum_ptr1[x]   = 0;
        img2_sq_sum_ptr1[x]   = 0;
        img12_mul_sum_ptr1[x] = 0;
    }

    //the main loop
    for (y = 1; y < height; y++)
    {
        img1 += stride_img1;
        img2 += stride_img2;

        temp = (y - 1) % 9 * width;

        img1_sum_ptr1      = img1_sum + temp;
        img2_sum_ptr1      = img2_sum + temp;
        img1_sq_sum_ptr1   = img1_sq_sum + temp;
        img2_sq_sum_ptr1   = img2_sq_sum + temp;
        img12_mul_sum_ptr1 = img12_mul_sum + temp;

        temp = y % 9 * width;

        img1_sum_ptr2      = img1_sum + temp;
        img2_sum_ptr2      = img2_sum + temp;
        img1_sq_sum_ptr2   = img1_sq_sum + temp;
        img2_sq_sum_ptr2   = img2_sq_sum + temp;
        img12_mul_sum_ptr2 = img12_mul_sum + temp;

        for (x = 0; x < width; x++)
        {
            img1_sum_ptr2[x]      = img1_sum_ptr1[x] + img1[x];
            img2_sum_ptr2[x]      = img2_sum_ptr1[x] + img2[x];
            img1_sq_sum_ptr2[x]   = img1_sq_sum_ptr1[x] + img1[x] * img1[x];
            img2_sq_sum_ptr2[x]   = img2_sq_sum_ptr1[x] + img2[x] * img2[x];
            img12_mul_sum_ptr2[x] = img12_mul_sum_ptr1[x] + img1[x] * img2[x];
        }

        if (y > 6)
        {
            //calculate the sum of the last 8 lines by subtracting the total sum of 8 lines back from the present sum
            temp = (y + 1) % 9 * width;

            img1_sum_ptr1      = img1_sum + temp;
            img2_sum_ptr1      = img2_sum + temp;
            img1_sq_sum_ptr1   = img1_sq_sum + temp;
            img2_sq_sum_ptr1   = img2_sq_sum + temp;
            img12_mul_sum_ptr1 = img12_mul_sum + temp;

            for (x = 0; x < width; x++)
            {
                img1_sum_ptr1[x]      = img1_sum_ptr2[x] - img1_sum_ptr1[x];
                img2_sum_ptr1[x]      = img2_sum_ptr2[x] - img2_sum_ptr1[x];
                img1_sq_sum_ptr1[x]   = img1_sq_sum_ptr2[x] - img1_sq_sum_ptr1[x];
                img2_sq_sum_ptr1[x]   = img2_sq_sum_ptr2[x] - img2_sq_sum_ptr1[x];
                img12_mul_sum_ptr1[x] = img12_mul_sum_ptr2[x] - img12_mul_sum_ptr1[x];
            }

            //here we calculate the sum over the 8x8 block of pixels
            //this is done by sliding a window across the column sums for the last 8 lines
            //each time adding the new column sum, and subtracting the one which fell out of the window
            img1_block      = 0;
            img2_block      = 0;
            img1_sq_block   = 0;
            img2_sq_block   = 0;
            img12_mul_block = 0;

            //prologue, and calculation of simularity measure from the first 8 column sums
            for (x = 0; x < 8; x++)
            {
                img1_block      += img1_sum_ptr1[x];
                img2_block      += img2_sum_ptr1[x];
                img1_sq_block   += img1_sq_sum_ptr1[x];
                img2_sq_block   += img2_sq_sum_ptr1[x];
                img12_mul_block += img12_mul_sum_ptr1[x];
            }

            if (lumimask)
            {
                y2 = y - 7;
                x2 = 0;

                if (luminance)
                {
                    mean = (img2_block + img1_block) / 128.0f;

                    if (!(y2 % 2 || x2 % 2))
                        *(img12_sum_block + y2 / 2 * width_uv + x2 / 2) = img2_block + img1_block;
                }
                else
                {
                    mean = *(img12_sum_block + y2 * width_uv + x2);
                    mean += *(img12_sum_block + y2 * width_uv + x2 + 4);
                    mean += *(img12_sum_block + (y2 + 4) * width_uv + x2);
                    mean += *(img12_sum_block + (y2 + 4) * width_uv + x2 + 4);

                    mean /= 512.0f;
                }

                weight = mean < 40 ? 0.0f :
                         (mean < 50 ? (mean - 40.0f) / 10.0f : 1.0f);
                plane_summed_weights += weight;

                plane_quality += weight * vp8_similarity(img1_block, img2_block, img1_sq_block, img2_sq_block, img12_mul_block);
            }
            else
                plane_quality += vp8_similarity(img1_block, img2_block, img1_sq_block, img2_sq_block, img12_mul_block);

            //and for the rest
            for (x = 8; x < width; x++)
            {
                img1_block      = img1_block + img1_sum_ptr1[x] - img1_sum_ptr1[x - 8];
                img2_block      = img2_block + img2_sum_ptr1[x] - img2_sum_ptr1[x - 8];
                img1_sq_block   = img1_sq_block + img1_sq_sum_ptr1[x] - img1_sq_sum_ptr1[x - 8];
                img2_sq_block   = img2_sq_block + img2_sq_sum_ptr1[x] - img2_sq_sum_ptr1[x - 8];
                img12_mul_block = img12_mul_block + img12_mul_sum_ptr1[x] - img12_mul_sum_ptr1[x - 8];

                if (lumimask)
                {
                    y2 = y - 7;
                    x2 = x - 7;

                    if (luminance)
                    {
                        mean = (img2_block + img1_block) / 128.0f;

                        if (!(y2 % 2 || x2 % 2))
                            *(img12_sum_block + y2 / 2 * width_uv + x2 / 2) = img2_block + img1_block;
                    }
                    else
                    {
                        mean = *(img12_sum_block + y2 * width_uv + x2);
                        mean += *(img12_sum_block + y2 * width_uv + x2 + 4);
                        mean += *(img12_sum_block + (y2 + 4) * width_uv + x2);
                        mean += *(img12_sum_block + (y2 + 4) * width_uv + x2 + 4);

                        mean /= 512.0f;
                    }

                    weight = mean < 40 ? 0.0f :
                             (mean < 50 ? (mean - 40.0f) / 10.0f : 1.0f);
                    plane_summed_weights += weight;

                    plane_quality += weight * vp8_similarity(img1_block, img2_block, img1_sq_block, img2_sq_block, img12_mul_block);
                }
                else
                    plane_quality += vp8_similarity(img1_block, img2_block, img1_sq_block, img2_sq_block, img12_mul_block);
            }
        }
    }

    if (plane_summed_weights == 0)
        return 1.0f;
    else
        return plane_quality / plane_summed_weights;
}

double vp8_calc_ssim
(
    YV12_BUFFER_CONFIG *source,
    YV12_BUFFER_CONFIG *dest,
    int lumamask,
    double *weight
)
{
    double a, b, c;
    double frame_weight;
    double ssimv;

    width_y = source->y_width;
    height_y = source->y_height;
    height_uv = source->uv_height;
    width_uv = source->uv_width;
    stride_uv = dest->uv_stride;
    stride = dest->y_stride;

    lumimask = lumamask;

    luminance = 1;
    a = vp8_ssim(source->y_buffer, dest->y_buffer,
                 source->y_stride, dest->y_stride, source->y_width, source->y_height);
    luminance = 0;

    frame_weight = plane_summed_weights / ((width_y - 7) * (height_y - 7));

    if (frame_weight == 0)
        a = b = c = 1.0f;
    else
    {
        b = vp8_ssim(source->u_buffer, dest->u_buffer,
                     source->uv_stride, dest->uv_stride, source->uv_width, source->uv_height);

        c = vp8_ssim(source->v_buffer, dest->v_buffer,
                     source->uv_stride, dest->uv_stride, source->uv_width, source->uv_height);
    }

    ssimv = a * .8 + .1 * (b + c);

    *weight = frame_weight;

    return ssimv;
}

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
