/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "mcomp.h"
#include "vpx_mem/vpx_mem.h"

#include <stdio.h>
#include <limits.h>
#include <math.h>

#ifdef ENTROPY_STATS
static int mv_ref_ct [31] [4] [2];
static int mv_mode_cts [4] [2];
#endif

static int mv_bits_sadcost[256];

void vp8cx_init_mv_bits_sadcost()
{
    int i;

    for (i = 0; i < 256; i++)
    {
        mv_bits_sadcost[i] = (int)sqrt(i * 16);
    }
}


int vp8_mv_bit_cost(MV *mv, MV *ref, int *mvcost[2], int Weight)
{
    // MV costing is based on the distribution of vectors in the previous frame and as such will tend to
    // over state the cost of vectors. In addition coding a new vector can have a knock on effect on the
    // cost of subsequent vectors and the quality of prediction from NEAR and NEAREST for subsequent blocks.
    // The "Weight" parameter allows, to a limited extent, for some account to be taken of these factors.
    return ((mvcost[0][(mv->row - ref->row) >> 1] + mvcost[1][(mv->col - ref->col) >> 1]) * Weight) >> 7;
}

static int mv_err_cost(MV *mv, MV *ref, int *mvcost[2], int error_per_bit)
{
    //int i;
    //return ((mvcost[0][(mv->row - ref->row)>>1] + mvcost[1][(mv->col - ref->col)>>1] + 128) * error_per_bit) >> 8;
    //return ( (vp8_mv_bit_cost(mv,  ref, mvcost, 100) + 128) * error_per_bit) >> 8;

    //i = (vp8_mv_bit_cost(mv,  ref, mvcost, 100) * error_per_bit + 128) >> 8;
    return ((mvcost[0][(mv->row - ref->row) >> 1] + mvcost[1][(mv->col - ref->col) >> 1]) * error_per_bit + 128) >> 8;
    //return (vp8_mv_bit_cost(mv,  ref, mvcost, 128) * error_per_bit + 128) >> 8;
}


static int mv_bits(MV *mv, MV *ref, int *mvcost[2])
{
    // get the estimated number of bits for a motion vector, to be used for costing in SAD based
    // motion estimation
    return ((mvcost[0][(mv->row - ref->row) >> 1]  +  mvcost[1][(mv->col - ref->col)>> 1]) + 128) >> 8;
}

void vp8_init_dsmotion_compensation(MACROBLOCK *x, int stride)
{
    int Len;
    int search_site_count = 0;


    // Generate offsets for 4 search sites per step.
    Len = MAX_FIRST_STEP;
    x->ss[search_site_count].mv.col = 0;
    x->ss[search_site_count].mv.row = 0;
    x->ss[search_site_count].offset = 0;
    search_site_count++;

    while (Len > 0)
    {

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = 0;
        x->ss[search_site_count].mv.row = -Len;
        x->ss[search_site_count].offset = -Len * stride;
        search_site_count++;

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = 0;
        x->ss[search_site_count].mv.row = Len;
        x->ss[search_site_count].offset = Len * stride;
        search_site_count++;

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = -Len;
        x->ss[search_site_count].mv.row = 0;
        x->ss[search_site_count].offset = -Len;
        search_site_count++;

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = Len;
        x->ss[search_site_count].mv.row = 0;
        x->ss[search_site_count].offset = Len;
        search_site_count++;

        // Contract.
        Len /= 2;
    }

    x->ss_count = search_site_count;
    x->searches_per_step = 4;
}

void vp8_init3smotion_compensation(MACROBLOCK *x, int stride)
{
    int Len;
    int search_site_count = 0;

    // Generate offsets for 8 search sites per step.
    Len = MAX_FIRST_STEP;
    x->ss[search_site_count].mv.col = 0;
    x->ss[search_site_count].mv.row = 0;
    x->ss[search_site_count].offset = 0;
    search_site_count++;

    while (Len > 0)
    {

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = 0;
        x->ss[search_site_count].mv.row = -Len;
        x->ss[search_site_count].offset = -Len * stride;
        search_site_count++;

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = 0;
        x->ss[search_site_count].mv.row = Len;
        x->ss[search_site_count].offset = Len * stride;
        search_site_count++;

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = -Len;
        x->ss[search_site_count].mv.row = 0;
        x->ss[search_site_count].offset = -Len;
        search_site_count++;

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = Len;
        x->ss[search_site_count].mv.row = 0;
        x->ss[search_site_count].offset = Len;
        search_site_count++;

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = -Len;
        x->ss[search_site_count].mv.row = -Len;
        x->ss[search_site_count].offset = -Len * stride - Len;
        search_site_count++;

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = Len;
        x->ss[search_site_count].mv.row = -Len;
        x->ss[search_site_count].offset = -Len * stride + Len;
        search_site_count++;

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = -Len;
        x->ss[search_site_count].mv.row = Len;
        x->ss[search_site_count].offset = Len * stride - Len;
        search_site_count++;

        // Compute offsets for search sites.
        x->ss[search_site_count].mv.col = Len;
        x->ss[search_site_count].mv.row = Len;
        x->ss[search_site_count].offset = Len * stride + Len;
        search_site_count++;


        // Contract.
        Len /= 2;
    }

    x->ss_count = search_site_count;
    x->searches_per_step = 8;
}


#define MVC(r,c) (((mvcost[0][(r)-rr] + mvcost[1][(c) - rc]) * error_per_bit + 128 )>>8 ) // estimated cost of a motion vector (r,c)
#define PRE(r,c) (*(d->base_pre) + d->pre + ((r)>>2) * d->pre_stride + ((c)>>2)) // pointer to predictor base of a motionvector
#define SP(x) (((x)&3)<<1) // convert motion vector component to offset for svf calc
#define DIST(r,c) vfp->svf( PRE(r,c), d->pre_stride, SP(c),SP(r), z,b->src_stride,&sse) // returns subpixel variance error function.
#define IFMVCV(r,c,s,e) if ( c >= minc && c <= maxc && r >= minr && r <= maxr) s else e;
#define ERR(r,c) (MVC(r,c)+DIST(r,c)) // returns distortion + motion vector cost
#define CHECK_BETTER(v,r,c) IFMVCV(r,c,{if((v = ERR(r,c)) < besterr) { besterr = v; br=r; bc=c; }}, v=INT_MAX;)// checks if (r,c) has better score than previous best
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

//#define CHECK_BETTER(v,r,c) if((v = ERR(r,c)) < besterr) { besterr = v; br=r; bc=c; }

int vp8_find_best_sub_pixel_step_iteratively(MACROBLOCK *x, BLOCK *b, BLOCKD *d, MV *bestmv, MV *ref_mv, int error_per_bit, const vp8_variance_fn_ptr_t *vfp, int *mvcost[2])
{
    unsigned char *y = *(d->base_pre) + d->pre + (bestmv->row) * d->pre_stride + bestmv->col;
    unsigned char *z = (*(b->base_src) + b->src);

    int rr = ref_mv->row >> 1, rc = ref_mv->col >> 1;
    int br = bestmv->row << 2, bc = bestmv->col << 2;
    int tr = br, tc = bc;
    unsigned int besterr = INT_MAX;
    unsigned int left, right, up, down, diag;
    unsigned int sse;
    unsigned int whichdir;
    unsigned int halfiters = 4;
    unsigned int quarteriters = 4;

    int minc = MAX(x->mv_col_min << 2, (ref_mv->col >> 1) - ((1 << mvlong_width) - 1));
    int maxc = MIN(x->mv_col_max << 2, (ref_mv->col >> 1) + ((1 << mvlong_width) - 1));
    int minr = MAX(x->mv_row_min << 2, (ref_mv->row >> 1) - ((1 << mvlong_width) - 1));
    int maxr = MIN(x->mv_row_max << 2, (ref_mv->row >> 1) + ((1 << mvlong_width) - 1));

    // central mv
    bestmv->row <<= 3;
    bestmv->col <<= 3;

    // calculate central point error
    besterr = vfp->vf(y, d->pre_stride, z, b->src_stride, &sse);
    besterr += mv_err_cost(bestmv, ref_mv, mvcost, error_per_bit);

    // TODO: Each subsequent iteration checks at least one point in common with the last iteration could be 2 ( if diag selected)
    while (--halfiters)
    {
        // 1/2 pel
        CHECK_BETTER(left, tr, tc - 2);
        CHECK_BETTER(right, tr, tc + 2);
        CHECK_BETTER(up, tr - 2, tc);
        CHECK_BETTER(down, tr + 2, tc);

        whichdir = (left < right ? 0 : 1) + (up < down ? 0 : 2);

        switch (whichdir)
        {
        case 0:
            CHECK_BETTER(diag, tr - 2, tc - 2);
            break;
        case 1:
            CHECK_BETTER(diag, tr - 2, tc + 2);
            break;
        case 2:
            CHECK_BETTER(diag, tr + 2, tc - 2);
            break;
        case 3:
            CHECK_BETTER(diag, tr + 2, tc + 2);
            break;
        }

        // no reason to check the same one again.
        if (tr == br && tc == bc)
            break;

        tr = br;
        tc = bc;
    }

    // TODO: Each subsequent iteration checks at least one point in common with the last iteration could be 2 ( if diag selected)
    // 1/4 pel
    while (--quarteriters)
    {
        CHECK_BETTER(left, tr, tc - 1);
        CHECK_BETTER(right, tr, tc + 1);
        CHECK_BETTER(up, tr - 1, tc);
        CHECK_BETTER(down, tr + 1, tc);

        whichdir = (left < right ? 0 : 1) + (up < down ? 0 : 2);

        switch (whichdir)
        {
        case 0:
            CHECK_BETTER(diag, tr - 1, tc - 1);
            break;
        case 1:
            CHECK_BETTER(diag, tr - 1, tc + 1);
            break;
        case 2:
            CHECK_BETTER(diag, tr + 1, tc - 1);
            break;
        case 3:
            CHECK_BETTER(diag, tr + 1, tc + 1);
            break;
        }

        // no reason to check the same one again.
        if (tr == br && tc == bc)
            break;

        tr = br;
        tc = bc;
    }

    bestmv->row = br << 1;
    bestmv->col = bc << 1;

    if ((abs(bestmv->col - ref_mv->col) > MAX_FULL_PEL_VAL) || (abs(bestmv->row - ref_mv->row) > MAX_FULL_PEL_VAL))
        return INT_MAX;

    return besterr;
}
#undef MVC
#undef PRE
#undef SP
#undef DIST
#undef ERR
#undef CHECK_BETTER
#undef MIN
#undef MAX
int vp8_find_best_sub_pixel_step(MACROBLOCK *x, BLOCK *b, BLOCKD *d, MV *bestmv, MV *ref_mv, int error_per_bit, const vp8_variance_fn_ptr_t *vfp, int *mvcost[2])
{
    int bestmse = INT_MAX;
    MV startmv;
    //MV this_mv;
    MV this_mv;
    unsigned char *y = *(d->base_pre) + d->pre + (bestmv->row) * d->pre_stride + bestmv->col;
    unsigned char *z = (*(b->base_src) + b->src);
    int left, right, up, down, diag;
    unsigned int sse;
    int whichdir ;


    // Trap uncodable vectors
    if ((abs((bestmv->col << 3) - ref_mv->col) > MAX_FULL_PEL_VAL) || (abs((bestmv->row << 3) - ref_mv->row) > MAX_FULL_PEL_VAL))
    {
        bestmv->row <<= 3;
        bestmv->col <<= 3;
        return INT_MAX;
    }

    // central mv
    bestmv->row <<= 3;
    bestmv->col <<= 3;
    startmv = *bestmv;

    // calculate central point error
    bestmse = vfp->vf(y, d->pre_stride, z, b->src_stride, &sse);
    bestmse += mv_err_cost(bestmv, ref_mv, mvcost, error_per_bit);

    // go left then right and check error
    this_mv.row = startmv.row;
    this_mv.col = ((startmv.col - 8) | 4);
    left = vfp->svf_halfpix_h(y - 1, d->pre_stride, z, b->src_stride, &sse);
    left += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (left < bestmse)
    {
        *bestmv = this_mv;
        bestmse = left;
    }

    this_mv.col += 8;
    right = vfp->svf_halfpix_h(y, d->pre_stride, z, b->src_stride, &sse);
    right += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (right < bestmse)
    {
        *bestmv = this_mv;
        bestmse = right;
    }

    // go up then down and check error
    this_mv.col = startmv.col;
    this_mv.row = ((startmv.row - 8) | 4);
    up = vfp->svf_halfpix_v(y - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
    up += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (up < bestmse)
    {
        *bestmv = this_mv;
        bestmse = up;
    }

    this_mv.row += 8;
    down = vfp->svf_halfpix_v(y, d->pre_stride, z, b->src_stride, &sse);
    down += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (down < bestmse)
    {
        *bestmv = this_mv;
        bestmse = down;
    }


    // now check 1 more diagonal
    whichdir = (left < right ? 0 : 1) + (up < down ? 0 : 2);
    //for(whichdir =0;whichdir<4;whichdir++)
    //{
    this_mv = startmv;

    switch (whichdir)
    {
    case 0:
        this_mv.col = (this_mv.col - 8) | 4;
        this_mv.row = (this_mv.row - 8) | 4;
        diag = vfp->svf_halfpix_hv(y - 1 - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
        break;
    case 1:
        this_mv.col += 4;
        this_mv.row = (this_mv.row - 8) | 4;
        diag = vfp->svf_halfpix_hv(y - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
        break;
    case 2:
        this_mv.col = (this_mv.col - 8) | 4;
        this_mv.row += 4;
        diag = vfp->svf_halfpix_hv(y - 1, d->pre_stride, z, b->src_stride, &sse);
        break;
    case 3:
    default:
        this_mv.col += 4;
        this_mv.row += 4;
        diag = vfp->svf_halfpix_hv(y, d->pre_stride, z, b->src_stride, &sse);
        break;
    }

    diag += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

//  }


    // time to check quarter pels.
    if (bestmv->row < startmv.row)
        y -= d->pre_stride;

    if (bestmv->col < startmv.col)
        y--;

    startmv = *bestmv;



    // go left then right and check error
    this_mv.row = startmv.row;

    if (startmv.col & 7)
    {
        this_mv.col = startmv.col - 2;
        left = vfp->svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
    }
    else
    {
        this_mv.col = (startmv.col - 8) | 6;
        left = vfp->svf(y - 1, d->pre_stride, 6, this_mv.row & 7, z, b->src_stride, &sse);
    }

    left += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (left < bestmse)
    {
        *bestmv = this_mv;
        bestmse = left;
    }

    this_mv.col += 4;
    right = vfp->svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
    right += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (right < bestmse)
    {
        *bestmv = this_mv;
        bestmse = right;
    }

    // go up then down and check error
    this_mv.col = startmv.col;

    if (startmv.row & 7)
    {
        this_mv.row = startmv.row - 2;
        up = vfp->svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
    }
    else
    {
        this_mv.row = (startmv.row - 8) | 6;
        up = vfp->svf(y - d->pre_stride, d->pre_stride, this_mv.col & 7, 6, z, b->src_stride, &sse);
    }

    up += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (up < bestmse)
    {
        *bestmv = this_mv;
        bestmse = up;
    }

    this_mv.row += 4;
    down = vfp->svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
    down += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (down < bestmse)
    {
        *bestmv = this_mv;
        bestmse = down;
    }


    // now check 1 more diagonal
    whichdir = (left < right ? 0 : 1) + (up < down ? 0 : 2);

//  for(whichdir=0;whichdir<4;whichdir++)
//  {
    this_mv = startmv;

    switch (whichdir)
    {
    case 0:

        if (startmv.row & 7)
        {
            this_mv.row -= 2;

            if (startmv.col & 7)
            {
                this_mv.col -= 2;
                diag = vfp->svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
            }
            else
            {
                this_mv.col = (startmv.col - 8) | 6;
                diag = vfp->svf(y - 1, d->pre_stride, 6, this_mv.row & 7, z, b->src_stride, &sse);;
            }
        }
        else
        {
            this_mv.row = (startmv.row - 8) | 6;

            if (startmv.col & 7)
            {
                this_mv.col -= 2;
                diag = vfp->svf(y - d->pre_stride, d->pre_stride, this_mv.col & 7, 6, z, b->src_stride, &sse);
            }
            else
            {
                this_mv.col = (startmv.col - 8) | 6;
                diag = vfp->svf(y - d->pre_stride - 1, d->pre_stride, 6, 6, z, b->src_stride, &sse);
            }
        }

        break;
    case 1:
        this_mv.col += 2;

        if (startmv.row & 7)
        {
            this_mv.row -= 2;
            diag = vfp->svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
        }
        else
        {
            this_mv.row = (startmv.row - 8) | 6;
            diag = vfp->svf(y - d->pre_stride, d->pre_stride, this_mv.col & 7, 6, z, b->src_stride, &sse);
        }

        break;
    case 2:
        this_mv.row += 2;

        if (startmv.col & 7)
        {
            this_mv.col -= 2;
            diag = vfp->svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
        }
        else
        {
            this_mv.col = (startmv.col - 8) | 6;
            diag = vfp->svf(y - 1, d->pre_stride, 6, this_mv.row & 7, z, b->src_stride, &sse);;
        }

        break;
    case 3:
        this_mv.col += 2;
        this_mv.row += 2;
        diag = vfp->svf(y, d->pre_stride,  this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
        break;
    }

    diag += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

//  }

    return bestmse;
}

int vp8_find_best_half_pixel_step(MACROBLOCK *mb, BLOCK *b, BLOCKD *d, MV *bestmv, MV *ref_mv, int error_per_bit, const vp8_variance_fn_ptr_t *vfp, int *mvcost[2])
{
    int bestmse = INT_MAX;
    MV startmv;
    //MV this_mv;
    MV this_mv;
    unsigned char *y = *(d->base_pre) + d->pre + (bestmv->row) * d->pre_stride + bestmv->col;
    unsigned char *z = (*(b->base_src) + b->src);
    int left, right, up, down, diag;
    unsigned int sse;

    // Trap uncodable vectors
    if ((abs((bestmv->col << 3) - ref_mv->col) > MAX_FULL_PEL_VAL) || (abs((bestmv->row << 3) - ref_mv->row) > MAX_FULL_PEL_VAL))
    {
        bestmv->row <<= 3;
        bestmv->col <<= 3;
        return INT_MAX;
    }

    // central mv
    bestmv->row <<= 3;
    bestmv->col <<= 3;
    startmv = *bestmv;

    // calculate central point error
    bestmse = vfp->vf(y, d->pre_stride, z, b->src_stride, &sse);
    bestmse += mv_err_cost(bestmv, ref_mv, mvcost, error_per_bit);

    // go left then right and check error
    this_mv.row = startmv.row;
    this_mv.col = ((startmv.col - 8) | 4);
    left = vfp->svf_halfpix_h(y - 1, d->pre_stride, z, b->src_stride, &sse);
    left += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (left < bestmse)
    {
        *bestmv = this_mv;
        bestmse = left;
    }

    this_mv.col += 8;
    right = vfp->svf_halfpix_h(y, d->pre_stride, z, b->src_stride, &sse);
    right += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (right < bestmse)
    {
        *bestmv = this_mv;
        bestmse = right;
    }

    // go up then down and check error
    this_mv.col = startmv.col;
    this_mv.row = ((startmv.row - 8) | 4);
    up = vfp->svf_halfpix_v(y - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
    up += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (up < bestmse)
    {
        *bestmv = this_mv;
        bestmse = up;
    }

    this_mv.row += 8;
    down = vfp->svf_halfpix_v(y, d->pre_stride, z, b->src_stride, &sse);
    down += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (down < bestmse)
    {
        *bestmv = this_mv;
        bestmse = down;
    }

    // somewhat strangely not doing all the diagonals for half pel is slower than doing them.
#if 0
    // now check 1 more diagonal -
    whichdir = (left < right ? 0 : 1) + (up < down ? 0 : 2);
    this_mv = startmv;

    switch (whichdir)
    {
    case 0:
        this_mv.col = (this_mv.col - 8) | 4;
        this_mv.row = (this_mv.row - 8) | 4;
        diag = vfp->svf(y - 1 - d->pre_stride, d->pre_stride, 4, 4, z, b->src_stride, &sse);
        break;
    case 1:
        this_mv.col += 4;
        this_mv.row = (this_mv.row - 8) | 4;
        diag = vfp->svf(y - d->pre_stride, d->pre_stride, 4, 4, z, b->src_stride, &sse);
        break;
    case 2:
        this_mv.col = (this_mv.col - 8) | 4;
        this_mv.row += 4;
        diag = vfp->svf(y - 1, d->pre_stride, 4, 4, z, b->src_stride, &sse);
        break;
    case 3:
        this_mv.col += 4;
        this_mv.row += 4;
        diag = vfp->svf(y, d->pre_stride, 4, 4, z, b->src_stride, &sse);
        break;
    }

    diag += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

#else
    this_mv.col = (this_mv.col - 8) | 4;
    this_mv.row = (this_mv.row - 8) | 4;
    diag = vfp->svf_halfpix_hv(y - 1 - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
    diag += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

    this_mv.col += 8;
    diag = vfp->svf_halfpix_hv(y - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
    diag += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

    this_mv.col = (this_mv.col - 8) | 4;
    this_mv.row = startmv.row + 4;
    diag = vfp->svf_halfpix_hv(y - 1, d->pre_stride, z, b->src_stride, &sse);
    diag += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

    this_mv.col += 8;
    diag = vfp->svf_halfpix_hv(y, d->pre_stride, z, b->src_stride, &sse);
    diag += mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

#endif
    return bestmse;
}


#define MVC(r,c) (((mvsadcost[0][((r)<<2)-rr] + mvsadcost[1][((c)<<2) - rc]) * error_per_bit + 128 )>>8 ) // estimated cost of a motion vector (r,c)
#define PRE(r,c) (*(d->base_pre) + d->pre + (r) * d->pre_stride + (c)) // pointer to predictor base of a motionvector
#define DIST(r,c,v) vfp->sdf( src,src_stride,PRE(r,c),d->pre_stride, v) // returns sad error score.
#define ERR(r,c,v) (MVC(r,c)+DIST(r,c,v)) // returns distortion + motion vector cost
#define CHECK_BETTER(v,r,c) if ((v = ERR(r,c,besterr)) < besterr) { besterr = v; br=r; bc=c; } // checks if (r,c) has better score than previous best
static const MV next_chkpts[6][3] =
{
    {{ -2, 0}, { -1, -2}, {1, -2}},
    {{ -1, -2}, {1, -2}, {2, 0}},
    {{1, -2}, {2, 0}, {1, 2}},
    {{2, 0}, {1, 2}, { -1, 2}},
    {{1, 2}, { -1, 2}, { -2, 0}},
    {{ -1, 2}, { -2, 0}, { -1, -2}}
};
int vp8_hex_search
(
    MACROBLOCK *x,
    BLOCK *b,
    BLOCKD *d,
    MV *ref_mv,
    MV *best_mv,
    int search_param,
    int error_per_bit,
    int *num00,
    const vp8_variance_fn_ptr_t *vfp,
    int *mvsadcost[2],
    int *mvcost[2],
    MV *center_mv
)
{
    MV hex[6] = { { -1, -2}, {1, -2}, {2, 0}, {1, 2}, { -1, 2}, { -2, 0} } ;
    MV neighbors[8] = { { -1, -1}, {0, -1}, {1, -1}, { -1, 0}, {1, 0}, { -1, 1}, {0, 1}, {1, 1} } ;
    int i, j;
    unsigned char *src = (*(b->base_src) + b->src);
    int src_stride = b->src_stride;
    int rr = center_mv->row, rc = center_mv->col;
    int br = ref_mv->row >> 3, bc = ref_mv->col >> 3, tr, tc;
    unsigned int besterr, thiserr = 0x7fffffff;
    int k = -1, tk;

    if (bc < x->mv_col_min) bc = x->mv_col_min;

    if (bc > x->mv_col_max) bc = x->mv_col_max;

    if (br < x->mv_row_min) br = x->mv_row_min;

    if (br > x->mv_row_max) br = x->mv_row_max;

    rr >>= 1;
    rc >>= 1;

    besterr = ERR(br, bc, thiserr);

    // hex search
    //j=0
    tr = br;
    tc = bc;

    for (i = 0; i < 6; i++)
    {
        int nr = tr + hex[i].row, nc = tc + hex[i].col;

        if (nc < x->mv_col_min) continue;

        if (nc > x->mv_col_max) continue;

        if (nr < x->mv_row_min) continue;

        if (nr > x->mv_row_max) continue;

        //CHECK_BETTER(thiserr,nr,nc);
        if ((thiserr = ERR(nr, nc, besterr)) < besterr)
        {
            besterr = thiserr;
            br = nr;
            bc = nc;
            k = i;
        }
    }

    if (tr == br && tc == bc)
        goto cal_neighbors;

    for (j = 1; j < 127; j++)
    {
        tr = br;
        tc = bc;
        tk = k;

        for (i = 0; i < 3; i++)
        {
            int nr = tr + next_chkpts[tk][i].row, nc = tc + next_chkpts[tk][i].col;

            if (nc < x->mv_col_min) continue;

            if (nc > x->mv_col_max) continue;

            if (nr < x->mv_row_min) continue;

            if (nr > x->mv_row_max) continue;

            //CHECK_BETTER(thiserr,nr,nc);
            if ((thiserr = ERR(nr, nc, besterr)) < besterr)
            {
                besterr = thiserr;
                br = nr;
                bc = nc; //k=(tk+5+i)%6;}
                k = tk + 5 + i;

                if (k >= 12) k -= 12;
                else if (k >= 6) k -= 6;
            }
        }

        if (tr == br && tc == bc)
            break;
    }

    // check 8 1 away neighbors
cal_neighbors:
    tr = br;
    tc = bc;

    for (i = 0; i < 8; i++)
    {
        int nr = tr + neighbors[i].row, nc = tc + neighbors[i].col;

        if (nc < x->mv_col_min) continue;

        if (nc > x->mv_col_max) continue;

        if (nr < x->mv_row_min) continue;

        if (nr > x->mv_row_max) continue;

        CHECK_BETTER(thiserr, nr, nc);
    }

    best_mv->row = br;
    best_mv->col = bc;

    return vfp->vf(src, src_stride, PRE(br, bc), d->pre_stride, &thiserr) + mv_err_cost(best_mv, center_mv, mvcost, error_per_bit) ;
}
#undef MVC
#undef PRE
#undef SP
#undef DIST
#undef ERR
#undef CHECK_BETTER


int vp8_diamond_search_sad
(
    MACROBLOCK *x,
    BLOCK *b,
    BLOCKD *d,
    MV *ref_mv,
    MV *best_mv,
    int search_param,
    int error_per_bit,
    int *num00,
    vp8_variance_fn_ptr_t *fn_ptr,
    int *mvsadcost[2],
    int *mvcost[2],
    MV *center_mv
)
{
    int i, j, step;

    unsigned char *what = (*(b->base_src) + b->src);
    int what_stride = b->src_stride;
    unsigned char *in_what;
    int in_what_stride = d->pre_stride;
    unsigned char *best_address;

    int tot_steps;
    MV this_mv;

    int bestsad = INT_MAX;
    int best_site = 0;
    int last_site = 0;

    int ref_row = ref_mv->row >> 3;
    int ref_col = ref_mv->col >> 3;
    int this_row_offset;
    int this_col_offset;
    search_site *ss;

    unsigned char *check_here;
    int thissad;

    *num00 = 0;

    // Work out the start point for the search
    in_what = (unsigned char *)(*(d->base_pre) + d->pre + (ref_row * (d->pre_stride)) + ref_col);
    best_address = in_what;

    // We need to check that the starting point for the search (as indicated by ref_mv) is within the buffer limits
    if ((ref_col > x->mv_col_min) && (ref_col < x->mv_col_max) &&
    (ref_row > x->mv_row_min) && (ref_row < x->mv_row_max))
    {
        // Check the starting position
        bestsad = fn_ptr->sdf(what, what_stride, in_what, in_what_stride, 0x7fffffff) + mv_err_cost(ref_mv, center_mv, mvsadcost, error_per_bit);
    }

    // search_param determines the length of the initial step and hence the number of iterations
    // 0 = initial step (MAX_FIRST_STEP) pel : 1 = (MAX_FIRST_STEP/2) pel, 2 = (MAX_FIRST_STEP/4) pel... etc.
    ss = &x->ss[search_param * x->searches_per_step];
    tot_steps = (x->ss_count / x->searches_per_step) - search_param;

    i = 1;
    best_mv->row = ref_row;
    best_mv->col = ref_col;

    for (step = 0; step < tot_steps ; step++)
    {
        for (j = 0 ; j < x->searches_per_step ; j++)
        {
            // Trap illegal vectors
            this_row_offset = best_mv->row + ss[i].mv.row;
            this_col_offset = best_mv->col + ss[i].mv.col;

            if ((this_col_offset > x->mv_col_min) && (this_col_offset < x->mv_col_max) &&
            (this_row_offset > x->mv_row_min) && (this_row_offset < x->mv_row_max))

            {
                check_here = ss[i].offset + best_address;
                thissad = fn_ptr->sdf(what, what_stride, check_here , in_what_stride, bestsad);

                if (thissad < bestsad)
                {
                    this_mv.row = this_row_offset << 3;
                    this_mv.col = this_col_offset << 3;
                    thissad += mv_err_cost(&this_mv, center_mv, mvsadcost, error_per_bit);

                    if (thissad < bestsad)
                    {
                        bestsad = thissad;
                        best_site = i;
                    }
                }
            }

            i++;
        }

        if (best_site != last_site)
        {
            best_mv->row += ss[best_site].mv.row;
            best_mv->col += ss[best_site].mv.col;
            best_address += ss[best_site].offset;
            last_site = best_site;
        }
        else if (best_address == in_what)
            (*num00)++;
    }

    this_mv.row = best_mv->row << 3;
    this_mv.col = best_mv->col << 3;

    if (bestsad == INT_MAX)
        return INT_MAX;

    return fn_ptr->vf(what, what_stride, best_address, in_what_stride, (unsigned int *)(&thissad))
    + mv_err_cost(&this_mv, center_mv, mvcost, error_per_bit);
}

int vp8_diamond_search_sadx4
(
    MACROBLOCK *x,
    BLOCK *b,
    BLOCKD *d,
    MV *ref_mv,
    MV *best_mv,
    int search_param,
    int error_per_bit,
    int *num00,
    vp8_variance_fn_ptr_t *fn_ptr,
    int *mvsadcost[2],
    int *mvcost[2],
    MV *center_mv
)
{
    int i, j, step;

    unsigned char *what = (*(b->base_src) + b->src);
    int what_stride = b->src_stride;
    unsigned char *in_what;
    int in_what_stride = d->pre_stride;
    unsigned char *best_address;

    int tot_steps;
    MV this_mv;

    int bestsad = INT_MAX;
    int best_site = 0;
    int last_site = 0;

    int ref_row = ref_mv->row >> 3;
    int ref_col = ref_mv->col >> 3;
    int this_row_offset;
    int this_col_offset;
    search_site *ss;

    unsigned char *check_here;
    unsigned int thissad;

    *num00 = 0;

    // Work out the start point for the search
    in_what = (unsigned char *)(*(d->base_pre) + d->pre + (ref_row * (d->pre_stride)) + ref_col);
    best_address = in_what;

    // We need to check that the starting point for the search (as indicated by ref_mv) is within the buffer limits
    if ((ref_col > x->mv_col_min) && (ref_col < x->mv_col_max) &&
    (ref_row > x->mv_row_min) && (ref_row < x->mv_row_max))
    {
        // Check the starting position
        bestsad = fn_ptr->sdf(what, what_stride, in_what, in_what_stride, 0x7fffffff) + mv_err_cost(ref_mv, center_mv, mvsadcost, error_per_bit);
    }

    // search_param determines the length of the initial step and hence the number of iterations
    // 0 = initial step (MAX_FIRST_STEP) pel : 1 = (MAX_FIRST_STEP/2) pel, 2 = (MAX_FIRST_STEP/4) pel... etc.
    ss = &x->ss[search_param * x->searches_per_step];
    tot_steps = (x->ss_count / x->searches_per_step) - search_param;

    i = 1;
    best_mv->row = ref_row;
    best_mv->col = ref_col;

    for (step = 0; step < tot_steps ; step++)
    {
        int all_in = 1, t;

        // To know if all neighbor points are within the bounds, 4 bounds checking are enough instead of
        // checking 4 bounds for each points.
        all_in &= ((best_mv->row + ss[i].mv.row)> x->mv_row_min);
        all_in &= ((best_mv->row + ss[i+1].mv.row) < x->mv_row_max);
        all_in &= ((best_mv->col + ss[i+2].mv.col) > x->mv_col_min);
        all_in &= ((best_mv->col + ss[i+3].mv.col) < x->mv_col_max);

        if (all_in)
        {
            unsigned int sad_array[4];

            for (j = 0 ; j < x->searches_per_step ; j += 4)
            {
                unsigned char *block_offset[4];

                for (t = 0; t < 4; t++)
                    block_offset[t] = ss[i+t].offset + best_address;

                fn_ptr->sdx4df(what, what_stride, block_offset, in_what_stride, sad_array);

                for (t = 0; t < 4; t++, i++)
                {
                    if (sad_array[t] < bestsad)
                    {
                        this_mv.row = (best_mv->row + ss[i].mv.row) << 3;
                        this_mv.col = (best_mv->col + ss[i].mv.col) << 3;
                        sad_array[t] += mv_err_cost(&this_mv, center_mv, mvsadcost, error_per_bit);

                        if (sad_array[t] < bestsad)
                        {
                            bestsad = sad_array[t];
                            best_site = i;
                        }
                    }
                }
            }
        }
        else
        {
            for (j = 0 ; j < x->searches_per_step ; j++)
            {
                // Trap illegal vectors
                this_row_offset = best_mv->row + ss[i].mv.row;
                this_col_offset = best_mv->col + ss[i].mv.col;

                if ((this_col_offset > x->mv_col_min) && (this_col_offset < x->mv_col_max) &&
                (this_row_offset > x->mv_row_min) && (this_row_offset < x->mv_row_max))
                {
                    check_here = ss[i].offset + best_address;
                    thissad = fn_ptr->sdf(what, what_stride, check_here , in_what_stride, bestsad);

                    if (thissad < bestsad)
                    {
                        this_mv.row = this_row_offset << 3;
                        this_mv.col = this_col_offset << 3;
                        thissad += mv_err_cost(&this_mv, center_mv, mvsadcost, error_per_bit);

                        if (thissad < bestsad)
                        {
                            bestsad = thissad;
                            best_site = i;
                        }
                    }
                }
                i++;
            }
        }

        if (best_site != last_site)
        {
            best_mv->row += ss[best_site].mv.row;
            best_mv->col += ss[best_site].mv.col;
            best_address += ss[best_site].offset;
            last_site = best_site;
        }
        else if (best_address == in_what)
            (*num00)++;
    }

    this_mv.row = best_mv->row << 3;
    this_mv.col = best_mv->col << 3;

    if (bestsad == INT_MAX)
        return INT_MAX;

    return fn_ptr->vf(what, what_stride, best_address, in_what_stride, (unsigned int *)(&thissad))
    + mv_err_cost(&this_mv, center_mv, mvcost, error_per_bit);
}


#if !(CONFIG_REALTIME_ONLY)
int vp8_full_search_sad(MACROBLOCK *x, BLOCK *b, BLOCKD *d, MV *ref_mv, int error_per_bit, int distance, vp8_variance_fn_ptr_t *fn_ptr, int *mvcost[2], int *mvsadcost[2], MV *center_mv)
{
    unsigned char *what = (*(b->base_src) + b->src);
    int what_stride = b->src_stride;
    unsigned char *in_what;
    int in_what_stride = d->pre_stride;
    int mv_stride = d->pre_stride;
    unsigned char *bestaddress;
    MV *best_mv = &d->bmi.mv.as_mv;
    MV this_mv;
    int bestsad = INT_MAX;
    int r, c;

    unsigned char *check_here;
    int thissad;

    int ref_row = ref_mv->row >> 3;
    int ref_col = ref_mv->col >> 3;

    int row_min = ref_row - distance;
    int row_max = ref_row + distance;
    int col_min = ref_col - distance;
    int col_max = ref_col + distance;

    // Work out the mid point for the search
    in_what = *(d->base_pre) + d->pre;
    bestaddress = in_what + (ref_row * d->pre_stride) + ref_col;

    best_mv->row = ref_row;
    best_mv->col = ref_col;

    // We need to check that the starting point for the search (as indicated by ref_mv) is within the buffer limits
    if ((ref_col > x->mv_col_min) && (ref_col < x->mv_col_max) &&
    (ref_row > x->mv_row_min) && (ref_row < x->mv_row_max))
    {
        // Baseline value at the centre

        //bestsad = fn_ptr->sf( what,what_stride,bestaddress,in_what_stride) + (int)sqrt(mv_err_cost(ref_mv,ref_mv, mvcost,error_per_bit*14));
        bestsad = fn_ptr->sdf(what, what_stride, bestaddress, in_what_stride, 0x7fffffff) + mv_err_cost(ref_mv, center_mv, mvsadcost, error_per_bit);
    }

    // Apply further limits to prevent us looking using vectors that stretch beyiond the UMV border
    if (col_min < x->mv_col_min)
        col_min = x->mv_col_min;

    if (col_max > x->mv_col_max)
        col_max = x->mv_col_max;

    if (row_min < x->mv_row_min)
        row_min = x->mv_row_min;

    if (row_max > x->mv_row_max)
        row_max = x->mv_row_max;

    for (r = row_min; r < row_max ; r++)
    {
        this_mv.row = r << 3;
        check_here = r * mv_stride + in_what + col_min;

        for (c = col_min; c < col_max; c++)
        {
            thissad = fn_ptr->sdf(what, what_stride, check_here , in_what_stride, bestsad);

            this_mv.col = c << 3;
            //thissad += (int)sqrt(mv_err_cost(&this_mv,ref_mv, mvcost,error_per_bit*14));
            //thissad  += error_per_bit * mv_bits_sadcost[mv_bits(&this_mv, ref_mv, mvcost)];
            thissad  += mv_err_cost(&this_mv, center_mv, mvsadcost, error_per_bit); //mv_bits(error_per_bit, &this_mv, ref_mv, mvsadcost);

            if (thissad < bestsad)
            {
                bestsad = thissad;
                best_mv->row = r;
                best_mv->col = c;
                bestaddress = check_here;
            }

            check_here++;
        }
    }

    this_mv.row = best_mv->row << 3;
    this_mv.col = best_mv->col << 3;

    if (bestsad < INT_MAX)
        return fn_ptr->vf(what, what_stride, bestaddress, in_what_stride, (unsigned int *)(&thissad))
        + mv_err_cost(&this_mv, center_mv, mvcost, error_per_bit);
    else
        return INT_MAX;
}

int vp8_full_search_sadx3(MACROBLOCK *x, BLOCK *b, BLOCKD *d, MV *ref_mv, int error_per_bit, int distance, vp8_variance_fn_ptr_t *fn_ptr, int *mvcost[2], int *mvsadcost[2], MV *center_mv)
{
    unsigned char *what = (*(b->base_src) + b->src);
    int what_stride = b->src_stride;
    unsigned char *in_what;
    int in_what_stride = d->pre_stride;
    int mv_stride = d->pre_stride;
    unsigned char *bestaddress;
    MV *best_mv = &d->bmi.mv.as_mv;
    MV this_mv;
    int bestsad = INT_MAX;
    int r, c;

    unsigned char *check_here;
    unsigned int thissad;

    int ref_row = ref_mv->row >> 3;
    int ref_col = ref_mv->col >> 3;

    int row_min = ref_row - distance;
    int row_max = ref_row + distance;
    int col_min = ref_col - distance;
    int col_max = ref_col + distance;

    unsigned int sad_array[3];

    // Work out the mid point for the search
    in_what = *(d->base_pre) + d->pre;
    bestaddress = in_what + (ref_row * d->pre_stride) + ref_col;

    best_mv->row = ref_row;
    best_mv->col = ref_col;

    // We need to check that the starting point for the search (as indicated by ref_mv) is within the buffer limits
    if ((ref_col > x->mv_col_min) && (ref_col < x->mv_col_max) &&
    (ref_row > x->mv_row_min) && (ref_row < x->mv_row_max))
    {
        // Baseline value at the centre
        bestsad = fn_ptr->sdf(what, what_stride, bestaddress, in_what_stride, 0x7fffffff) + mv_err_cost(ref_mv, center_mv, mvsadcost, error_per_bit);
    }

    // Apply further limits to prevent us looking using vectors that stretch beyiond the UMV border
    if (col_min < x->mv_col_min)
        col_min = x->mv_col_min;

    if (col_max > x->mv_col_max)
        col_max = x->mv_col_max;

    if (row_min < x->mv_row_min)
        row_min = x->mv_row_min;

    if (row_max > x->mv_row_max)
        row_max = x->mv_row_max;

    for (r = row_min; r < row_max ; r++)
    {
        this_mv.row = r << 3;
        check_here = r * mv_stride + in_what + col_min;
        c = col_min;

        while ((c + 2) < col_max)
        {
            int i;

            fn_ptr->sdx3f(what, what_stride, check_here , in_what_stride, sad_array);

            for (i = 0; i < 3; i++)
            {
                thissad = sad_array[i];

                if (thissad < bestsad)
                {
                    this_mv.col = c << 3;
                    thissad  += mv_err_cost(&this_mv, center_mv, mvsadcost, error_per_bit);

                    if (thissad < bestsad)
                    {
                        bestsad = thissad;
                        best_mv->row = r;
                        best_mv->col = c;
                        bestaddress = check_here;
                    }
                }

                check_here++;
                c++;
            }
        }

        while (c < col_max)
        {
            thissad = fn_ptr->sdf(what, what_stride, check_here , in_what_stride, bestsad);

            if (thissad < bestsad)
            {
                this_mv.col = c << 3;
                thissad  += mv_err_cost(&this_mv, center_mv, mvsadcost, error_per_bit);

                if (thissad < bestsad)
                {
                    bestsad = thissad;
                    best_mv->row = r;
                    best_mv->col = c;
                    bestaddress = check_here;
                }
            }

            check_here ++;
            c ++;
        }

    }

    this_mv.row = best_mv->row << 3;
    this_mv.col = best_mv->col << 3;

    if (bestsad < INT_MAX)
        return fn_ptr->vf(what, what_stride, bestaddress, in_what_stride, (unsigned int *)(&thissad))
        + mv_err_cost(&this_mv, center_mv, mvcost, error_per_bit);
    else
        return INT_MAX;
}

int vp8_full_search_sadx8(MACROBLOCK *x, BLOCK *b, BLOCKD *d, MV *ref_mv, int error_per_bit, int distance, vp8_variance_fn_ptr_t *fn_ptr, int *mvcost[2], int *mvsadcost[2], MV *center_mv)
{
    unsigned char *what = (*(b->base_src) + b->src);
    int what_stride = b->src_stride;
    unsigned char *in_what;
    int in_what_stride = d->pre_stride;
    int mv_stride = d->pre_stride;
    unsigned char *bestaddress;
    MV *best_mv = &d->bmi.mv.as_mv;
    MV this_mv;
    int bestsad = INT_MAX;
    int r, c;

    unsigned char *check_here;
    unsigned int thissad;

    int ref_row = ref_mv->row >> 3;
    int ref_col = ref_mv->col >> 3;

    int row_min = ref_row - distance;
    int row_max = ref_row + distance;
    int col_min = ref_col - distance;
    int col_max = ref_col + distance;

    DECLARE_ALIGNED_ARRAY(16, unsigned short, sad_array8, 8);
    unsigned int sad_array[3];

    // Work out the mid point for the search
    in_what = *(d->base_pre) + d->pre;
    bestaddress = in_what + (ref_row * d->pre_stride) + ref_col;

    best_mv->row = ref_row;
    best_mv->col = ref_col;

    // We need to check that the starting point for the search (as indicated by ref_mv) is within the buffer limits
    if ((ref_col > x->mv_col_min) && (ref_col < x->mv_col_max) &&
    (ref_row > x->mv_row_min) && (ref_row < x->mv_row_max))
    {
        // Baseline value at the centre
        bestsad = fn_ptr->sdf(what, what_stride, bestaddress, in_what_stride, 0x7fffffff) + mv_err_cost(ref_mv, center_mv, mvsadcost, error_per_bit);
    }

    // Apply further limits to prevent us looking using vectors that stretch beyiond the UMV border
    if (col_min < x->mv_col_min)
        col_min = x->mv_col_min;

    if (col_max > x->mv_col_max)
        col_max = x->mv_col_max;

    if (row_min < x->mv_row_min)
        row_min = x->mv_row_min;

    if (row_max > x->mv_row_max)
        row_max = x->mv_row_max;

    for (r = row_min; r < row_max ; r++)
    {
        this_mv.row = r << 3;
        check_here = r * mv_stride + in_what + col_min;
        c = col_min;

        while ((c + 7) < col_max)
        {
            int i;

            fn_ptr->sdx8f(what, what_stride, check_here , in_what_stride, sad_array8);

            for (i = 0; i < 8; i++)
            {
                thissad = (unsigned int)sad_array8[i];

                if (thissad < bestsad)
                {
                    this_mv.col = c << 3;
                    thissad  += mv_err_cost(&this_mv, center_mv, mvsadcost, error_per_bit);

                    if (thissad < bestsad)
                    {
                        bestsad = thissad;
                        best_mv->row = r;
                        best_mv->col = c;
                        bestaddress = check_here;
                    }
                }

                check_here++;
                c++;
            }
        }

        while ((c + 2) < col_max)
        {
            int i;

            fn_ptr->sdx3f(what, what_stride, check_here , in_what_stride, sad_array);

            for (i = 0; i < 3; i++)
            {
                thissad = sad_array[i];

                if (thissad < bestsad)
                {
                    this_mv.col = c << 3;
                    thissad  += mv_err_cost(&this_mv, center_mv, mvsadcost, error_per_bit);

                    if (thissad < bestsad)
                    {
                        bestsad = thissad;
                        best_mv->row = r;
                        best_mv->col = c;
                        bestaddress = check_here;
                    }
                }

                check_here++;
                c++;
            }
        }

        while (c < col_max)
        {
            thissad = fn_ptr->sdf(what, what_stride, check_here , in_what_stride, bestsad);

            if (thissad < bestsad)
            {
                this_mv.col = c << 3;
                thissad  += mv_err_cost(&this_mv, center_mv, mvsadcost, error_per_bit);

                if (thissad < bestsad)
                {
                    bestsad = thissad;
                    best_mv->row = r;
                    best_mv->col = c;
                    bestaddress = check_here;
                }
            }

            check_here ++;
            c ++;
        }
    }

    this_mv.row = best_mv->row << 3;
    this_mv.col = best_mv->col << 3;

    if (bestsad < INT_MAX)
        return fn_ptr->vf(what, what_stride, bestaddress, in_what_stride, (unsigned int *)(&thissad))
        + mv_err_cost(&this_mv, center_mv, mvcost, error_per_bit);
    else
        return INT_MAX;
}
#endif /* !(CONFIG_REALTIME_ONLY) */

#ifdef ENTROPY_STATS
void print_mode_context(void)
{
    FILE *f = fopen("modecont.c", "w");
    int i, j;

    fprintf(f, "#include \"entropy.h\"\n");
    fprintf(f, "const int vp8_mode_contexts[6][4] =\n");
    fprintf(f, "{\n");

    for (j = 0; j < 6; j++)
    {
        fprintf(f, "  { // %d \n", j);
        fprintf(f, "    ");

        for (i = 0; i < 4; i++)
        {
            int overal_prob;
            int this_prob;
            int count; // = mv_ref_ct[j][i][0]+mv_ref_ct[j][i][1];

            // Overall probs
            count = mv_mode_cts[i][0] + mv_mode_cts[i][1];

            if (count)
                overal_prob = 256 * mv_mode_cts[i][0] / count;
            else
                overal_prob = 128;

            if (overal_prob == 0)
                overal_prob = 1;

            // context probs
            count = mv_ref_ct[j][i][0] + mv_ref_ct[j][i][1];

            if (count)
                this_prob = 256 * mv_ref_ct[j][i][0] / count;
            else
                this_prob = 128;

            if (this_prob == 0)
                this_prob = 1;

            fprintf(f, "%5d, ", this_prob);
            //fprintf(f,"%5d, %5d, %8d,", this_prob, overal_prob, (this_prob << 10)/overal_prob);
            //fprintf(f,"%8d, ", (this_prob << 10)/overal_prob);
        }

        fprintf(f, "  },\n");
    }

    fprintf(f, "};\n");
    fclose(f);
}

/* MV ref count ENTROPY_STATS stats code */
#ifdef ENTROPY_STATS
void init_mv_ref_counts()
{
    vpx_memset(mv_ref_ct, 0, sizeof(mv_ref_ct));
    vpx_memset(mv_mode_cts, 0, sizeof(mv_mode_cts));
}

void accum_mv_refs(MB_PREDICTION_MODE m, const int ct[4])
{
    if (m == ZEROMV)
    {
        ++mv_ref_ct [ct[0]] [0] [0];
        ++mv_mode_cts[0][0];
    }
    else
    {
        ++mv_ref_ct [ct[0]] [0] [1];
        ++mv_mode_cts[0][1];

        if (m == NEARESTMV)
        {
            ++mv_ref_ct [ct[1]] [1] [0];
            ++mv_mode_cts[1][0];
        }
        else
        {
            ++mv_ref_ct [ct[1]] [1] [1];
            ++mv_mode_cts[1][1];

            if (m == NEARMV)
            {
                ++mv_ref_ct [ct[2]] [2] [0];
                ++mv_mode_cts[2][0];
            }
            else
            {
                ++mv_ref_ct [ct[2]] [2] [1];
                ++mv_mode_cts[2][1];

                if (m == NEWMV)
                {
                    ++mv_ref_ct [ct[3]] [3] [0];
                    ++mv_mode_cts[3][0];
                }
                else
                {
                    ++mv_ref_ct [ct[3]] [3] [1];
                    ++mv_mode_cts[3][1];
                }
            }
        }
    }
}

#endif/* END MV ref count ENTROPY_STATS stats code */

#endif
