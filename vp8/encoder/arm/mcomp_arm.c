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

extern unsigned int vp8_sub_pixel_variance16x16s_neon
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    int  xoffset,
    int  yoffset,
    unsigned char *dst_ptr,
    int dst_pixels_per_line,
    unsigned int *sse
);
extern unsigned int vp8_sub_pixel_variance16x16s_4_0_neon
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    unsigned char *dst_ptr,
    int dst_pixels_per_line,
    unsigned int *sse
);
extern unsigned int vp8_sub_pixel_variance16x16s_0_4_neon
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    unsigned char *dst_ptr,
    int dst_pixels_per_line,
    unsigned int *sse
);
extern unsigned int vp8_sub_pixel_variance16x16s_4_4_neon
(
    unsigned char  *src_ptr,
    int  src_pixels_per_line,
    unsigned char *dst_ptr,
    int dst_pixels_per_line,
    unsigned int *sse
);


int vp8_find_best_sub_pixel_step(MACROBLOCK *x, BLOCK *b, BLOCKD *d, MV *bestmv, MV *ref_mv, int error_per_bit, vp8_subpixvariance_fn_t svf, vp8_variance_fn_t vf, int *mvcost[2])
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
    bestmse = vf(y, d->pre_stride, z, b->src_stride, &sse);
    bestmse += vp8_mv_err_cost(bestmv, ref_mv, mvcost, error_per_bit);

    // go left then right and check error
    this_mv.row = startmv.row;
    this_mv.col = ((startmv.col - 8) | 4);
    left = vp8_sub_pixel_variance16x16s_4_0_neon(y - 1, d->pre_stride, z, b->src_stride, &sse);
    left += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (left < bestmse)
    {
        *bestmv = this_mv;
        bestmse = left;
    }

    this_mv.col += 8;
    right = vp8_sub_pixel_variance16x16s_4_0_neon(y, d->pre_stride, z, b->src_stride, &sse);
    right += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (right < bestmse)
    {
        *bestmv = this_mv;
        bestmse = right;
    }

    // go up then down and check error
    this_mv.col = startmv.col;
    this_mv.row = ((startmv.row - 8) | 4);
    up = vp8_sub_pixel_variance16x16s_0_4_neon(y - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
    up += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (up < bestmse)
    {
        *bestmv = this_mv;
        bestmse = up;
    }

    this_mv.row += 8;
    down = vp8_sub_pixel_variance16x16s_0_4_neon(y, d->pre_stride, z, b->src_stride, &sse);
    down += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

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
        diag = vp8_sub_pixel_variance16x16s_4_4_neon(y - 1 - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
        break;
    case 1:
        this_mv.col += 4;
        this_mv.row = (this_mv.row - 8) | 4;
        diag = vp8_sub_pixel_variance16x16s_4_4_neon(y - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
        break;
    case 2:
        this_mv.col = (this_mv.col - 8) | 4;
        this_mv.row += 4;
        diag = vp8_sub_pixel_variance16x16s_4_4_neon(y - 1, d->pre_stride, z, b->src_stride, &sse);
        break;
    case 3:
        this_mv.col += 4;
        this_mv.row += 4;
        diag = vp8_sub_pixel_variance16x16s_4_4_neon(y, d->pre_stride, z, b->src_stride, &sse);
        break;
    }

    diag += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

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
        left = svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
    }
    else
    {
        this_mv.col = (startmv.col - 8) | 6;
        left = svf(y - 1, d->pre_stride, 6, this_mv.row & 7, z, b->src_stride, &sse);
    }

    left += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (left < bestmse)
    {
        *bestmv = this_mv;
        bestmse = left;
    }

    this_mv.col += 4;
    right = svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
    right += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

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
        up = svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
    }
    else
    {
        this_mv.row = (startmv.row - 8) | 6;
        up = svf(y - d->pre_stride, d->pre_stride, this_mv.col & 7, 6, z, b->src_stride, &sse);
    }

    up += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (up < bestmse)
    {
        *bestmv = this_mv;
        bestmse = up;
    }

    this_mv.row += 4;
    down = svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
    down += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

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
                diag = svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
            }
            else
            {
                this_mv.col = (startmv.col - 8) | 6;
                diag = svf(y - 1, d->pre_stride, 6, this_mv.row & 7, z, b->src_stride, &sse);;
            }
        }
        else
        {
            this_mv.row = (startmv.row - 8) | 6;

            if (startmv.col & 7)
            {
                this_mv.col -= 2;
                diag = svf(y - d->pre_stride, d->pre_stride, this_mv.col & 7, 6, z, b->src_stride, &sse);
            }
            else
            {
                this_mv.col = (startmv.col - 8) | 6;
                diag = svf(y - d->pre_stride - 1, d->pre_stride, 6, 6, z, b->src_stride, &sse);
            }
        }

        break;
    case 1:
        this_mv.col += 2;

        if (startmv.row & 7)
        {
            this_mv.row -= 2;
            diag = svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
        }
        else
        {
            this_mv.row = (startmv.row - 8) | 6;
            diag = svf(y - d->pre_stride, d->pre_stride, this_mv.col & 7, 6, z, b->src_stride, &sse);
        }

        break;
    case 2:
        this_mv.row += 2;

        if (startmv.col & 7)
        {
            this_mv.col -= 2;
            diag = svf(y, d->pre_stride, this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
        }
        else
        {
            this_mv.col = (startmv.col - 8) | 6;
            diag = svf(y - 1, d->pre_stride, 6, this_mv.row & 7, z, b->src_stride, &sse);;
        }

        break;
    case 3:
        this_mv.col += 2;
        this_mv.row += 2;
        diag = svf(y, d->pre_stride,  this_mv.col & 7, this_mv.row & 7, z, b->src_stride, &sse);
        break;
    }

    diag += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

//  }

    return bestmse;
}

int vp8_find_best_half_pixel_step(MACROBLOCK *mb, BLOCK *b, BLOCKD *d, MV *bestmv, MV *ref_mv, int error_per_bit, vp8_subpixvariance_fn_t svf, vp8_variance_fn_t vf, int *mvcost[2])
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
    bestmse = vf(y, d->pre_stride, z, b->src_stride, &sse);
    bestmse += vp8_mv_err_cost(bestmv, ref_mv, mvcost, error_per_bit);

    // go left then right and check error
    this_mv.row = startmv.row;
    this_mv.col = ((startmv.col - 8) | 4);
    left = vp8_sub_pixel_variance16x16s_4_0_neon(y - 1, d->pre_stride, z, b->src_stride, &sse);
    left += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (left < bestmse)
    {
        *bestmv = this_mv;
        bestmse = left;
    }

    this_mv.col += 8;
    right = vp8_sub_pixel_variance16x16s_4_0_neon(y, d->pre_stride, z, b->src_stride, &sse);
    right += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (right < bestmse)
    {
        *bestmv = this_mv;
        bestmse = right;
    }

    // go up then down and check error
    this_mv.col = startmv.col;
    this_mv.row = ((startmv.row - 8) | 4);
    up = vp8_sub_pixel_variance16x16s_0_4_neon(y - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
    up += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (up < bestmse)
    {
        *bestmv = this_mv;
        bestmse = up;
    }

    this_mv.row += 8;
    down = vp8_sub_pixel_variance16x16s_0_4_neon(y, d->pre_stride, z, b->src_stride, &sse);
    down += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

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
        diag = svf(y - 1 - d->pre_stride, d->pre_stride, 4, 4, z, b->src_stride, &sse);
        break;
    case 1:
        this_mv.col += 4;
        this_mv.row = (this_mv.row - 8) | 4;
        diag = svf(y - d->pre_stride, d->pre_stride, 4, 4, z, b->src_stride, &sse);
        break;
    case 2:
        this_mv.col = (this_mv.col - 8) | 4;
        this_mv.row += 4;
        diag = svf(y - 1, d->pre_stride, 4, 4, z, b->src_stride, &sse);
        break;
    case 3:
        this_mv.col += 4;
        this_mv.row += 4;
        diag = svf(y, d->pre_stride, 4, 4, z, b->src_stride, &sse);
        break;
    }

    diag += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

#else
    this_mv.col = (this_mv.col - 8) | 4;
    this_mv.row = (this_mv.row - 8) | 4;
    diag = vp8_sub_pixel_variance16x16s_4_4_neon(y - 1 - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
    diag += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

    this_mv.col += 8;
    diag = vp8_sub_pixel_variance16x16s_4_4_neon(y - d->pre_stride, d->pre_stride, z, b->src_stride, &sse);
    diag += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

    this_mv.col = (this_mv.col - 8) | 4;
    this_mv.row = startmv.row + 4;
    diag = vp8_sub_pixel_variance16x16s_4_4_neon(y - 1, d->pre_stride, z, b->src_stride, &sse);
    diag += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

    if (diag < bestmse)
    {
        *bestmv = this_mv;
        bestmse = diag;
    }

    this_mv.col += 8;
    diag = vp8_sub_pixel_variance16x16s_4_4_neon(y, d->pre_stride, z, b->src_stride, &sse);
    diag += vp8_mv_err_cost(&this_mv, ref_mv, mvcost, error_per_bit);

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
#define DIST(r,c,v) sf( src,src_stride,PRE(r,c),d->pre_stride, v) // returns sad error score.
#define ERR(r,c,v) (MVC(r,c)+DIST(r,c,v)) // returns distortion + motion vector cost
#define CHECK_BETTER(v,r,c) if ((v = ERR(r,c,besterr)) < besterr) { besterr = v; br=r; bc=c; } // checks if (r,c) has better score than previous best
const MV next_chkpts[6][3] =
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
    vp8_variance_fn_t vf,
    vp8_sad_fn_t      sf,
    int *mvsadcost[2],
    int *mvcost[2]
)
{
    MV hex[6] = { { -1, -2}, {1, -2}, {2, 0}, {1, 2}, { -1, 2}, { -2, 0} } ;
    MV neighbors[8] = { { -1, -1}, { -1, 0}, { -1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1} } ;
    int i, j;
    unsigned char *src = (*(b->base_src) + b->src);
    int src_stride = b->src_stride;
    int rr = ref_mv->row, rc = ref_mv->col, br = rr >> 3, bc = rc >> 3, tr, tc;
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

    return vf(src, src_stride, PRE(br, bc), d->pre_stride, &thiserr) + MVC(br, bc) ;
}
#undef MVC
#undef PRE
#undef SP
#undef DIST
#undef ERR
#undef CHECK_BETTER


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
