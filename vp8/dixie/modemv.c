/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "dixie.h"
#include <stdlib.h>
#include <assert.h>

static const unsigned char kf_y_mode_probs[]  = { 145, 156, 163, 128};
static const unsigned char kf_uv_mode_probs[] = { 142, 114, 183};
static const unsigned char kf_b_mode_probs[10][10][9] =
{
    { /* above mode 0 */
        { /* left mode 0 */ 231, 120,  48,  89, 115, 113, 120, 152, 112},
        { /* left mode 1 */ 152, 179,  64, 126, 170, 118,  46,  70,  95},
        { /* left mode 2 */ 175,  69, 143,  80,  85,  82,  72, 155, 103},
        { /* left mode 3 */  56,  58,  10, 171, 218, 189,  17,  13, 152},
        { /* left mode 4 */ 144,  71,  10,  38, 171, 213, 144,  34,  26},
        { /* left mode 5 */ 114,  26,  17, 163,  44, 195,  21,  10, 173},
        { /* left mode 6 */ 121,  24,  80, 195,  26,  62,  44,  64,  85},
        { /* left mode 7 */ 170,  46,  55,  19, 136, 160,  33, 206,  71},
        { /* left mode 8 */  63,  20,   8, 114, 114, 208,  12,   9, 226},
        { /* left mode 9 */  81,  40,  11,  96, 182,  84,  29,  16,  36}
    },
    { /* above mode 1 */
        { /* left mode 0 */ 134, 183,  89, 137,  98, 101, 106, 165, 148},
        { /* left mode 1 */  72, 187, 100, 130, 157, 111,  32,  75,  80},
        { /* left mode 2 */  66, 102, 167,  99,  74,  62,  40, 234, 128},
        { /* left mode 3 */  41,  53,   9, 178, 241, 141,  26,   8, 107},
        { /* left mode 4 */ 104,  79,  12,  27, 217, 255,  87,  17,   7},
        { /* left mode 5 */  74,  43,  26, 146,  73, 166,  49,  23, 157},
        { /* left mode 6 */  65,  38, 105, 160,  51,  52,  31, 115, 128},
        { /* left mode 7 */  87,  68,  71,  44, 114,  51,  15, 186,  23},
        { /* left mode 8 */  47,  41,  14, 110, 182, 183,  21,  17, 194},
        { /* left mode 9 */  66,  45,  25, 102, 197, 189,  23,  18,  22}
    },
    { /* above mode 2 */
        { /* left mode 0 */  88,  88, 147, 150,  42,  46,  45, 196, 205},
        { /* left mode 1 */  43,  97, 183, 117,  85,  38,  35, 179,  61},
        { /* left mode 2 */  39,  53, 200,  87,  26,  21,  43, 232, 171},
        { /* left mode 3 */  56,  34,  51, 104, 114, 102,  29,  93,  77},
        { /* left mode 4 */ 107,  54,  32,  26,  51,   1,  81,  43,  31},
        { /* left mode 5 */  39,  28,  85, 171,  58, 165,  90,  98,  64},
        { /* left mode 6 */  34,  22, 116, 206,  23,  34,  43, 166,  73},
        { /* left mode 7 */  68,  25, 106,  22,  64, 171,  36, 225, 114},
        { /* left mode 8 */  34,  19,  21, 102, 132, 188,  16,  76, 124},
        { /* left mode 9 */  62,  18,  78,  95,  85,  57,  50,  48,  51}
    },
    { /* above mode 3 */
        { /* left mode 0 */ 193, 101,  35, 159, 215, 111,  89,  46, 111},
        { /* left mode 1 */  60, 148,  31, 172, 219, 228,  21,  18, 111},
        { /* left mode 2 */ 112, 113,  77,  85, 179, 255,  38, 120, 114},
        { /* left mode 3 */  40,  42,   1, 196, 245, 209,  10,  25, 109},
        { /* left mode 4 */ 100,  80,   8,  43, 154,   1,  51,  26,  71},
        { /* left mode 5 */  88,  43,  29, 140, 166, 213,  37,  43, 154},
        { /* left mode 6 */  61,  63,  30, 155,  67,  45,  68,   1, 209},
        { /* left mode 7 */ 142,  78,  78,  16, 255, 128,  34, 197, 171},
        { /* left mode 8 */  41,  40,   5, 102, 211, 183,   4,   1, 221},
        { /* left mode 9 */  51,  50,  17, 168, 209, 192,  23,  25,  82}
    },
    { /* above mode 4 */
        { /* left mode 0 */ 125,  98,  42,  88, 104,  85, 117, 175,  82},
        { /* left mode 1 */  95,  84,  53,  89, 128, 100, 113, 101,  45},
        { /* left mode 2 */  75,  79, 123,  47,  51, 128,  81, 171,   1},
        { /* left mode 3 */  57,  17,   5,  71, 102,  57,  53,  41,  49},
        { /* left mode 4 */ 115,  21,   2,  10, 102, 255, 166,  23,   6},
        { /* left mode 5 */  38,  33,  13, 121,  57,  73,  26,   1,  85},
        { /* left mode 6 */  41,  10,  67, 138,  77, 110,  90,  47, 114},
        { /* left mode 7 */ 101,  29,  16,  10,  85, 128, 101, 196,  26},
        { /* left mode 8 */  57,  18,  10, 102, 102, 213,  34,  20,  43},
        { /* left mode 9 */ 117,  20,  15,  36, 163, 128,  68,   1,  26}
    },
    { /* above mode 5 */
        { /* left mode 0 */ 138,  31,  36, 171,  27, 166,  38,  44, 229},
        { /* left mode 1 */  67,  87,  58, 169,  82, 115,  26,  59, 179},
        { /* left mode 2 */  63,  59,  90, 180,  59, 166,  93,  73, 154},
        { /* left mode 3 */  40,  40,  21, 116, 143, 209,  34,  39, 175},
        { /* left mode 4 */  57,  46,  22,  24, 128,   1,  54,  17,  37},
        { /* left mode 5 */  47,  15,  16, 183,  34, 223,  49,  45, 183},
        { /* left mode 6 */  46,  17,  33, 183,   6,  98,  15,  32, 183},
        { /* left mode 7 */  65,  32,  73, 115,  28, 128,  23, 128, 205},
        { /* left mode 8 */  40,   3,   9, 115,  51, 192,  18,   6, 223},
        { /* left mode 9 */  87,  37,   9, 115,  59,  77,  64,  21,  47}
    },
    { /* above mode 6 */
        { /* left mode 0 */ 104,  55,  44, 218,   9,  54,  53, 130, 226},
        { /* left mode 1 */  64,  90,  70, 205,  40,  41,  23,  26,  57},
        { /* left mode 2 */  54,  57, 112, 184,   5,  41,  38, 166, 213},
        { /* left mode 3 */  30,  34,  26, 133, 152, 116,  10,  32, 134},
        { /* left mode 4 */  75,  32,  12,  51, 192, 255, 160,  43,  51},
        { /* left mode 5 */  39,  19,  53, 221,  26, 114,  32,  73, 255},
        { /* left mode 6 */  31,   9,  65, 234,   2,  15,   1, 118,  73},
        { /* left mode 7 */  88,  31,  35,  67, 102,  85,  55, 186,  85},
        { /* left mode 8 */  56,  21,  23, 111,  59, 205,  45,  37, 192},
        { /* left mode 9 */  55,  38,  70, 124,  73, 102,   1,  34,  98}
    },
    { /* above mode 7 */
        { /* left mode 0 */ 102,  61,  71,  37,  34,  53,  31, 243, 192},
        { /* left mode 1 */  69,  60,  71,  38,  73, 119,  28, 222,  37},
        { /* left mode 2 */  68,  45, 128,  34,   1,  47,  11, 245, 171},
        { /* left mode 3 */  62,  17,  19,  70, 146,  85,  55,  62,  70},
        { /* left mode 4 */  75,  15,   9,   9,  64, 255, 184, 119,  16},
        { /* left mode 5 */  37,  43,  37, 154, 100, 163,  85, 160,   1},
        { /* left mode 6 */  63,   9,  92, 136,  28,  64,  32, 201,  85},
        { /* left mode 7 */  86,   6,  28,   5,  64, 255,  25, 248,   1},
        { /* left mode 8 */  56,   8,  17, 132, 137, 255,  55, 116, 128},
        { /* left mode 9 */  58,  15,  20,  82, 135,  57,  26, 121,  40}
    },
    { /* above mode 8 */
        { /* left mode 0 */ 164,  50,  31, 137, 154, 133,  25,  35, 218},
        { /* left mode 1 */  51, 103,  44, 131, 131, 123,  31,   6, 158},
        { /* left mode 2 */  86,  40,  64, 135, 148, 224,  45, 183, 128},
        { /* left mode 3 */  22,  26,  17, 131, 240, 154,  14,   1, 209},
        { /* left mode 4 */  83,  12,  13,  54, 192, 255,  68,  47,  28},
        { /* left mode 5 */  45,  16,  21,  91,  64, 222,   7,   1, 197},
        { /* left mode 6 */  56,  21,  39, 155,  60, 138,  23, 102, 213},
        { /* left mode 7 */  85,  26,  85,  85, 128, 128,  32, 146, 171},
        { /* left mode 8 */  18,  11,   7,  63, 144, 171,   4,   4, 246},
        { /* left mode 9 */  35,  27,  10, 146, 174, 171,  12,  26, 128}
    },
    { /* above mode 9 */
        { /* left mode 0 */ 190,  80,  35,  99, 180,  80, 126,  54,  45},
        { /* left mode 1 */  85, 126,  47,  87, 176,  51,  41,  20,  32},
        { /* left mode 2 */ 101,  75, 128, 139, 118, 146, 116, 128,  85},
        { /* left mode 3 */  56,  41,  15, 176, 236,  85,  37,   9,  62},
        { /* left mode 4 */ 146,  36,  19,  30, 171, 255,  97,  27,  20},
        { /* left mode 5 */  71,  30,  17, 119, 118, 255,  17,  18, 138},
        { /* left mode 6 */ 101,  38,  60, 138,  55,  70,  43,  26, 142},
        { /* left mode 7 */ 138,  45,  61,  62, 219,   1,  81, 188,  64},
        { /* left mode 8 */  32,  41,  20, 117, 151, 142,  20,  21, 163},
        { /* left mode 9 */ 112,  19,  12,  61, 195, 128,  48,   4,  24}
    }
};
static const int kf_y_mode_tree[] =
{
    -B_PRED, 2,
    4, 6,
    -DC_PRED, -V_PRED,
    -H_PRED, -TM_PRED
};
static const int y_mode_tree[] =
{
    -DC_PRED, 2,
    4, 6,
    -V_PRED, -H_PRED,
    -TM_PRED, -B_PRED
};
static const int uv_mode_tree[6] =
{
    -DC_PRED, 2,
    -V_PRED, 4,
    -H_PRED, -TM_PRED
};
static const int b_mode_tree[18] =
{
    -B_DC_PRED, 2,                             /* 0 = DC_NODE */
    -B_TM_PRED, 4,                            /* 1 = TM_NODE */
    -B_VE_PRED, 6,                           /* 2 = VE_NODE */
    8, 12,                                  /* 3 = COM_NODE */
    -B_HE_PRED, 10,                        /* 4 = HE_NODE */
    -B_RD_PRED, -B_VR_PRED,               /* 5 = RD_NODE */
    -B_LD_PRED, 14,                        /* 6 = LD_NODE */
    -B_VL_PRED, 16,                      /* 7 = VL_NODE */
    -B_HD_PRED, -B_HU_PRED             /* 8 = HD_NODE */
};
static const int small_mv_tree[14] =
{
    2, 8,
    4, 6,
    -0, -1,
    -2, -3,
    10, 12,
    -4, -5,
    -6, -7
};
static const int mv_ref_tree[8] =
{
    -ZEROMV, 2,
    -NEARESTMV, 4,
    -NEARMV, 6,
    -NEWMV, -SPLITMV
};
static const int submv_ref_tree[6] =
{
    -LEFT4X4, 2,
    -ABOVE4X4, 4,
    -ZERO4X4, -NEW4X4
};
static const int split_mv_tree[6] =
{
    -3, 2,
    -2, 4,
    -0, -1
};
static const unsigned char default_b_mode_probs[] =
{ 120,  90,  79, 133,  87,  85,  80, 111, 151};
static const unsigned char mv_counts_to_probs[6][4] =
{
    {   7,   1,   1, 143 },
    {  14,  18,  14, 107 },
    { 135,  64,  57,  68 },
    {  60,  56, 128,  65 },
    { 159, 134, 128,  34 },
    { 234, 188, 128,  28 }

};
static const unsigned char split_mv_probs[3] =
{ 110, 111, 150};
static const unsigned char submv_ref_probs2[5][3] =
{
    { 147, 136, 18 },
    { 106, 145,  1 },
    { 179, 121,  1 },
    { 223,   1, 34 },
    { 208,   1,  1 }
};

const static int mv_partitions[4][16] =
{
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  1,  1,  1,  1,  1,  1 },
    {0, 0, 1, 1, 0, 0, 1, 1, 0, 0,  1,  1,  0,  0,  1,  1 },
    {0, 0, 1, 1, 0, 0, 1, 1, 2, 2,  3,  3,  2,  2,  3,  3 },
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }
};


struct mv_clamp_rect
{
    int to_left, to_right, to_top, to_bottom;
};


static union mv
        clamp_mv(union mv raw, const struct mv_clamp_rect *bounds)
{
    union mv newmv;

    newmv.d.x = (raw.d.x < bounds->to_left) ? bounds->to_left : raw.d.x;
    newmv.d.x = (raw.d.x > bounds->to_right) ? bounds->to_right : newmv.d.x;
    newmv.d.y = (raw.d.y < bounds->to_top) ? bounds->to_top : raw.d.y;
    newmv.d.y = (raw.d.y > bounds->to_bottom) ? bounds->to_bottom : newmv.d.y;
    return newmv;
}


static int
read_segment_id(struct bool_decoder *bool, struct vp8_segment_hdr *seg)
{
    return bool_get(bool, seg->tree_probs[0])
           ? 2 + bool_get(bool, seg->tree_probs[2])
           : bool_get(bool, seg->tree_probs[1]);
}


static enum prediction_mode
above_block_mode(const struct mb_info *this,
                 const struct mb_info *above,
                 unsigned int b)
{
    if (b < 4)
    {
        switch (above->base.y_mode)
        {
        case DC_PRED:
            return B_DC_PRED;
        case V_PRED:
            return B_VE_PRED;
        case H_PRED:
            return B_HE_PRED;
        case TM_PRED:
            return B_TM_PRED;
        case B_PRED:
            return above->split.modes[b+12];
        default:
            assert(0);
        }
    }

    return this->split.modes[b-4];
}


static enum prediction_mode
left_block_mode(const struct mb_info *this,
                const struct mb_info *left,
                unsigned int b)
{
    if (!(b & 3))
    {
        switch (left->base.y_mode)
        {
        case DC_PRED:
            return B_DC_PRED;
        case V_PRED:
            return B_VE_PRED;
        case H_PRED:
            return B_HE_PRED;
        case TM_PRED:
            return B_TM_PRED;
        case B_PRED:
            return left->split.modes[b+3];
        default:
            assert(0);
        }
    }

    return this->split.modes[b-1];
}


static void
decode_kf_mb_mode(struct mb_info      *this,
                  struct mb_info      *left,
                  struct mb_info      *above,
                  struct bool_decoder *bool)
{
    int y_mode, uv_mode;

    y_mode = bool_read_tree(bool, kf_y_mode_tree, kf_y_mode_probs);

    if (y_mode == B_PRED)
    {
        unsigned int i;

        for (i = 0; i < 16; i++)
        {
            enum prediction_mode a = above_block_mode(this, above, i);
            enum prediction_mode l = left_block_mode(this, left, i);
            enum prediction_mode b;

            b = bool_read_tree(bool, b_mode_tree, kf_b_mode_probs[a][l]);
            this->split.modes[i] = b;
        }
    }

    uv_mode = bool_read_tree(bool, uv_mode_tree, kf_uv_mode_probs);

    this->base.y_mode = y_mode;
    this->base.uv_mode = uv_mode;
    this->base.mv.raw = 0;
    this->base.ref_frame = 0;
}


static void
decode_intra_mb_mode(struct mb_info         *this,
                     struct vp8_entropy_hdr *hdr,
                     struct bool_decoder    *bool)
{
    /* Like decode_kf_mb_mode, but with probabilities transmitted in the
     * bitstream and no context on the above/left block mode.
     */
    int y_mode, uv_mode;

    y_mode = bool_read_tree(bool, y_mode_tree, hdr->y_mode_probs);

    if (y_mode == B_PRED)
    {
        unsigned int i;

        for (i = 0; i < 16; i++)
        {
            enum prediction_mode b;

            b = bool_read_tree(bool, b_mode_tree, default_b_mode_probs);
            this->split.modes[i] = b;
        }
    }

    uv_mode = bool_read_tree(bool, uv_mode_tree, hdr->uv_mode_probs);

    this->base.y_mode = y_mode;
    this->base.uv_mode = uv_mode;
    this->base.mv.raw = 0;
    this->base.ref_frame = CURRENT_FRAME;
}


static int
read_mv_component(struct bool_decoder *bool,
                  const unsigned char  mvc[MV_PROB_CNT])
{
    enum {IS_SHORT, SIGN, SHORT, BITS = SHORT + 8 - 1, LONG_WIDTH = 10};
    int x = 0;

    if (bool_get(bool, mvc[IS_SHORT])) /* Large */
    {
        int i = 0;

        for (i = 0; i < 3; i++)
            x += bool_get(bool, mvc[BITS + i]) << i;

        /* Skip bit 3, which is sometimes implicit */
        for (i = LONG_WIDTH - 1; i > 3; i--)
            x += bool_get(bool, mvc[BITS + i]) << i;

        if (!(x & 0xFFF0)  ||  bool_get(bool, mvc[BITS + 3]))
            x += 8;
    }
    else   /* small */
        x = bool_read_tree(bool, small_mv_tree, mvc + SHORT);

    if (x && bool_get(bool, mvc[SIGN]))
        x = -x;

    return x << 1;
}


static mv_t
above_block_mv(const struct mb_info *this,
               const struct mb_info *above,
               unsigned int          b)
{
    if (b < 4)
    {
        if (above->base.y_mode == SPLITMV)
            return above->split.mvs[b+12];

        return above->base.mv;
    }

    return this->split.mvs[b-4];
}


static mv_t
left_block_mv(const struct mb_info *this,
              const struct mb_info *left,
              unsigned int          b)
{
    if (!(b & 3))
    {
        if (left->base.y_mode == SPLITMV)
            return left->split.mvs[b+3];

        return left->base.mv;
    }

    return this->split.mvs[b-1];
}


static enum prediction_mode
submv_ref(struct bool_decoder *bool, union mv l, union mv a)
{
    enum subblock_mv_ref
    {
        SUBMVREF_NORMAL,
        SUBMVREF_LEFT_ZED,
        SUBMVREF_ABOVE_ZED,
        SUBMVREF_LEFT_ABOVE_SAME,
        SUBMVREF_LEFT_ABOVE_ZED
    };

    int lez = !(l.raw);
    int aez = !(a.raw);
    int lea = l.raw == a.raw;
    enum subblock_mv_ref ctx = SUBMVREF_NORMAL;

    if (lea && lez)
        ctx = SUBMVREF_LEFT_ABOVE_ZED;
    else if (lea)
        ctx = SUBMVREF_LEFT_ABOVE_SAME;
    else if (aez)
        ctx = SUBMVREF_ABOVE_ZED;
    else if (lez)
        ctx = SUBMVREF_LEFT_ZED;

    return bool_read_tree(bool, submv_ref_tree, submv_ref_probs2[ctx]);
}


static void
read_mv(struct bool_decoder  *bool,
        union mv             *mv,
        mv_component_probs_t  mvc[2])
{
    mv->d.y = read_mv_component(bool, mvc[0]);
    mv->d.x = read_mv_component(bool, mvc[1]);
}


static void
mv_bias(const struct mb_info *mb,
        const unsigned int   sign_bias[3],
        enum reference_frame ref_frame,
        union mv             *mv)
{
    if (sign_bias[mb->base.ref_frame] ^ sign_bias[ref_frame])
    {
        mv->d.x *= -1;
        mv->d.y *= -1;
    }
}


enum near_mv_v
{
    CNT_BEST = 0,
    CNT_ZEROZERO = 0,
    CNT_NEAREST,
    CNT_NEAR,
    CNT_SPLITMV
};


static void
find_near_mvs(const struct mb_info   *this,
              const struct mb_info   *left,
              const struct mb_info   *above,
              const unsigned int      sign_bias[3],
              union  mv               near_mvs[4],
              int                     cnt[4])
{
    const struct mb_info *aboveleft = above - 1;
    union  mv             *mv = near_mvs;
    int                   *cntx = cnt;

    /* Zero accumulators */
    mv[0].raw = mv[1].raw = mv[2].raw = 0;
    cnt[0] = cnt[1] = cnt[2] = cnt[3] = 0;

    /* Process above */
    if (above->base.ref_frame != CURRENT_FRAME)
    {
        if (above->base.mv.raw)
        {
            (++mv)->raw = above->base.mv.raw;
            mv_bias(above, sign_bias, this->base.ref_frame, mv);
            ++cntx;
        }

        *cntx += 2;
    }

    /* Process left */
    if (left->base.ref_frame != CURRENT_FRAME)
    {
        if (left->base.mv.raw)
        {
            union mv this_mv;

            this_mv.raw = left->base.mv.raw;
            mv_bias(left, sign_bias, this->base.ref_frame, &this_mv);

            if (this_mv.raw != mv->raw)
            {
                (++mv)->raw = this_mv.raw;
                ++cntx;
            }

            *cntx += 2;
        }
        else
            cnt[CNT_ZEROZERO] += 2;
    }

    /* Process above left */
    if (aboveleft->base.ref_frame != CURRENT_FRAME)
    {
        if (aboveleft->base.mv.raw)
        {
            union mv this_mv;

            this_mv.raw = aboveleft->base.mv.raw;
            mv_bias(aboveleft, sign_bias, this->base.ref_frame, &this_mv);

            if (this_mv.raw != mv->raw)
            {
                (++mv)->raw = this_mv.raw;
                ++cntx;
            }

            *cntx += 1;
        }
        else
            cnt[CNT_ZEROZERO] += 1;
    }

    /* If we have three distinct MV's ... */
    if (cnt[CNT_SPLITMV])
    {
        /* See if above-left MV can be merged with NEAREST */
        if (mv->raw == near_mvs[CNT_NEAREST].raw)
            cnt[CNT_NEAREST] += 1;
    }

    cnt[CNT_SPLITMV] = ((above->base.y_mode == SPLITMV)
                        + (left->base.y_mode == SPLITMV)) * 2
                       + (aboveleft->base.y_mode == SPLITMV);

    /* Swap near and nearest if necessary */
    if (cnt[CNT_NEAR] > cnt[CNT_NEAREST])
    {
        int tmp;
        tmp = cnt[CNT_NEAREST];
        cnt[CNT_NEAREST] = cnt[CNT_NEAR];
        cnt[CNT_NEAR] = tmp;
        tmp = near_mvs[CNT_NEAREST].raw;
        near_mvs[CNT_NEAREST].raw = near_mvs[CNT_NEAR].raw;
        near_mvs[CNT_NEAR].raw = tmp;
    }

    /* Use near_mvs[CNT_BEST] to store the "best" MV. Note that this storage
     * shares the same address as near_mvs[CNT_ZEROZERO].
     */
    if (cnt[CNT_NEAREST] >= cnt[CNT_BEST])
        near_mvs[CNT_BEST] = near_mvs[CNT_NEAREST];
}


static void
decode_split_mv(struct mb_info         *this,
                const struct mb_info   *left,
                const struct mb_info   *above,
                struct vp8_entropy_hdr *hdr,
                union  mv              *best_mv,
                struct bool_decoder    *bool)
{
    const int *partition;
    int        j, k, mask, partition_id;

    partition_id = bool_read_tree(bool, split_mv_tree, split_mv_probs);
    partition = mv_partitions[partition_id];

    for (j = 0, mask = 0; mask < 65535; j++)
    {
        union mv mv, left_mv, above_mv;
        enum prediction_mode subblock_mode;

        /* Find the first subblock in this partition. */
        for (k = 0; j != partition[k]; k++);

        /* Decode the next MV */
        left_mv = left_block_mv(this, left, k);
        above_mv = above_block_mv(this, above, k);
        subblock_mode = submv_ref(bool, left_mv,  above_mv);

        switch (subblock_mode)
        {
        case LEFT4X4:
            mv = left_mv;
            break;
        case ABOVE4X4:
            mv = above_mv;
            break;
        case ZERO4X4:
            mv.raw = 0;
            break;
        case NEW4X4:
            read_mv(bool, &mv, hdr->mv_probs);
            mv.d.x += best_mv->d.x;
            mv.d.y += best_mv->d.y;
            break;
        default:
            assert(0);
        }

        /* Fill the MV's for this partition */
        for (; k < 16; k++)
            if (j == partition[k])
            {
                this->split.mvs[k] = mv;
                mask |= 1 << k;
            }
    }
}


static void
decode_mvs(struct vp8_decoder_ctx       *ctx,
           struct mb_info               *this,
           const struct mb_info         *left,
           const struct mb_info         *above,
           const struct mv_clamp_rect   *bounds,
           struct bool_decoder          *bool)
{
    struct vp8_entropy_hdr *hdr = &ctx->entropy_hdr;
    union mv          near_mvs[4];
    union mv          clamped_best_mv;
    int               mv_cnts[4];
    unsigned char     probs[4];
    enum {BEST, NEAREST, NEAR};
    int ref_frame;

    this->base.ref_frame = bool_get(bool, hdr->prob_last)
                           ? 2 + bool_get(bool, hdr->prob_gf)
                           : 1;

    find_near_mvs(this, this - 1, above, ctx->reference_hdr.sign_bias,
                  near_mvs, mv_cnts);
    probs[0] = mv_counts_to_probs[mv_cnts[0]][0];
    probs[1] = mv_counts_to_probs[mv_cnts[1]][1];
    probs[2] = mv_counts_to_probs[mv_cnts[2]][2];
    probs[3] = mv_counts_to_probs[mv_cnts[3]][3];

    this->base.y_mode = bool_read_tree(bool, mv_ref_tree, probs);
    this->base.uv_mode = this->base.y_mode;

    switch (this->base.y_mode)
    {
    case NEARESTMV:
        this->base.mv = clamp_mv(near_mvs[NEAREST], bounds);
        break;
    case NEARMV:
        this->base.mv = clamp_mv(near_mvs[NEAR], bounds);
        break;
    case ZEROMV:
        this->base.mv.raw = 0;
        break;
    case NEWMV:
        clamped_best_mv = clamp_mv(near_mvs[BEST], bounds);
        read_mv(bool, &this->base.mv, hdr->mv_probs);
        this->base.mv.d.x += clamped_best_mv.d.x;
        this->base.mv.d.y += clamped_best_mv.d.y;
        break;
    case SPLITMV:
        clamped_best_mv = clamp_mv(near_mvs[BEST], bounds);
        decode_split_mv(this, left, above, hdr, &clamped_best_mv, bool);
        this->base.mv = this->split.mvs[15];
        break;
    default:
        assert(0);
    }
}


void
vp8_dixie_modemv_process_row(struct vp8_decoder_ctx *ctx,
                             struct bool_decoder    *bool,
                             unsigned int            row,
                             unsigned int            start_col,
                             unsigned int            num_cols)
{
    struct mb_info       *above, *this;
    unsigned int          col;
    struct mv_clamp_rect  bounds;

    this = ctx->mb_info_rows[row];
    above = ctx->mb_info_rows[row - 1];

    /* Calculate the eighth-pel MV bounds using a 1 MB border. */
    bounds.to_left   = -((start_col + 1) << 7);
    bounds.to_right  = (ctx->mb_cols - start_col) << 7;
    bounds.to_top    = -((row + 1) << 7);
    bounds.to_bottom = (ctx->mb_rows - row) << 7;

    for (col = start_col; col < start_col + num_cols; col++)
    {
        if (ctx->segment_hdr.update_map)
            this->base.segment_id = read_segment_id(bool,
                                                    &ctx->segment_hdr);

        if (ctx->entropy_hdr.coeff_skip_enabled)
            this->base.skip_coeff = bool_get(bool,
                                             ctx->entropy_hdr.coeff_skip_prob);

        if (ctx->frame_hdr.is_keyframe)
        {
            if (!ctx->segment_hdr.update_map)
                this->base.segment_id = 0;

            decode_kf_mb_mode(this, this - 1, above, bool);
        }
        else
        {
            if (bool_get(bool, ctx->entropy_hdr.prob_inter))
                decode_mvs(ctx, this, this - 1, above, &bounds, bool);
            else
                decode_intra_mb_mode(this, &ctx->entropy_hdr, bool);

            bounds.to_left -= 16 << 3;
            bounds.to_right -= 16 << 3;
        }

        /* Advance to next mb */
        this++;
        above++;
    }
}


void
vp8_dixie_modemv_init(struct vp8_decoder_ctx *ctx)
{
    unsigned int    mbi_w, mbi_h, i;
    struct mb_info *mbi;

    mbi_w = ctx->mb_cols + 1; /* For left border col */
    mbi_h = ctx->mb_rows + 1; /* For above border row */

    /* TODO: Handle buffer size changes */

    if (!ctx->mb_info_storage)
        ctx->mb_info_storage = calloc(mbi_w * mbi_h,
                                      sizeof(*ctx->mb_info_storage));

    if (!ctx->mb_info_rows_storage)
        ctx->mb_info_rows_storage = calloc(mbi_h,
                                           sizeof(*ctx->mb_info_rows_storage));

    /* Set up row pointers */
    mbi = ctx->mb_info_storage + 1;

    for (i = 0; i < mbi_h; i++)
    {
        ctx->mb_info_rows_storage[i] = mbi;
        mbi += mbi_w;
    }

    ctx->mb_info_rows = ctx->mb_info_rows_storage + 1;
}
