/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_RDOPT_H
#define __INC_RDOPT_H

#define RDCOST(RM,DM,R,D) ( ((128+(R)*(RM)) >> 8) + (DM)*(D) )

static void insertsortmv(int arr[], int len)
{
    int i, j, k;

    for ( i = 1 ; i <= len-1 ; i++ )
    {
        for ( j = 0 ; j < i ; j++ )
        {
            if ( arr[j] > arr[i] )
            {
                int temp;

                temp = arr[i];

                for ( k = i; k >j; k--)
                    arr[k] = arr[k - 1] ;

                arr[j] = temp ;
            }
        }
    }
}

static void insertsortsad(int arr[],int idx[], int len)
{
    int i, j, k;

    for ( i = 1 ; i <= len-1 ; i++ )
    {
        for ( j = 0 ; j < i ; j++ )
        {
            if ( arr[j] > arr[i] )
            {
                int temp, tempi;

                temp = arr[i];
                tempi = idx[i];

                for ( k = i; k >j; k--)
                {
                    arr[k] = arr[k - 1] ;
                    idx[k] = idx[k - 1];
                }

                arr[j] = temp ;
                idx[j] = tempi;
            }
        }
    }
}

extern void vp8_initialize_rd_consts(VP8_COMP *cpi, int Qvalue);
extern void vp8_rd_pick_inter_mode(VP8_COMP *cpi, MACROBLOCK *x, int recon_yoffset, int recon_uvoffset, int *returnrate, int *returndistortion, int *returnintra);
extern void vp8_rd_pick_intra_mode(VP8_COMP *cpi, MACROBLOCK *x, int *rate);

extern void vp8_mv_pred
(
    VP8_COMP *cpi,
    MACROBLOCKD *xd,
    const MODE_INFO *here,
    int_mv *mvp,
    int refframe,
    int *ref_frame_sign_bias,
    int *sr,
    int near_sadidx[]
);
void vp8_cal_sad(VP8_COMP *cpi, MACROBLOCKD *xd, MACROBLOCK *x, int recon_yoffset, int near_sadidx[]);

#endif
