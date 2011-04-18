#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
#pragma OPENCL EXTENSION cl_amd_printf : enable

typedef unsigned char uc;
typedef signed char sc;

__inline signed char vp8_filter_mask(sc, sc, uc, uc, uc, uc, uc, uc, uc, uc);
__inline signed char vp8_simple_filter_mask(signed char, signed char, uc, uc, uc, uc);
__inline signed char vp8_hevmask(signed char, uc, uc, uc, uc);
__inline signed char vp8_signed_char_clamp(int);

__inline void vp8_mbfilter(signed char mask,signed char hev,global uc *op2,
    global uc *op1,global uc *op0,global uc *oq0,global uc *oq1,global uc *oq2);

void vp8_simple_filter(signed char mask,global uc *base, int op1_off,int op0_off,int oq0_off,int oq1_off);


typedef struct
{
    signed char lim[16];
    signed char flim[16];
    signed char thr[16];
    signed char mbflim[16];
    signed char mbthr[16];
    signed char uvlim[16];
    signed char uvflim[16];
    signed char uvthr[16];
    signed char uvmbflim[16];
    signed char uvmbthr[16];
} loop_filter_info;




void vp8_filter(
    signed char mask,
    signed char hev,
    global uc *base,
    int op1_off,
    int op0_off,
    int oq0_off,
    int oq1_off
)
{

    global uc *op1 = &base[op1_off];
    global uc *op0 = &base[op0_off];
    global uc *oq0 = &base[oq0_off];
    global uc *oq1 = &base[oq1_off];

    signed char ps0, qs0;
    signed char ps1, qs1;
    signed char vp8_filter, Filter1, Filter2;
    signed char u;

    ps1 = (signed char) * op1 ^ 0x80;
    ps0 = (signed char) * op0 ^ 0x80;
    qs0 = (signed char) * oq0 ^ 0x80;
    qs1 = (signed char) * oq1 ^ 0x80;

    /* add outer taps if we have high edge variance */
    vp8_filter = vp8_signed_char_clamp(ps1 - qs1);
    vp8_filter &= hev;

    /* inner taps */
    vp8_filter = vp8_signed_char_clamp(vp8_filter + 3 * (qs0 - ps0));
    vp8_filter &= mask;

    /* save bottom 3 bits so that we round one side +4 and the other +3
     * if it equals 4 we'll set to adjust by -1 to account for the fact
     * we'd round 3 the other way
     */
    Filter1 = vp8_signed_char_clamp(vp8_filter + 4);
    Filter2 = vp8_signed_char_clamp(vp8_filter + 3);
    Filter1 >>= 3;
    Filter2 >>= 3;
    u = vp8_signed_char_clamp(qs0 - Filter1);
    *oq0 = u ^ 0x80;
    u = vp8_signed_char_clamp(ps0 + Filter2);
    *op0 = u ^ 0x80;
    vp8_filter = Filter1;

    /* outer tap adjustments */
    vp8_filter += 1;
    vp8_filter >>= 1;
    vp8_filter &= ~hev;

    u = vp8_signed_char_clamp(qs1 - vp8_filter);
    *oq1 = u ^ 0x80;
    u = vp8_signed_char_clamp(ps1 + vp8_filter);
    *op1 = u ^ 0x80;
}


kernel void vp8_loop_filter_horizontal_edge_kernel
(
    global unsigned char *s_base,
    int s_off,
    int p, /* pitch */
    global signed char *flimit,
    global signed char *limit,
    global signed char *thresh,
    int off_stride
)
{
    int  hev = 0; /* high edge variance */
    signed char mask = 0;
    int i = get_global_id(0);

    if (i < get_global_size(0)){
        s_off += i;

        mask = vp8_filter_mask(limit[i], flimit[i], s_base[s_off - 4*p],
                s_base[s_off - 3*p], s_base[s_off - 2*p], s_base[s_off - p],
                s_base[s_off], s_base[s_off + p], s_base[s_off + 2*p],
                s_base[s_off + 3*p]);

        hev = vp8_hevmask(thresh[i], s_base[s_off - 2*p], s_base[s_off - p],
                s_base[s_off], s_base[s_off+p]);

        vp8_filter(mask, hev, s_base, s_off - 2 * p, s_off - p, s_off,
                s_off + p);
    }
}


kernel void vp8_loop_filter_vertical_edge_kernel
(
    global unsigned char *s_base,
    int s_off,
    int p,
    global signed char *flimit,
    global signed char *limit,
    global signed char *thresh,
    int off_stride
)
{

    int  hev = 0; /* high edge variance */
    signed char mask = 0;
    int i = get_global_id(0);

    if ( i < get_global_size(0) ){
        s_off += p * i;
        mask = vp8_filter_mask(limit[i], flimit[i],
                s_base[s_off-4], s_base[s_off-3], s_base[s_off-2],
                s_base[s_off-1], s_base[s_off], s_base[s_off+1],
                s_base[s_off+2], s_base[s_off+3]);

        hev = vp8_hevmask(thresh[i], s_base[s_off-2], s_base[s_off-1],
                s_base[s_off], s_base[s_off+1]);

        vp8_filter(mask, hev, s_base, s_off - 2, s_off - 1, s_off, s_off + 1);

    }
}


kernel void vp8_mbloop_filter_horizontal_edge_kernel
(
    global unsigned char *s_base,
    int s_off,
    int p,
    global signed char *flimit,
    global signed char *limit,
    global signed char *thresh,
    int off_stride
)
{

    global uc *s = s_base+s_off;

    signed char hev = 0; /* high edge variance */
    signed char mask = 0;
    int i = get_global_id(0);

    if (i < get_global_size(0)){
        s += i;

        mask = vp8_filter_mask(limit[i], flimit[i],
                               s[-4*p], s[-3*p], s[-2*p], s[-1*p],
                               s[0*p], s[1*p], s[2*p], s[3*p]);

        hev = vp8_hevmask(thresh[i], s[-2*p], s[-1*p], s[0*p], s[1*p]);

        vp8_mbfilter(mask, hev, s - 3 * p, s - 2 * p, s - 1 * p, s, s + 1 * p, s + 2 * p);

    }
}


kernel void vp8_mbloop_filter_vertical_edge_kernel
(
    global unsigned char *s_base,
    int s_off,
    int p,
    global signed char *flimit,
    global signed char *limit,
    global signed char *thresh,
    int off_stride
)
{

    global uc *s = s_base + s_off;

    signed char hev = 0; /* high edge variance */
    signed char mask = 0;
    int i = get_global_id(0);

    if (i < get_global_size(0)){
        s += p * i;

        mask = vp8_filter_mask(limit[i], flimit[i],
                               s[-4], s[-3], s[-2], s[-1], s[0], s[1], s[2], s[3]);

        hev = vp8_hevmask(thresh[i], s[-2], s[-1], s[0], s[1]);

        vp8_mbfilter(mask, hev, s - 3, s - 2, s - 1, s, s + 1, s + 2);

    }
}


kernel void vp8_loop_filter_simple_horizontal_edge_kernel
(
    global unsigned char *s_base,
    int s_off,
    int p,
    global const signed char *flimit,
    global const signed char *limit,
    global const signed char *thresh,
    int off_stride
)
{

    signed char mask = 0;
    int i = get_global_id(0);
    (void) thresh;

    if (i < get_global_size(0))
    {
        s_off += i;
        mask = vp8_simple_filter_mask(limit[i], flimit[i], s_base[s_off-2*p], s_base[s_off-p], s_base[s_off], s_base[s_off+p]);
        vp8_simple_filter(mask, s_base, s_off - 2 * p, s_off - 1 * p, s_off, s_off + 1 * p);
    }
}


kernel void vp8_loop_filter_simple_vertical_edge_kernel
(
    global unsigned char *s_base,
    int s_off,
    int p,
    global signed char *flimit,
    global signed char *limit,
    global signed char *thresh,
    int off_stride
)
{

    signed char mask = 0;
    int i = get_global_id(0);
    (void) thresh;

    if (i < get_global_size(0)){
        s_off += p * i;
        mask = vp8_simple_filter_mask(limit[i], flimit[i], s_base[s_off-2], s_base[s_off-1], s_base[s_off], s_base[s_off+1]);
        vp8_simple_filter(mask, s_base, s_off - 2, s_off - 1, s_off, s_off + 1);
    }

}



//Inline and non-kernel functions follow.

__inline void vp8_mbfilter(
    signed char mask,
    signed char hev,
    global uc *op2,
    global uc *op1,
    global uc *op0,
    global uc *oq0,
    global uc *oq1,
    global uc *oq2
)
{
    signed char s, u;
    signed char vp8_filter, Filter1, Filter2;
    signed char ps2 = (signed char) * op2 ^ 0x80;
    signed char ps1 = (signed char) * op1 ^ 0x80;
    signed char ps0 = (signed char) * op0 ^ 0x80;
    signed char qs0 = (signed char) * oq0 ^ 0x80;
    signed char qs1 = (signed char) * oq1 ^ 0x80;
    signed char qs2 = (signed char) * oq2 ^ 0x80;

    /* add outer taps if we have high edge variance */
    vp8_filter = vp8_signed_char_clamp(ps1 - qs1);
    vp8_filter = vp8_signed_char_clamp(vp8_filter + 3 * (qs0 - ps0));
    vp8_filter &= mask;

    Filter2 = vp8_filter;
    Filter2 &= hev;

    /* save bottom 3 bits so that we round one side +4 and the other +3 */
    Filter1 = vp8_signed_char_clamp(Filter2 + 4);
    Filter2 = vp8_signed_char_clamp(Filter2 + 3);
    Filter1 >>= 3;
    Filter2 >>= 3;
    qs0 = vp8_signed_char_clamp(qs0 - Filter1);
    ps0 = vp8_signed_char_clamp(ps0 + Filter2);


    /* only apply wider filter if not high edge variance */
    vp8_filter &= ~hev;
    Filter2 = vp8_filter;

    /* roughly 3/7th difference across boundary */
    u = vp8_signed_char_clamp((63 + Filter2 * 27) >> 7);
    s = vp8_signed_char_clamp(qs0 - u);
    *oq0 = s ^ 0x80;
    s = vp8_signed_char_clamp(ps0 + u);
    *op0 = s ^ 0x80;

    /* roughly 2/7th difference across boundary */
    u = vp8_signed_char_clamp((63 + Filter2 * 18) >> 7);
    s = vp8_signed_char_clamp(qs1 - u);
    *oq1 = s ^ 0x80;
    s = vp8_signed_char_clamp(ps1 + u);
    *op1 = s ^ 0x80;

    /* roughly 1/7th difference across boundary */
    u = vp8_signed_char_clamp((63 + Filter2 * 9) >> 7);
    s = vp8_signed_char_clamp(qs2 - u);
    *oq2 = s ^ 0x80;
    s = vp8_signed_char_clamp(ps2 + u);
    *op2 = s ^ 0x80;
}


__inline signed char vp8_signed_char_clamp(int t)
{
    t = (t < -128 ? -128 : t);
    t = (t > 127 ? 127 : t);
    return (signed char) t;
}


/* is there high variance internal edge ( 11111111 yes, 00000000 no) */
__inline signed char vp8_hevmask(signed char thresh, uc p1, uc p0, uc q0, uc q1)
{
    signed char hev = 0;
    hev  |= (abs(p1 - p0) > thresh) * -1;
    hev  |= (abs(q1 - q0) > thresh) * -1;
    return hev;
}


/* should we apply any filter at all ( 11111111 yes, 00000000 no) */
__inline signed char vp8_filter_mask(
    signed char limit,
    signed char flimit,
     uc p3, uc p2, uc p1, uc p0, uc q0, uc q1, uc q2, uc q3)
{
    signed char mask = 0;
    mask |= (abs(p3 - p2) > limit) * -1;
    mask |= (abs(p2 - p1) > limit) * -1;
    mask |= (abs(p1 - p0) > limit) * -1;
    mask |= (abs(q1 - q0) > limit) * -1;
    mask |= (abs(q2 - q1) > limit) * -1;
    mask |= (abs(q3 - q2) > limit) * -1;
    mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2  > flimit * 2 + limit) * -1;
    mask = ~mask;
    return mask;
}

/* should we apply any filter at all ( 11111111 yes, 00000000 no) */
__inline signed char vp8_simple_filter_mask(
    signed char limit,
    signed char flimit,
    uc p1,
    uc p0,
    uc q0,
    uc q1
)
{
    signed char mask = (abs(p0 - q0) * 2 + abs(p1 - q1) / 2  <= flimit * 2 + limit) * -1;
    return mask;
}

void vp8_simple_filter(
    signed char mask,
    global uc *base,
    int op1_off,
    int op0_off,
    int oq0_off,
    int oq1_off
)
{

    global uc *op1 = base + op1_off;
    global uc *op0 = base + op0_off;
    global uc *oq0 = base + oq0_off;
    global uc *oq1 = base + oq1_off;

    signed char vp8_filter, Filter1, Filter2;
    signed char p1 = (signed char) * op1 ^ 0x80;
    signed char p0 = (signed char) * op0 ^ 0x80;
    signed char q0 = (signed char) * oq0 ^ 0x80;
    signed char q1 = (signed char) * oq1 ^ 0x80;
    signed char u;

    vp8_filter = vp8_signed_char_clamp(p1 - q1);
    vp8_filter = vp8_signed_char_clamp(vp8_filter + 3 * (q0 - p0));
    vp8_filter &= mask;

    /* save bottom 3 bits so that we round one side +4 and the other +3 */
    Filter1 = vp8_signed_char_clamp(vp8_filter + 4);
    Filter1 >>= 3;
    u = vp8_signed_char_clamp(q0 - Filter1);
    *oq0  = u ^ 0x80;

    Filter2 = vp8_signed_char_clamp(vp8_filter + 3);
    Filter2 >>= 3;
    u = vp8_signed_char_clamp(p0 + Filter2);
    *op0 = u ^ 0x80;
}