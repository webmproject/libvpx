/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */



/* Arithmetic bool coder with largish probability range.
   Timothy S Murphy  6 August 2004 */

#include <assert.h>
#include <math.h>

#include "bool_coder.h"

#if tim_vp8
    extern "C" {
#       include "VP8cx/treewriter.h"
    }
#endif

int_types::~int_types() {}

void bool_coder_spec::check_prec() const {
    assert( w  &&  (r==Up || w > 1)  &&  w < 24  &&  (ebias || w < 17));
}

bool bool_coder_spec::float_init( uint Ebits, uint Mbits) {
    uint b = (ebits = Ebits) + (mbits = Mbits);
    if( b) {
        assert( ebits < 6  &&  w + mbits < 31);
        assert( ebits + mbits  <  sizeof(Index) * 8);
        ebias = (1 << ebits) + 1 + mbits;
        mmask = (1 << mbits) - 1;
        max_index = ( ( half_index = 1 << b ) << 1) - 1;
    } else {
        ebias = 0;
        max_index = 255;
        half_index = 128;
    }
    check_prec();
    return b? 1:0;
}

void bool_coder_spec::cost_init()
{
    static cdouble c = -(1 << 20)/log( 2.);

    FILE *f = fopen( "costs.txt", "w");
    assert( f);

    assert( sizeof(int) >= 4);  /* for C interface */
    assert( max_index <= 255);   /* size of Ctbl */
    uint i = 0;  do {
        cdouble p = ( *this)( (Index) i);
        Ctbl[i] = (uint32) ( log( p) * c);
        fprintf(
            f, "cost( %d -> %10.7f) = %10d = %12.5f bits\n",
            i, p, Ctbl[i], (double) Ctbl[i] / (1<<20)
        );
    } while( ++i <= max_index);
    fclose( f);
}

bool_coder_spec_explicit_table::bool_coder_spec_explicit_table(
    cuint16 tbl[256], Rounding rr, uint prec
)
  : bool_coder_spec( prec, rr)
{
    check_prec();
    uint i = 0;
    if( tbl)
        do { Ptbl[i] = tbl[i];}  while( ++i < 256);
    else
        do { Ptbl[i] = i << 8;}  while( ++i < 256);
    cost_init();
}


bool_coder_spec_exponential_table::bool_coder_spec_exponential_table(
    uint x, Rounding rr, uint prec
)
  : bool_coder_spec( prec, rr)
{
    assert( x > 1  &&  x <= 16);
    check_prec();
    Ptbl[128] = 32768u;
    Ptbl[0] = (uint16) pow( 2., 16. - x);
    --x;
    int i=1;  do {
        cdouble d = pow( .5, 1. + (1. - i/128.)*x) * 65536.;
        uint16 v = (uint16) d;
        if( v < i)
            v = i;
        Ptbl[256-i] = (uint16) ( 65536U - (Ptbl[i] = v));
    } while( ++i < 128);
    cost_init();
}

bool_coder_spec::bool_coder_spec( FILE *fp) {
    fscanf( fp, "%d", &w);
    int v;
    fscanf( fp, "%d", &v);
    assert( 0 <= v  &&  v <= 2);
    r = (Rounding) v;
    fscanf( fp, "%d", &ebits);
    fscanf( fp, "%d", &mbits);
    if( float_init( ebits, mbits))
        return;
    int i=0;  do {
        uint v;
        fscanf( fp, "%d", &v);
        assert( 0 <=v  &&  v <= 65535U);
        Ptbl[i] = v;
    } while( ++i < 256);
    cost_init();
}

void bool_coder_spec::dump( FILE *fp) const {
    fprintf( fp, "%d %d %d %d\n", w, (int) r, ebits, mbits);
    if( ebits  ||  mbits)
        return;
    int i=0;  do { fprintf( fp, "%d\n", Ptbl[i]);}  while( ++i < 256);
}

vp8bc_index_t bool_coder_spec::operator()( double p) const
{
    if( p <= 0.)
        return 0;
    if( p >= 1.)
        return max_index;
    if( ebias) {
        if( p > .5)
            return max_index - ( *this)( 1. - p);
        int e;
        uint m = (uint) ldexp( frexp( p, &e), mbits + 2);
        uint x = 1 << (mbits + 1);
        assert( x <= m  &&  m < x<<1);
        if( (m = (m >> 1) + (m & 1)) >= x) {
            m = x >> 1;
            ++e;
        }
        int y = 1 << ebits;
        if( (e += y) >= y)
            return half_index - 1;
        if( e < 0)
            return 0;
        return (Index) ( (e << mbits) + (m & mmask));
    }

    cuint16 v = (uint16) (p * 65536.);
    int i = 128;
    int j = 128;
    uint16 w;
    while( w = Ptbl[i], j >>= 1) {
        if( w < v)
            i += j;
        else if( w == v)
            return (uchar) i;
        else
            i -= j;
    }
    if( w > v) {
        cuint16 x = Ptbl[i-1];
        if( v <= x  ||  w - v > v - x)
            --i;
    } else if( w < v  &&  i < 255) {
        cuint16 x = Ptbl[i+1];
        if( x <= v  ||  x - v < v - w)
            ++i;
    }
    return (Index) i;
}

double bool_coder_spec::operator()( Index i) const {
    if( !ebias)
        return Ptbl[i]/65536.;
    if( i >= half_index)
        return 1. - ( *this)( (Index) (max_index - i));
    return ldexp( (double)mantissa( i), - (int) exponent( i));
}



void bool_writer::carry() {
    uchar *p = B;
    assert( p > Bstart);
    while( *--p == 255) { assert( p > Bstart);  *p = 0;}
    ++*p;
}


bool_writer::bool_writer( c_spec& s, uchar *Dest, size_t Len)
  : bool_coder( s),
    Bstart( Dest),
    Bend( Len? Dest+Len : 0),
    B( Dest)
{
    assert( Dest);
    reset();
}

bool_writer::~bool_writer() { flush();}

#if 1
    extern "C" { int bc_v = 0;}
#else
#   define bc_v 0
#endif


void bool_writer::raw( bool value, uint32 s) {
    uint32 L = Low;

    assert( Range >= min_range  &&  Range <= spec.max_range());
    assert( !is_toast  &&  s  &&  s < Range);

    if( bc_v) printf(
        "Writing a %d, B %x  Low %x  Range %x  s %x   blag %d ...\n",
        value? 1:0, B-Bstart, Low, Range, s, bit_lag
    );
    if( value) {
        L += s;
        s = Range - s;
    } else
        s -= rinc;
    if( s < min_range) {
        int ct = bit_lag;  do {
            if( !--ct) {
                ct = 8;
                if( L & (1 << 31))
                    carry();
                assert( !Bend  ||  B < Bend);
                *B++ = (uchar) (L >> 23);
                L &= (1<<23) - 1;
            }
        } while( L += L, (s += s + rinc) < min_range);
        bit_lag = ct;
    }
    Low = L;
    Range = s;
    if( bc_v)
        printf(
            "...done, B %x  Low %x  Range %x  blag %d \n",
                B-Bstart, Low, Range, bit_lag
        );
}

bool_writer& bool_writer::flush() {
    if( is_toast)
        return *this;
    int b = bit_lag;
    uint32 L = Low;
    assert( b);
    if( L & (1 << (32 - b)))
        carry();
    L <<= b & 7;
    b >>= 3;
    while( --b >= 0)
        L <<= 8;
    b = 4;
    assert( !Bend  ||  B + 4 <= Bend);
    do {
        *B++ = (uchar) (L >> 24);
        L <<= 8;
    } while( --b);
    is_toast = 1;
    return *this;
}


bool_reader::bool_reader( c_spec& s, cuchar *src, size_t Len)
  : bool_coder( s),
    Bstart( src),
    B( src),
    Bend( Len? src+Len : 0),
    shf( 32 - s.w),
    bct( 8)
{
    int i = 4;  do { Low <<= 8;  Low |= *B++;}  while( --i);
}


bool bool_reader::raw( uint32 s) {

    bool val = 0;
    uint32 L = Low;
    cuint32 S = s << shf;

    assert( Range >= min_range  &&  Range <= spec.max_range());
    assert( s  &&  s < Range  &&  (L >> shf) < Range);

    if( bc_v)
        printf(
            "Reading, B %x  Low %x  Range %x  s %x  bct %d ...\n",
            B-Bstart, Low, Range, s, bct
        );

    if( L >= S) {
        L -= S;
        s = Range - s;
        assert( L < (s << shf));
        val = 1;
    } else
        s -= rinc;
    if( s < min_range) {
        int ct = bct;
        do {
            assert( ~L & (1 << 31));
            L += L;
            if( !--ct) {
                ct = 8;
                if( !Bend  ||  B < Bend)
                    L |= *B++;
            }
        } while( (s += s + rinc) < min_range);
        bct = ct;
    }
    Low = L;
    Range = s;
    if( bc_v)
        printf(
            "...done, val %d  B %x  Low %x  Range %x  bct %d\n",
            val? 1:0, B-Bstart, Low, Range, bct
        );
    return val;
}


/* C interfaces */

// spec interface

struct NS : bool_coder_namespace {
    static Rounding r( vp8bc_c_prec *p, Rounding rr =down_full) {
        return p? (Rounding) p->r : rr;
    }
};

bool_coder_spec *vp8bc_vp6spec() {
    return new bool_coder_spec_explicit_table( 0, bool_coder_namespace::Down, 8);
}
bool_coder_spec *vp8bc_float_spec(
    unsigned int Ebits, unsigned int Mbits, vp8bc_c_prec *p
) {
    return new bool_coder_spec_float( Ebits, Mbits, NS::r( p), p? p->prec : 12);
}
bool_coder_spec *vp8bc_literal_spec(
    const unsigned short m[256], vp8bc_c_prec *p
) {
    return new bool_coder_spec_explicit_table( m, NS::r( p), p? p->prec : 16);
}
bool_coder_spec *vp8bc_exponential_spec( unsigned int x, vp8bc_c_prec *p)
{
    return new bool_coder_spec_exponential_table( x, NS::r( p), p? p->prec : 16);
}
bool_coder_spec *vp8bc_spec_from_file( FILE *fp) {
    return new bool_coder_spec( fp);
}
void vp8bc_destroy_spec( c_bool_coder_spec *p) { delete p;}

void vp8bc_spec_to_file( c_bool_coder_spec *p, FILE *fp) { p->dump( fp);}

vp8bc_index_t vp8bc_index( c_bool_coder_spec *p, double x) {
    return ( *p)( x);
}

vp8bc_index_t vp8bc_index_from_counts(
    c_bool_coder_spec *p, unsigned int L, unsigned int R
) {
    return ( *p)( (R += L)? (double) L/R : .5);
}

double vp8bc_probability( c_bool_coder_spec *p, vp8bc_index_t i) {
    return ( *p)( i);
}

vp8bc_index_t vp8bc_complement( c_bool_coder_spec *p, vp8bc_index_t i) {
    return p->complement( i);
}
unsigned int vp8bc_cost_zero( c_bool_coder_spec *p, vp8bc_index_t i) {
    return p->cost_zero( i);
}
unsigned int vp8bc_cost_one( c_bool_coder_spec *p, vp8bc_index_t i) {
    return p->cost_one( i);
}
unsigned int vp8bc_cost_bit( c_bool_coder_spec *p, vp8bc_index_t i, int v) {
    return p->cost_bit( i, v);
}

#if tim_vp8
    extern "C" int tok_verbose;

#   define dbg_l 1000000

    static vp8bc_index_t dbg_i [dbg_l];
    static char dbg_v [dbg_l];
    static size_t dbg_w = 0, dbg_r = 0;
#endif

// writer interface

bool_writer *vp8bc_create_writer(
    c_bool_coder_spec *p, unsigned char *D, size_t L
) {
    return new bool_writer( *p, D, L);
}

size_t vp8bc_destroy_writer( bool_writer *p) {
    const size_t s = p->flush().bytes_written();
    delete p;
    return s;
}

void vp8bc_write_bool( bool_writer *p, int v, vp8bc_index_t i)
{
#   if tim_vp8
        // bc_v = dbg_w < 10;
        if( bc_v = tok_verbose)
            printf( " writing %d at prob %d\n", v? 1:0, i);
        accum_entropy_bc( &p->Spec(), i, v);

        ( *p)( i, (bool) v);

        if( dbg_w < dbg_l) {
            dbg_i [dbg_w] = i;
            dbg_v [dbg_w++] = v? 1:0;
        }
#   else
        ( *p)( i, (bool) v);
#   endif
}

void vp8bc_write_bits( bool_writer *p, unsigned int v, int n)
{
#   if tim_vp8
        {
            c_bool_coder_spec * const s = & p->Spec();
            const vp8bc_index_t i = s->half_index();
            int m = n;
            while( --m >= 0)
                accum_entropy_bc( s, i, (v>>m) & 1);
        }
#   endif

    p->write_bits( n, v);
}

c_bool_coder_spec *vp8bc_writer_spec( c_bool_writer *w) { return & w->Spec();}

// reader interface

bool_reader *vp8bc_create_reader(
    c_bool_coder_spec *p, const unsigned char *S, size_t L
) {
    return new bool_reader( *p, S, L);
}

void vp8bc_destroy_reader( bool_reader * p) { delete p;}

int vp8bc_read_bool( bool_reader *p, vp8bc_index_t i)
{
#   if tim_vp8
        // bc_v = dbg_r < 10;
        bc_v = tok_verbose;
        const int v = ( *p)( i)? 1:0;
        if( tok_verbose)
            printf( " reading %d at prob %d\n", v, i);
        if( dbg_r < dbg_l) {
            assert( dbg_r <= dbg_w);
            if( i != dbg_i[dbg_r]  ||  v != dbg_v[dbg_r]) {
                printf(
        "Position %d: INCORRECTLY READING %d  prob %d, wrote %d  prob %d\n",
                    dbg_r, v, i, dbg_v[dbg_r], dbg_i[dbg_r]
                );
            }
            ++dbg_r;
        }
        return v;
#   else
        return ( *p)( i)? 1:0;
#   endif
}

unsigned int vp8bc_read_bits( bool_reader *p, int n) { return p->read_bits( n);}

c_bool_coder_spec *vp8bc_reader_spec( c_bool_reader *r) { return & r->Spec();}

#undef bc_v
