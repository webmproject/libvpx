/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#ifndef bool_coder_h
#define bool_coder_h 1

/* Arithmetic bool coder with largish probability range.
   Timothy S Murphy  6 August 2004 */

/* So as not to force users to drag in too much of my idiosyncratic C++ world,
   I avoid fancy storage management. */

#include <assert.h>

#include <stddef.h>
#include <stdio.h>

typedef unsigned char vp8bc_index_t; // probability index

/* There are a couple of slight variants in the details of finite-precision
   arithmetic coding.  May be safely ignored by most users. */

enum vp8bc_rounding
{
    vp8bc_down = 0,     // just like VP8
    vp8bc_down_full = 1, // handles minimum probability correctly
    vp8bc_up = 2
};

#if _MSC_VER

/* Note that msvc by default does not inline _anything_ (regardless of the
   setting of inline_depth) and that a command-line option (-Ob1 or -Ob2)
   is required to inline even the smallest functions. */

#   pragma inline_depth( 255)           // I mean it when I inline something
#   pragma warning( disable : 4099)     // No class vs. struct harassment
#   pragma warning( disable : 4250)     // dominance complaints
#   pragma warning( disable : 4284)     // operator-> in templates
#   pragma warning( disable : 4800)     // bool conversion

// don't let prefix ++,-- stand in for postfix, disaster would ensue

#   pragma warning( error : 4620 4621)

#endif  // _MSC_VER


#if __cplusplus

// Sometimes one wishes to be definite about integer lengths.

struct int_types
{
    typedef const bool cbool;
    typedef const signed char cchar;
    typedef const short cshort;
    typedef const int cint;
    typedef const int clong;

    typedef const double cdouble;
    typedef const size_t csize_t;

    typedef unsigned char uchar;    // 8 bits
    typedef const uchar cuchar;

    typedef short int16;
    typedef unsigned short uint16;
    typedef const int16 cint16;
    typedef const uint16 cuint16;

    typedef int int32;
    typedef unsigned int uint32;
    typedef const int32 cint32;
    typedef const uint32 cuint32;

    typedef unsigned int uint;
    typedef unsigned int ulong;
    typedef const uint cuint;
    typedef const ulong culong;


    // All structs consume space, may as well have a vptr.

    virtual ~int_types();
};


struct bool_coder_spec;
struct bool_coder;
struct bool_writer;
struct bool_reader;


struct bool_coder_namespace : int_types
{
    typedef vp8bc_index_t Index;
    typedef bool_coder_spec Spec;
    typedef const Spec c_spec;

    enum Rounding
    {
        Down = vp8bc_down,
        down_full = vp8bc_down_full,
        Up = vp8bc_up
    };
};


// Archivable specification of a bool coder includes rounding spec
// and probability mapping table.  The latter replaces a uchar j
// (0 <= j < 256) with an arbitrary uint16 tbl[j] = p.
// p/65536 is then the probability of a zero.

struct bool_coder_spec : bool_coder_namespace
{
    friend struct bool_coder;
    friend struct bool_writer;
    friend struct bool_reader;
    friend struct bool_coder_spec_float;
    friend struct bool_coder_spec_explicit_table;
    friend struct bool_coder_spec_exponential_table;
    friend struct BPsrc;
private:
    uint w;                 // precision
    Rounding r;

    uint ebits, mbits, ebias;
    uint32 mmask;

    Index max_index, half_index;

    uint32 mantissa(Index i) const
    {
        assert(i < half_index);
        return (1 << mbits) + (i & mmask);
    }
    uint exponent(Index i) const
    {
        assert(i < half_index);
        return ebias - (i >> mbits);
    }

    uint16 Ptbl[256];       // kinda clunky, but so is storage management.

    /* Cost in bits of encoding a zero at every probability, scaled by 2^20.
       Assumes that index is at most 8 bits wide. */

    uint32 Ctbl[256];

    uint32 split(Index i, uint32 R) const    // 1 <= split <= max( 1, R-1)
    {
        if (!ebias)
            return 1 + (((R - 1) * Ptbl[i]) >> 16);

        if (i >= half_index)
            return R - split(max_index - i, R);

        return 1 + (((R - 1) * mantissa(i)) >> exponent(i));
    }

    uint32 max_range() const
    {
        return (1 << w) - (r == down_full ? 0 : 1);
    }
    uint32 min_range() const
    {
        return (1 << (w - 1)) + (r == down_full ? 1 : 0);
    }
    uint32 Rinc() const
    {
        return r == Up ? 1 : 0;
    }

    void check_prec() const;

    bool float_init(uint Ebits, uint Mbits);

    void cost_init();

    bool_coder_spec(
        uint prec, Rounding rr, uint Ebits = 0, uint Mbits = 0
    )
        : w(prec), r(rr)
    {
        float_init(Ebits, Mbits);
    }
public:
    // Read complete spec from file.
    bool_coder_spec(FILE *);

    // Write spec to file.
    void dump(FILE *) const;

    // return probability index best approximating prob.
    Index operator()(double prob) const;

    // probability corresponding to index
    double operator()(Index i) const;

    Index complement(Index i) const
    {
        return max_index - i;
    }

    Index max_index() const
    {
        return max_index;
    }
    Index half_index() const
    {
        return half_index;
    }

    uint32 cost_zero(Index i) const
    {
        return Ctbl[i];
    }
    uint32 cost_one(Index i) const
    {
        return Ctbl[ max_index - i];
    }
    uint32 cost_bit(Index i, bool b) const
    {
        return Ctbl[b? max_index-i:i];
    }
};


/* Pseudo floating-point probability specification.

   At least one of Ebits and Mbits must be nonzero.

   Since all arithmetic is done at 32 bits, Ebits is at most 5.

   Total significant bits in index is Ebits + Mbits + 1.

   Below the halfway point (i.e. when the top significant bit is 0),
   the index is (e << Mbits) + m.

   The exponent e is between 0 and (2**Ebits) - 1,
   the mantissa m is between 0 and (2**Mbits) - 1.

   Prepending an implicit 1 to the mantissa, the probability is then

        (2**Mbits + m) >> (e - 2**Ebits - 1 - Mbits),

   which has (1/2)**(2**Ebits + 1) as a minimum
   and (1/2) * [1 - 2**(Mbits + 1)] as a maximum.

   When the index is above the halfway point, the probability is the
   complement of the probability associated to the complement of the index.

   Note that the probability increases with the index and that, because of
   the symmetry, we cannot encode probability exactly 1/2; though we
   can get as close to 1/2 as we like, provided we have enough Mbits.

   The latter is of course not a problem in practice, one never has
   exact probabilities and entropy errors are second order, that is, the
   "overcoding" of a zero will be largely compensated for by the
   "undercoding" of a one (or vice-versa).

   Compared to arithmetic probability specs (a la VP8), this will do better
   at very high and low probabilities and worse at probabilities near 1/2,
   as well as facilitating the usage of wider or narrower probability indices.
*/

struct bool_coder_spec_float : bool_coder_spec
{
    bool_coder_spec_float(
        uint Ebits = 3, uint Mbits = 4, Rounding rr = down_full, uint prec = 12
    )
        : bool_coder_spec(prec, rr, Ebits, Mbits)
    {
        cost_init();
    }
};


struct bool_coder_spec_explicit_table : bool_coder_spec
{
    bool_coder_spec_explicit_table(
        cuint16 probability_table[256] = 0,  // default is tbl[i] = i << 8.
        Rounding = down_full,
        uint precision = 16
    );
};

// Contruct table via multiplicative interpolation between
// p[128] = 1/2  and p[0] = (1/2)^x.
// Since we are working with 16-bit precision, x is at most 16.
// For probabilities to increase with i, we must have x > 1.
// For 0 <= i <= 128, p[i] = (1/2)^{ 1 + [1 - (i/128)]*[x-1] }.
// Finally, p[128+i] = 1 - p[128 - i].

struct bool_coder_spec_exponential_table : bool_coder_spec
{
    bool_coder_spec_exponential_table(uint x, Rounding = down_full, uint prec = 16);
};


// Commonalities between writer and reader.

struct bool_coder : bool_coder_namespace
{
    friend struct bool_writer;
    friend struct bool_reader;
    friend struct BPsrc;
private:
    uint32 Low, Range;
    cuint32 min_range;
    cuint32 rinc;
    c_spec spec;

    void _reset()
    {
        Low = 0;
        Range = spec.max_range();
    }

    bool_coder(c_spec &s)
        :  min_range(s.min_range()),
           rinc(s.Rinc()),
           spec(s)
    {
        _reset();
    }

    uint32 half() const
    {
        return 1 + ((Range - 1) >> 1);
    }
public:
    c_spec &Spec() const
    {
        return spec;
    }
};


struct bool_writer : bool_coder
{
    friend struct BPsrc;
private:
    uchar *Bstart, *Bend, *B;
    int bit_lag;
    bool is_toast;
    void carry();
    void reset()
    {
        _reset();
        bit_lag = 32 - spec.w;
        is_toast = 0;
    }
    void raw(bool value, uint32 split);
public:
    bool_writer(c_spec &, uchar *Dest, size_t Len);
    virtual ~bool_writer();

    void operator()(Index p, bool v)
    {
        raw(v, spec.split(p, Range));
    }

    uchar *buf() const
    {
        return Bstart;
    }
    size_t bytes_written() const
    {
        return B - Bstart;
    }

    // Call when done with input, flushes internal state.
    // DO NOT write any more data after calling this.

    bool_writer &flush();

    void write_bits(int n, uint val)
    {
        if (n)
        {
            uint m = 1 << (n - 1);

            do
            {
                raw((bool)(val & m), half());
            }
            while (m >>= 1);
        }
    }

#   if 0
    // We are agnostic about storage management.
    // By default, overflows throw an assert but user can
    // override to provide an expanding buffer using ...

    virtual void overflow(uint Len) const;

    // ... this function copies already-written data into new buffer
    // and retains new buffer location.

    void new_buffer(uchar *dest, uint Len);

    // Note that storage management is the user's responsibility.
#   endif
};


// This could be adjusted to use a little less lookahead.

struct bool_reader : bool_coder
{
    friend struct BPsrc;
private:
    cuchar *const Bstart;   // for debugging
    cuchar *B;
    cuchar *const Bend;
    cuint shf;
    uint bct;
    bool raw(uint32 split);
public:
    bool_reader(c_spec &s, cuchar *src, size_t Len);

    bool operator()(Index p)
    {
        return raw(spec.split(p, Range));
    }

    uint read_bits(int num_bits)
    {
        uint v = 0;

        while (--num_bits >= 0)
            v += v + (raw(half()) ? 1 : 0);

        return v;
    }
};

extern "C" {

#endif /*  __cplusplus */


    /* C interface */

    typedef struct bool_coder_spec bool_coder_spec;
    typedef struct bool_writer bool_writer;
    typedef struct bool_reader bool_reader;

    typedef const bool_coder_spec c_bool_coder_spec;
    typedef const bool_writer c_bool_writer;
    typedef const bool_reader c_bool_reader;


    /* Optionally override default precision when constructing coder_specs.
       Just pass a zero pointer if you don't care.
       Precision is at most 16 bits for table specs, at most 23 otherwise. */

    struct vp8bc_prec
    {
        enum vp8bc_rounding r;      /* see top header file for def */
        unsigned int prec;          /* range precision in bits */
    };

    typedef const struct vp8bc_prec vp8bc_c_prec;

    /* bool_coder_spec contains mapping of uchars to actual probabilities
       (16 bit uints) as well as (usually immaterial) selection of
       exact finite-precision algorithm used (for now, the latter can only
       be overridden using the C++ interface).
       See comments above the corresponding C++ constructors for discussion,
       especially of exponential probability table generation. */

    bool_coder_spec *vp8bc_vp8spec(); // just like vp8

    bool_coder_spec *vp8bc_literal_spec(
        const unsigned short prob_map[256],  // 0 is like vp8 w/more precision
        vp8bc_c_prec*
    );

    bool_coder_spec *vp8bc_float_spec(
        unsigned int exponent_bits, unsigned int mantissa_bits, vp8bc_c_prec*
    );

    bool_coder_spec *vp8bc_exponential_spec(unsigned int min_exp, vp8bc_c_prec *);

    bool_coder_spec *vp8bc_spec_from_file(FILE *);


    void vp8bc_destroy_spec(c_bool_coder_spec *);

    void vp8bc_spec_to_file(c_bool_coder_spec *, FILE *);


    /* Nearest index to supplied probability of zero, 0 <= prob <= 1. */

    vp8bc_index_t vp8bc_index(c_bool_coder_spec *, double prob);

    vp8bc_index_t vp8bc_index_from_counts(
        c_bool_coder_spec *p, unsigned int zero_ct, unsigned int one_ct
    );

    /* In case you want to look */

    double vp8bc_probability(c_bool_coder_spec *, vp8bc_index_t);

    /* Opposite index */

    vp8bc_index_t vp8bc_complement(c_bool_coder_spec *, vp8bc_index_t);

    /* Cost in bits of encoding a zero at given probability, scaled by 2^20.
       (assumes that an int holds at least 32 bits). */

    unsigned int vp8bc_cost_zero(c_bool_coder_spec *, vp8bc_index_t);

    unsigned int vp8bc_cost_one(c_bool_coder_spec *, vp8bc_index_t);
    unsigned int vp8bc_cost_bit(c_bool_coder_spec *, vp8bc_index_t, int);


    /* bool_writer interface */

    /* Length = 0 disables checking for writes beyond buffer end. */

    bool_writer *vp8bc_create_writer(
        c_bool_coder_spec *, unsigned char *Destination, size_t Length
    );

    /* Flushes out any buffered data and returns total # of bytes written. */

    size_t vp8bc_destroy_writer(bool_writer *);

    void vp8bc_write_bool(bool_writer *, int boolean_val, vp8bc_index_t false_prob);

    void vp8bc_write_bits(
        bool_writer *, unsigned int integer_value, int number_of_bits
    );

    c_bool_coder_spec *vp8bc_writer_spec(c_bool_writer *);


    /* bool_reader interface */

    /* Length = 0 disables checking for reads beyond buffer end. */

    bool_reader *vp8bc_create_reader(
        c_bool_coder_spec *, const unsigned char *Source, size_t Length
    );
    void vp8bc_destroy_reader(bool_reader *);

    int vp8bc_read_bool(bool_reader *, vp8bc_index_t false_prob);

    unsigned int vp8bc_read_bits(bool_reader *, int number_of_bits);

    c_bool_coder_spec *vp8bc_reader_spec(c_bool_reader *);

#if __cplusplus
}
#endif

#endif  /* bool_coder_h */
