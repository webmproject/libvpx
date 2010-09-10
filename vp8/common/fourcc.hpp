/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef FOURCC_HPP
#define FOURCC_HPP

#include <iosfwd>
#include <cstring>


#if defined(__POWERPC__) || defined(__APPLE__) || defined(__MERKS__)
using namespace std;
#endif

class four_cc
{
public:

    four_cc();
    four_cc(const char*);
    explicit four_cc(unsigned long);

    bool operator==(const four_cc&) const;
    bool operator!=(const four_cc&) const;

    bool operator==(const char*) const;
    bool operator!=(const char*) const;

    operator unsigned long() const;
    unsigned long as_long() const;

    four_cc& operator=(unsigned long);

    char operator[](int) const;

    std::ostream& put(std::ostream&) const;

    bool printable() const;

private:

    union
    {
        char code[4];
        unsigned long code_as_long;
    };

};


inline four_cc::four_cc()
{
}

inline four_cc::four_cc(unsigned long x)
    : code_as_long(x)
{
}

inline four_cc::four_cc(const char* str)
{
    memcpy(code, str, 4);
}


inline bool four_cc::operator==(const four_cc& rhs) const
{
    return code_as_long == rhs.code_as_long;
}

inline bool four_cc::operator!=(const four_cc& rhs) const
{
    return !operator==(rhs);
}

inline bool four_cc::operator==(const char* rhs) const
{
    return (memcmp(code, rhs, 4) == 0);
}

inline bool four_cc::operator!=(const char* rhs) const
{
    return !operator==(rhs);
}


inline four_cc::operator unsigned long() const
{
    return code_as_long;
}

inline unsigned long four_cc::as_long() const
{
    return code_as_long;
}

inline char four_cc::operator[](int i) const
{
    return code[i];
}

inline four_cc& four_cc::operator=(unsigned long val)
{
    code_as_long = val;
    return *this;
}

inline std::ostream& operator<<(std::ostream& os, const four_cc& rhs)
{
    return rhs.put(os);
}

#endif
