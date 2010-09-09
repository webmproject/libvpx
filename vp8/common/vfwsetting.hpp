/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#if !defined(VFWSETTING_HPP)
#define VFWSETTING_HPP
//______________________________________________________________________________
//
//  VFWSetting.hpp
//

#include "four_cc.hpp"
#include <iosfwd>

namespace vpxvp
{

    //--------------------------------------
    class VFWSetting
    {
        friend std::ostream& operator<<(std::ostream& os, const VFWSetting& vfws);

    public:

        enum Mode
        {
            m_setting,
            m_config
        };

        enum
        {
            header_size = 8,
            Size = 16
        };

        VFWSetting(four_cc fcc);
        ~VFWSetting();

        four_cc fcc() const;
        Mode mode() const;

        int setting() const;
        int value() const;
        void setting_value(int i_setting, int i_value);  //  Sets mode to m_setting

        long size() const;
        const void* data() const;
        int data(const void* p_data, unsigned long ul_size);

    private:

        VFWSetting(const VFWSetting& vfws);  //  Not implemented
        VFWSetting& operator=(const VFWSetting& vfws);  //  Not implemented

        int extract_(const void* p_data, unsigned long ul_size);
        void update_() const;

        four_cc m_fcc;
        Mode m_mode;
        int m_i_setting;
        int m_i_value;

        mutable unsigned char m_p_data[Size];
    };

}  //  namespace vpxvp

#endif  //  VFWSETTING_HPP
