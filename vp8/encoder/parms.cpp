/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#if 0

#include <map>
#include <string>
#include <fstream>
extern "C"
{
    #include "vp8/common/onyx.h"
}


using namespace std;

typedef map<string,int> Parms;

#define ALLPARMS(O,DOTHIS) \
    DOTHIS(O,  interquantizer             )\
    DOTHIS(O,  auto_gold                   )\
    DOTHIS(O,  auto_adjust_gold_quantizer    )\
    DOTHIS(O,  goldquantizer              )\
    DOTHIS(O,  goldfreq                   )\
    DOTHIS(O,  auto_key                    )\
    DOTHIS(O,  auto_adjust_key_quantizer     )\
    DOTHIS(O,  keyquantizer               )\
    DOTHIS(O,  keyfreq                    )\
    DOTHIS(O,  pass                       )\
    DOTHIS(O,  fixed_q                     )\
    DOTHIS(O,  target_bandwidth            )\
    DOTHIS(O,  auto_worst_q                 )\
    DOTHIS(O,  worst_quality               )\
    DOTHIS(O,  best_allowed_q               )\
    DOTHIS(O,  end_usage                   )\
    DOTHIS(O,  starting_buffer_level        )\
    DOTHIS(O,  optimal_buffer_level         )\
    DOTHIS(O,  maximum_buffer_size          )\
    DOTHIS(O,  under_shoot_pct              )\
    DOTHIS(O,  allow_df                    )\
    DOTHIS(O,  drop_frames_water_mark        )\
    DOTHIS(O,  max_allowed_datarate         )\
    DOTHIS(O,  two_pass_vbrbias             )\
    DOTHIS(O,  two_pass_vbrmin_section       )\
    DOTHIS(O,  two_pass_vbrmax_section       )\
    DOTHIS(O,  filter_type                 )\
    DOTHIS(O,  compressor_speed            )\
    DOTHIS(O,  mbpitch_feature             )\
    DOTHIS(O,  allow_spatial_resampling     )\
    DOTHIS(O,  resample_down_water_mark      )\
    DOTHIS(O,  resample_up_water_mark        )\
    DOTHIS(O,  noise_sensitivity           )\
    DOTHIS(O,  horiz_scale                 )\
    DOTHIS(O,  vert_scale                  )


#define GET(O,V) O->V = x[#V];
#define PUT(O,V) x[#V] = O->V;


extern "C" void get_parms(VP8_CONFIG *ocf,char *filename)
{

    Parms x;
    int value;
    string variable;
    string equal;

    ifstream config_file(filename);

    ALLPARMS(ocf, PUT);

    // store all the parms in a map (really simple parsing)
    while(!config_file.eof() && config_file.is_open())
    {
        config_file >> variable;
        config_file >> equal;

        if(equal != "=")
            continue;

        config_file >> value;

        x[variable] = value;
    }

    ALLPARMS(ocf, GET);

}

#define PRINT(O,V) debug_file<<#V <<" = " << O->V <<"\n";
extern "C" void print_parms(VP8_CONFIG *ocf,char *filename)
{
    ofstream debug_file(filename,ios_base::app);
    ALLPARMS(ocf, PRINT);
    debug_file << "=============================================="<<"\n";
}

#endif
