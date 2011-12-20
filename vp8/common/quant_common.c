/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "quant_common.h"


#if !CONFIG_EXTEND_QRANGE
static const int dc_qlookup[QINDEX_RANGE] =
{
    4,    5,    6,    7,    8,    9,    10,   10,   11,   12,   13,   14,   15,   16,   17,   17,
    18,   19,   20,   20,   21,   21,   22,   22,   23,   23,   24,   25,   25,   26,   27,   28,
    29,   30,   31,   32,   33,   34,   35,   36,   37,   37,   38,   39,   40,   41,   42,   43,
    44,   45,   46,   46,   47,   48,   49,   50,   51,   52,   53,   54,   55,   56,   57,   58,
    59,   60,   61,   62,   63,   64,   65,   66,   67,   68,   69,   70,   71,   72,   73,   74,
    75,   76,   76,   77,   78,   79,   80,   81,   82,   83,   84,   85,   86,   87,   88,   89,
    91,   93,   95,   96,   98,   100,  101,  102,  104,  106,  108,  110,  112,  114,  116,  118,
    122,  124,  126,  128,  130,  132,  134,  136,  138,  140,  143,  145,  148,  151,  154,  157,
};

static const int ac_qlookup[QINDEX_RANGE] =
{
    4,    5,    6,    7,    8,    9,    10,   11,   12,   13,   14,   15,   16,   17,   18,   19,
    20,   21,   22,   23,   24,   25,   26,   27,   28,   29,   30,   31,   32,   33,   34,   35,
    36,   37,   38,   39,   40,   41,   42,   43,   44,   45,   46,   47,   48,   49,   50,   51,
    52,   53,   54,   55,   56,   57,   58,   60,   62,   64,   66,   68,   70,   72,   74,   76,
    78,   80,   82,   84,   86,   88,   90,   92,   94,   96,   98,  100,  102,  104,  106,  108,
    110,  112,  114,  116,  119,  122,  125,  128,  131,  134,  137,  140,  143,  146,  149,  152,
    155,  158,  161,  164,  167,  170,  173,  177,  181,  185,  189,  193,  197,  201,  205,  209,
    213,  217,  221,  225,  229,  234,  239,  245,  249,  254,  259,  264,  269,  274,  279,  284,
};
#else
/*static int dc_qlookup[QINDEX_RANGE] =
{
      4,    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,   15,   16,   17,   18,   19,
     20,   22,   24,   26,   28,   30,   32,   34,   36,   38,   40,   38,   40,   42,   44,   47,
     50,   53,   56,   59,   62,   65,   68,   68,   72,   76,   80,   80,   84,   84,   86,   90,
     91,   96,   98,  102,  107,  112,  118,  124,  130,  136,  142,  148,  150,  156,  162,  168,
    174,  180,  182,  188,  196,  200,  204,  208,  212,  216,  220,  224,  228,  232,  236,  240,
    244,  249,  255,  260,  264,  269,  276,  281,  288,  293,  300,  303,  307,  312,  317,  323,
    328,  331,  338,  344,  347,  354,  366,  378,  386,  398,  401,  411,  422,  432,  441,  451,
    464,  483,  496,  507,  520,  529,  540,  551,  570,  585,  604,  624,  645,  667,  692,  718,

};
static int ac_qlookup[QINDEX_RANGE] =
{
    4,    5,    6,    7,    8,    9,    10,   11,   12,   13,   14,   15,   16,   17,   18,   19,
    20,   22,   24,   26,   28,   30,   32,   34,   36,   38,   40,   42,   44,   46,   48,   51,
    54,   57,   60,   63,   66,   69,   72,   76,   80,   84,   88,   92,   96,   100,  105,  110,
    115,  120,  125,  130,  135,  140,  146,  152,  158,  164,  170,  176,  182,  188,  194,  200,
    206,  212,  218,  224,  232,  240,  248,  256,  264,  272,  280,  288,  296,  304,  312,  320,
    330,  340,  350,  360,  370,  380,  392,  404,  416,  428,  440,  454,  468,  482,  496,  510,
    524,  540,  556,  572,  588,  604,  622,  640,  658,  676,  696,  716,  736,  756,  776,  796,
    820,  844,  868,  892,  916,  944,  972,  1000, 1032, 1064, 1096, 1128, 1168, 1208, 1252, 1300
};*/

static int dc_qlookup[QINDEX_RANGE] =
{
    4,    5,    6,    7,    8,    9,    10,   10,   11,   12,   13,   14,   15,   16,   17,   17,
    18,   19,   20,   20,   21,   21,   22,   22,   23,   23,   24,   25,   25,   26,   27,   28,
    29,   30,   31,   32,   33,   34,   35,   36,   37,   37,   38,   39,   40,   41,   42,   43,
    44,   45,   46,   46,   47,   48,   49,   50,   51,   52,   53,   54,   55,   56,   57,   58,
    59,   60,   61,   62,   63,   64,   65,   66,   67,   68,   69,   70,   71,   72,   73,   74,
    75,   76,   76,   77,   78,   79,   80,   81,   82,   83,   84,   85,   86,   87,   88,   89,
    91,   93,   95,   96,   98,   100,  101,  102,  104,  106,  108,  110,  112,  114,  116,  118,
    122,  124,  126,  128,  130,  132,  134,  136,  138,  140,  143,  145,  148,  151,  154,  157,
};

static int ac_qlookup[QINDEX_RANGE] =
{
    4,    5,    6,    7,    8,    9,    10,   11,   12,   13,   14,   15,   16,   17,   18,   19,
    20,   21,   22,   23,   24,   25,   26,   27,   28,   29,   30,   31,   32,   33,   34,   35,
    36,   37,   38,   39,   40,   41,   42,   43,   44,   45,   46,   47,   48,   49,   50,   51,
    52,   53,   54,   55,   56,   57,   58,   60,   62,   64,   66,   68,   70,   72,   74,   76,
    78,   80,   82,   84,   86,   88,   90,   92,   94,   96,   98,  100,  102,  104,  106,  108,
    110,  112,  114,  116,  119,  122,  125,  128,  131,  134,  137,  140,  143,  146,  149,  152,
    155,  158,  161,  164,  167,  170,  173,  177,  181,  185,  189,  193,  197,  201,  205,  209,
    213,  217,  221,  225,  229,  234,  239,  245,  249,  254,  259,  264,  269,  274,  279,  284,
};

//static int dc_qlookup[QINDEX_RANGE];
//static int ac_qlookup[QINDEX_RANGE];

#endif

#if CONFIG_EXTEND_QRANGE
#define ACDC_MIN 4
void vp8_init_quant_tables()
{
    int i;
    int current_val = 4;
    int last_val = 4;
    int ac_val;
    int dc_max;

    for ( i = 0; i < QINDEX_RANGE; i++ )
    {
        ac_qlookup[i] = ac_qlookup[i] << 2;
        dc_qlookup[i] = dc_qlookup[i] << 2;
    }

    // Not active by default for now.
    return;

    for ( i = 0; i < QINDEX_RANGE; i++ )
    {
        ac_qlookup[i] = current_val;
        current_val = (int)((double)current_val * 1.042);
        //current_val = (int)((double)current_val * 1.01765);
        if ( current_val == last_val )
            current_val++;
        last_val = current_val;

        ac_val = ac_qlookup[i];
        dc_max = (int)(((double)ac_val * 0.75) + 0.5);
        dc_qlookup[i] = (0.000000305 * ac_val * ac_val * ac_val) +
                        (-0.00065 * ac_val * ac_val) +
                        (0.9 * ac_val) + 0.5;
        if ( dc_qlookup[i] > dc_max )
            dc_qlookup[i] = dc_max;
        if ( dc_qlookup[i] < ACDC_MIN )
            dc_qlookup[i] = ACDC_MIN;
    }
}
#endif

int vp8_dc_quant(int QIndex, int Delta)
{
    int retval;

    QIndex = QIndex + Delta;

    if (QIndex > MAXQ)
        QIndex = MAXQ;
    else if (QIndex < 0)
        QIndex = 0;

    retval = dc_qlookup[ QIndex ];
    return retval;
}

int vp8_dc2quant(int QIndex, int Delta)
{
    int retval;

    QIndex = QIndex + Delta;

    if (QIndex > MAXQ)
        QIndex = MAXQ;
    else if (QIndex < 0)
        QIndex = 0;

#if !CONFIG_EXTEND_QRANGE
    retval = dc_qlookup[ QIndex ] * 2;
#else
    retval = dc_qlookup[ QIndex ];
#endif
    return retval;

}
int vp8_dc_uv_quant(int QIndex, int Delta)
{
    int retval;

    QIndex = QIndex + Delta;

    if (QIndex > MAXQ)
        QIndex = MAXQ;
    else if (QIndex < 0)
        QIndex = 0;

    retval = dc_qlookup[ QIndex ];

    return retval;
}

int vp8_ac_yquant(int QIndex)
{
    int retval;

    if (QIndex > MAXQ)
        QIndex = MAXQ;
    else if (QIndex < 0)
        QIndex = 0;

    retval = ac_qlookup[ QIndex ];
    return retval;
}

int vp8_ac2quant(int QIndex, int Delta)
{
    int retval;

    QIndex = QIndex + Delta;

    if (QIndex > MAXQ)
        QIndex = MAXQ;
    else if (QIndex < 0)
        QIndex = 0;
#if !CONFIG_EXTEND_QRANGE
    retval = (ac_qlookup[ QIndex ] * 155) / 100;
    if (retval < 8)
        retval = 8;
#else
    retval = (ac_qlookup[ QIndex ] * 775) / 1000;
    if (retval < 4)
        retval = 4;
#endif
    return retval;
}
int vp8_ac_uv_quant(int QIndex, int Delta)
{
    int retval;

    QIndex = QIndex + Delta;

    if (QIndex > MAXQ)
        QIndex = MAXQ;
    else if (QIndex < 0)
        QIndex = 0;

    retval = ac_qlookup[ QIndex ];
    return retval;
}
