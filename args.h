/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef ARGS_H
#define ARGS_H
#include <stdio.h>

struct arg
{
    char                 **argv;
    const char            *name;
    const char            *val;
    unsigned int           argv_step;
    const struct arg_def  *def;
};

typedef struct arg_def
{
    const char *short_name;
    const char *long_name;
    int         has_val;
    const char *desc;
} arg_def_t;
#define ARG_DEF(s,l,v,d) {s,l,v,d}
#define ARG_DEF_LIST_END {0}

struct arg arg_init(char **argv);
int arg_match(struct arg *arg_, const struct arg_def *def, char **argv);
const char *arg_next(struct arg *arg);
void arg_show_usage(FILE *fp, const struct arg_def *const *defs);
char **argv_dup(int argc, const char **argv);

unsigned int arg_parse_uint(const struct arg *arg);
int arg_parse_int(const struct arg *arg);
struct vpx_rational arg_parse_rational(const struct arg *arg);
#endif
