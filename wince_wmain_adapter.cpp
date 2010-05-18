/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


/* This program is created to take command arguments and pass
 * them to main() in example.c or example_xma.c, because the
 * correspending part in example.c or example_xma.c does not
 * work on Pocket PC platform.
 * To modify the command arguments, go to "Property" page and
 * fill in "Command arguments." For example:
 *  --codec vp6 --flipuv --progress _bnd.vp6
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_NUM_ARG 64
#define MAX_SIZ_ARG 512

extern "C"
{
    int main(int argc, char **argv);
}

int wmain(int argc, wchar_t **argv) {
    char *cargv[MAX_NUM_ARG];
    char chargv[MAX_SIZ_ARG];
    int ret;

    /* transform command line arguments from (wchar_t *) to (char *) */
    for(int i=0; i<argc; i++) {
        wcstombs( chargv, argv[i], sizeof(chargv));
        cargv[i] = _strdup(chargv);
    }

    ret = main(argc, (char **)cargv);

    //free the memory located by _strdup()
    for(int i=0; i<argc; i++)
        free(cargv[i]);

    return ret;
}
