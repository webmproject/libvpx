#define EXPERIMENTAL_C
#include <stdio.h>

#include "args.h"

/* Get argument definitions */
#include "experimental.h"

/* Build argument definition list */
static const arg_def_t *xxx_def_list[] = {
#include "experimental.h"
NULL
};

void xxx_show_usage(FILE *fp)
{
    arg_show_usage(fp, xxx_def_list);
}

int xxx_parse_arg(char **argi)
{
    struct arg arg;

    arg = arg_init(argi);
    if(0);
#include "experimental.h"
    else return 0;
    return 1;
}
