#if defined(EXPERIMENTAL_C)
/* The experimental.c file includes this file multiple times to build up the
 * required state.
 */
#if !defined(XXX_ARG_DEF)
#define XXX_ARG_DEF(sym, value) \
    static const arg_def_t xxx_arg_def_##sym = \
        ARG_DEF(NULL, #sym, 1, "Experimental");

#define XXX_DEFINE_INT(sym, value) \
    XXX_ARG_DEF(sym, value); int xxx_##sym = value;
#define XXX_DEFINE_UINT(sym, value) \
    XXX_ARG_DEF(sym, value); unsigned int xxx_##sym = value;

#elif !defined(XXX_ARG_DEF_LIST)
#define XXX_ARG_DEF_LIST(sym) &xxx_arg_def_##sym,

#undef  XXX_DEFINE_INT
#define XXX_DEFINE_INT(sym, value) XXX_ARG_DEF_LIST(sym)

#undef  XXX_DEFINE_UINT
#define XXX_DEFINE_UINT(sym, value) XXX_ARG_DEF_LIST(sym)

#elif !defined(XXX_ARG_MATCH)
#define XXX_ARG_MATCH

#undef  XXX_DEFINE_INT
#define XXX_DEFINE_INT(sym, value)\
    else if (arg_match(&arg, &xxx_arg_def_##sym, argi)) \
        xxx_##sym = arg_parse_int(&arg);

#undef  XXX_DEFINE_UINT
#define XXX_DEFINE_UINT(sym, value)\
    else if (arg_match(&arg, &xxx_arg_def_##sym, argi)) \
        xxx_##sym = arg_parse_uint(&arg);

#endif
#else
/* All other files just get the extern references to these symbols. */

#define XXX_DEFINE_INT(sym, value) extern int xxx_##sym;
#define XXX_DEFINE_UINT(sym, value) extern unsigned int xxx_##sym;


#include <stdio.h>
void xxx_show_usage(FILE *fp);
int xxx_parse_arg(char **argi);
#endif

/*
 * BEGIN EXPERIMENTS BELOW
 *
 * XXX_DEFINE_INT(knob, 0)
 */
XXX_DEFINE_INT(foo, 0)
XXX_DEFINE_INT(bar, 0)
