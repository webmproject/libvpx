#include "vpx_ports/config.h"

#include "../../common/opencl/vp8_opencl.h"
#include "vp8_decode_cl.h"

#include <stdio.h>

extern int cl_init_dequant();
extern int cl_destroy_dequant();

int cl_decode_destroy(){

#if ENABLE_CL_IDCT_DEQUANT
    int err;
    err = cl_destroy_dequant();
#endif
    
    return CL_SUCCESS;
}

int cl_decode_init()
{
#if ENABLE_CL_IDCT_DEQUANT
    int err;
#endif

    //Initialize programs to null value
    //Enables detection of if they've been initialized as well.
    cl_data.dequant_program = NULL;

#if ENABLE_CL_IDCT_DEQUANT
    err = cl_init_dequant();
    if (err != CL_SUCCESS)
        return err;
#endif

    return CL_SUCCESS;
}
