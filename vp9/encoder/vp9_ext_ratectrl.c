/*
 *  Copyright (c) 2020 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/encoder/vp9_ext_ratectrl.h"
#include "vp9/common/vp9_common.h"

void vp9_extrc_init(EXT_RATECTRL *ext_ratectrl) { vp9_zero(*ext_ratectrl); }

void vp9_extrc_create(vpx_rc_funcs_t funcs, vpx_rc_config_t ratectrl_config,
                      EXT_RATECTRL *ext_ratectrl) {
  vp9_extrc_delete(ext_ratectrl);
  ext_ratectrl->funcs = funcs;
  ext_ratectrl->ratectrl_config = ratectrl_config;
  ext_ratectrl->funcs.create_model(ext_ratectrl->funcs.priv,
                                   &ext_ratectrl->ratectrl_config,
                                   ext_ratectrl->model);
  ext_ratectrl->ready = 1;
}

void vp9_extrc_delete(EXT_RATECTRL *ext_ratectrl) {
  if (ext_ratectrl->ready) {
    ext_ratectrl->funcs.delete_model(ext_ratectrl->model);
  }
  vp9_extrc_init(ext_ratectrl);
}
