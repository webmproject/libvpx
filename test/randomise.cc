/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "test/randomise.h"

namespace libvpx_test {

// Add further specialisations as necessary

template<>
bool Randomise::uniform<bool>() {
  return rnd_.Rand8() & 1 ? true : false;
}

template<>
uint8_t Randomise::uniform<uint8_t>() {
  return rnd_.Rand8();
}

template<>
uint16_t Randomise::uniform<uint16_t>() {
  return rnd_.Rand16();
}

template<>
uint32_t Randomise::uniform<uint32_t>() {
  const uint32_t l = uniform<uint16_t>();
  const uint32_t h = uniform<uint16_t>();
  return h << 16 | l;
}

template<>
uint64_t Randomise::uniform<uint64_t>() {
  const uint64_t l = uniform<uint32_t>();
  const uint64_t h = uniform<uint32_t>();
  return h << 32 | l;
}

template<>
int8_t Randomise::uniform<int8_t>() { return uniform<uint8_t>(); }

template<>
int16_t Randomise::uniform<int16_t>() { return uniform<uint16_t>(); }

template<>
int32_t Randomise::uniform<int32_t>() { return uniform<uint32_t>(); }

template<>
int64_t Randomise::uniform<int64_t>() { return uniform<uint64_t>(); }

}  // namespace libvpx_test
