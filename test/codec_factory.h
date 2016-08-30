/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#ifndef TEST_CODEC_FACTORY_H_
#define TEST_CODEC_FACTORY_H_

#include "./aom_config.h"
#include "aom/aom_decoder.h"
#include "aom/aom_encoder.h"
#if CONFIG_AV1_ENCODER
#include "aom/aomcx.h"
#endif
#if CONFIG_AV1_DECODER
#include "aom/aomdx.h"
#endif

#include "test/decode_test_driver.h"
#include "test/encode_test_driver.h"
namespace libaom_test {

const int kCodecFactoryParam = 0;

class CodecFactory {
 public:
  CodecFactory() {}

  virtual ~CodecFactory() {}

  virtual Decoder *CreateDecoder(aom_codec_dec_cfg_t cfg,
                                 unsigned long deadline) const = 0;

  virtual Decoder *CreateDecoder(aom_codec_dec_cfg_t cfg,
                                 const aom_codec_flags_t flags,
                                 unsigned long deadline)  // NOLINT(runtime/int)
      const = 0;

  virtual Encoder *CreateEncoder(aom_codec_enc_cfg_t cfg,
                                 unsigned long deadline,
                                 const unsigned long init_flags,
                                 TwopassStatsStore *stats) const = 0;

  virtual aom_codec_err_t DefaultEncoderConfig(aom_codec_enc_cfg_t *cfg,
                                               int usage) const = 0;
};

/* Provide CodecTestWith<n>Params classes for a variable number of parameters
 * to avoid having to include a pointer to the CodecFactory in every test
 * definition.
 */
template <class T1>
class CodecTestWithParam
    : public ::testing::TestWithParam<
          std::tr1::tuple<const libaom_test::CodecFactory *, T1> > {};

template <class T1, class T2>
class CodecTestWith2Params
    : public ::testing::TestWithParam<
          std::tr1::tuple<const libaom_test::CodecFactory *, T1, T2> > {};

template <class T1, class T2, class T3>
class CodecTestWith3Params
    : public ::testing::TestWithParam<
          std::tr1::tuple<const libaom_test::CodecFactory *, T1, T2, T3> > {};

/*
 * AV1 Codec Definitions
 */
#if CONFIG_AV1
class AV1Decoder : public Decoder {
 public:
  AV1Decoder(aom_codec_dec_cfg_t cfg, unsigned long deadline)
      : Decoder(cfg, deadline) {}

  AV1Decoder(aom_codec_dec_cfg_t cfg, const aom_codec_flags_t flag,
             unsigned long deadline)  // NOLINT
      : Decoder(cfg, flag, deadline) {}

 protected:
  virtual aom_codec_iface_t *CodecInterface() const {
#if CONFIG_AV1_DECODER
    return &aom_codec_av1_dx_algo;
#else
    return NULL;
#endif
  }
};

class AV1Encoder : public Encoder {
 public:
  AV1Encoder(aom_codec_enc_cfg_t cfg, unsigned long deadline,
             const unsigned long init_flags, TwopassStatsStore *stats)
      : Encoder(cfg, deadline, init_flags, stats) {}

 protected:
  virtual aom_codec_iface_t *CodecInterface() const {
#if CONFIG_AV1_ENCODER
    return &aom_codec_av1_cx_algo;
#else
    return NULL;
#endif
  }
};

class AV1CodecFactory : public CodecFactory {
 public:
  AV1CodecFactory() : CodecFactory() {}

  virtual Decoder *CreateDecoder(aom_codec_dec_cfg_t cfg,
                                 unsigned long deadline) const {
    return CreateDecoder(cfg, 0, deadline);
  }

  virtual Decoder *CreateDecoder(aom_codec_dec_cfg_t cfg,
                                 const aom_codec_flags_t flags,
                                 unsigned long deadline) const {  // NOLINT
#if CONFIG_AV1_DECODER
    return new AV1Decoder(cfg, flags, deadline);
#else
    return NULL;
#endif
  }

  virtual Encoder *CreateEncoder(aom_codec_enc_cfg_t cfg,
                                 unsigned long deadline,
                                 const unsigned long init_flags,
                                 TwopassStatsStore *stats) const {
#if CONFIG_AV1_ENCODER
    return new AV1Encoder(cfg, deadline, init_flags, stats);
#else
    return NULL;
#endif
  }

  virtual aom_codec_err_t DefaultEncoderConfig(aom_codec_enc_cfg_t *cfg,
                                               int usage) const {
#if CONFIG_AV1_ENCODER
    return aom_codec_enc_config_default(&aom_codec_av1_cx_algo, cfg, usage);
#else
    return AOM_CODEC_INCAPABLE;
#endif
  }
};

const libaom_test::AV1CodecFactory kAV1;

#define AV1_INSTANTIATE_TEST_CASE(test, ...)                                \
  INSTANTIATE_TEST_CASE_P(                                                  \
      AV1, test,                                                            \
      ::testing::Combine(                                                   \
          ::testing::Values(static_cast<const libaom_test::CodecFactory *>( \
              &libaom_test::kAV1)),                                         \
          __VA_ARGS__))
#else
#define AV1_INSTANTIATE_TEST_CASE(test, ...)
#endif  // CONFIG_AV1

}  // namespace libaom_test
#endif  // TEST_CODEC_FACTORY_H_
