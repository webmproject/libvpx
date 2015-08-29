#include <stdlib.h>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "test/acm_random.h"
#include "test/util.h"
#include "vp10/encoder/dct.c"

using libvpx_test::ACMRandom;

namespace {
void reference_dct_1d(const double *in, double *out, int size) {
  const double PI = 3.141592653589793238462643383279502884;
  const double kInvSqrt2 = 0.707106781186547524400844362104;
  for (int k = 0; k < size; ++k) {
      out[k] = 0;  // initialize out[k]
      for (int n = 0; n < size; ++n) {
        out[k] += in[n] * cos(PI * (2 * n + 1) * k / (2*size));
      }
      if (k == 0)
          out[k] = out[k] * kInvSqrt2;
  }
}

typedef void (*FdctFuncRef)(const double *in, double *out, int size);
typedef void (*IdctFuncRef)(const double *in, double *out, int size);
typedef void (*FdctFunc)(const tran_low_t *in, tran_low_t *out);
typedef void (*IdctFunc)(const tran_low_t *in, tran_low_t *out);

class TransTestBase {
 public:
  virtual ~TransTestBase() {}

 protected:
  double max_error;
  int txfmSize;
  FdctFunc fwd_txfm;
  FdctFuncRef fwd_txfm_ref;
  virtual void RunFwdAccuracyCheck() {
    tran_low_t* input  = new tran_low_t[txfmSize];
    tran_low_t* output = new tran_low_t[txfmSize];
    double* refInput   = new double[txfmSize];
    double* refOutput  = new double[txfmSize];

    ACMRandom rnd(ACMRandom::DeterministicSeed());
    const int count_test_block = 5000;
    for (int ti =  0; ti < count_test_block; ++ti) {
      for (int ni = 0; ni < txfmSize; ++ni) {
        input[ni] = rnd.Rand8() - rnd.Rand8();
        refInput[ni] = (double)input[ni];
      }

      fwd_txfm(input, output);
      fwd_txfm_ref(refInput, refOutput, txfmSize);

      for (int ni = 0; ni < txfmSize; ++ni) {
        EXPECT_LE(abs(output[ni]-(tran_low_t)round(refOutput[ni])),
                  max_error);
      }
    }

    delete[] input;
    delete[] output;
    delete[] refInput;
    delete[] refOutput;
  }
};

typedef std::tr1::tuple<FdctFunc, FdctFuncRef, int, int> FdctParam;
class Vp10FwdDct
  : public TransTestBase,
    public ::testing::TestWithParam<FdctParam> {
 public:
      virtual void SetUp() {
        fwd_txfm = GET_PARAM(0);
        fwd_txfm_ref = GET_PARAM(1);
        txfmSize = GET_PARAM(2);
        max_error = GET_PARAM(3);
      }
      virtual void TearDown(){}
};

TEST_P(Vp10FwdDct, RunFwdAccuracyCheck) {
  RunFwdAccuracyCheck();
}

INSTANTIATE_TEST_CASE_P(
    C, Vp10FwdDct,
    ::testing::Values(
        FdctParam(&fdct4, &reference_dct_1d, 4, 1),
        FdctParam(&fdct8, &reference_dct_1d, 8, 1),
        FdctParam(&fdct16, &reference_dct_1d, 16, 2),
        FdctParam(&fdct32, &reference_dct_1d, 32, 4))
);
}  // namespace

