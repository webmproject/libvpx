#include <memory>
class SimpleEncode {
 public:
  SimpleEncode(int frame_width, int frame_height, vpx_rational_t frame_rate,
               int target_bitrate);
  ~SimpleEncode();
  SimpleEncode(SimpleEncode &&) = delete;
  SimpleEncode &operator=(SimpleEncode &&) = delete;

 private:
  class impl;
  std::unique_ptr<impl> pimpl;
};
