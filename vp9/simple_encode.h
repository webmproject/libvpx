#include <memory>
#include <vector>
class SimpleEncode {
 public:
  SimpleEncode(int frame_width, int frame_height, int frame_rate_num,
               int frame_rate_den, int target_bitrate, int num_frames,
               FILE *file);
  ~SimpleEncode();
  SimpleEncode(SimpleEncode &&) = delete;
  SimpleEncode &operator=(SimpleEncode &&) = delete;

  // Make encoder compute the first pass stats and store it internally for
  // future encode
  void ComputeFirstPassStats();

  // Output the first pass stats
  std::vector<std::vector<double>> ObserveFirstPassStats();

  // Initialize the encoder for actual encoding
  void StartEncode();

  // Free the encoder
  void EndEncode();

  // Encode a frame
  void EncodeFrame(char *cx_data, size_t *size, size_t max_size);

 private:
  class impl;
  int frame_width;
  int frame_height;
  int frame_rate_num;
  int frame_rate_den;
  int target_bitrate;
  int num_frames;
  FILE *file;
  std::unique_ptr<impl> pimpl;
};
