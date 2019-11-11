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
  // This funtion should be called after ComputeFirstPassStats()
  void StartEncode();

  // Free the encoder
  // This funtion should be called after StartEncode() or EncodeFrame()
  void EndEncode();

  // Encode a frame
  // This funtion should be called after StartEncode() before EndEncode()
  void EncodeFrame(char *cx_data, size_t *size, size_t max_size);

  // Get the number of coding frames for the video. The coding frames include
  // show frame and no show frame.
  // This funtion should be called after ComputeFirstPassStats()
  int GetCodingFrameNum();

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
