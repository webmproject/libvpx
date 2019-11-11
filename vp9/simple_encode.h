#include <memory>
#include <vector>

enum FrameType {
  kKeyFrame = 0,
  kInterFrame,
  kAlternateReference,
};

struct EncodeFrameResult {
  // TODO(angiebird): int show_index;
  // TODO(angiebird): FrameType frame_type;
  size_t coding_data_size;
  // The EncodeFrame will allocate a buffer, write the coding data into the
  // buffer and give the ownership of the buffer to coding_data
  std::unique_ptr<unsigned char[]> coding_data;
  // TODO(angiebird): double psnr;
  // TODO(angiebird): int quantize_index ;
};

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
  void EncodeFrame(EncodeFrameResult *encode_frame_result);

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
