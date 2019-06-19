/*The dataset is from
 *Computer Vision Group
 *TUM Department of Informatics Technical
 *University of Munich
 *https://vision.in.tum.de/data/datasets/rgbd-dataset/download#freiburg1_xyz
 */
PointCloud point_cloud;  // pointCloud object
Camera cam;              // build camera object
MotionField mf;          // motion field
void setup() {
  size(640, 480, P3D);
  // basic settings
  float focal = 525.0f;        // focal distance of camera
  int frame_no = 0;            // frame number
  float fov = PI / 3;          // field of view
  int block_size = 8;          // block size
  float normalizer = 5000.0f;  // normalizer
  // initialize
  point_cloud = new PointCloud();
  // synchronized rgb, depth and ground truth
  String head = "../data/";
  String[] rgb_depth_gt = loadStrings(head + "rgb_depth_groundtruth.txt");

  // read in rgb and depth image file paths as well as corresponding camera
  // posiiton and quaternion
  String[] info = split(rgb_depth_gt[frame_no], ' ');
  String rgb_path = head + info[1];
  String depth_path = head + info[3];
  float tx = float(info[7]), ty = float(info[8]),
        tz = float(info[9]);  // real camera position
  float qx = float(info[10]), qy = float(info[11]), qz = float(info[12]),
        qw = float(info[13]);  // quaternion

  // build transformer
  Transform trans = new Transform(tx, ty, tz, qx, qy, qz, qw, focal, width,
                                  height, normalizer);
  PImage rgb = loadImage(rgb_path);
  PImage depth = loadImage(depth_path);
  // generate point cloud
  point_cloud.generate(rgb, depth, trans);
  // get the center of cloud
  PVector cloud_center = point_cloud.getCloudCenter();
  // initialize camera
  cam =
      new Camera(fov, new PVector(0, 0, 0), cloud_center, new PVector(0, 1, 0));
  // initialize motion field
  mf = new MotionField(cam, point_cloud, block_size);
}
void draw() {
  background(0);
  // run camera dragged mouse to rotate camera
  // w: go forward
  // s: go backward
  // a: go left
  // d: go right
  // up arrow: go up
  // down arrow: go down
  //+ increase move speed
  //- decrease move speed
  // r: rotate the camera
  // b: reset to initial position
  cam.run();
  // render the point lists
  point_cloud.render();
  // update motion field
  mf.run();
  // draw motion field
  mf.showMotionField();
}
