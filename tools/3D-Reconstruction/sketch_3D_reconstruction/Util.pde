// show grids
void showGrids(int block_size) {
  ortho(-width, 0, -height, 0);
  camera(0, 0, 0, 0, 0, 1, 0, 1, 0);
  stroke(0, 0, 255);
  for (int i = 0; i < height; i += block_size) {
    line(0, i, width, i);
  }
  for (int i = 0; i < width; i += block_size) {
    line(i, 0, i, height);
  }
}

// save the point clould information
void savePointCloud(PointCloud point_cloud, String file_name) {
  String[] results = new String[point_cloud.points.size()];
  for (int i = 0; i < point_cloud.points.size(); i++) {
    PVector point = point_cloud.points.get(i);
    results[i] = str(point.x) + ' ' + str(point.y) + ' ' + str(point.z);
  }
  saveStrings(file_name, results);
}
