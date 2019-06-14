class Scene {
  Camera camera;
  PointCloud point_cloud;
  MotionField motion_field;
  ArrayList<PVector> last_positions;
  ArrayList<PVector> current_positions;
  int[] render_list;

  Scene(Camera camera, PointCloud point_cloud, MotionField motion_field) {
    this.camera = camera;
    this.point_cloud = point_cloud;
    this.motion_field = motion_field;
    last_positions = project2Camera();
    current_positions = project2Camera();
    render_list = new int[width * height];
    updateRenderList();
  }

  ArrayList<PVector> project2Camera() {
    ArrayList<PVector> projs = new ArrayList<PVector>();
    for (int i = 0; i < point_cloud.size(); i++) {
      PVector proj = camera.project(point_cloud.getPosition(i));
      projs.add(proj);
    }
    return projs;
  }
  // update render list by using depth test
  void updateRenderList() {
    // clear render list
    for (int i = 0; i < width * height; i++) render_list[i] = -1;
    // depth test and get render list
    float[] depth = new float[width * height];
    for (int i = 0; i < width * height; i++) depth[i] = Float.POSITIVE_INFINITY;
    for (int i = 0; i < current_positions.size(); i++) {
      PVector pos = current_positions.get(i);
      int row = int(pos.y + height / 2);
      int col = int(pos.x + width / 2);
      int idx = row * width + col;
      if (row >= 0 && row < height && col >= 0 && col < width) {
        if (render_list[idx] == -1 || pos.z < depth[idx]) {
          depth[idx] = pos.z;
          render_list[idx] = i;
        }
      }
    }
  }

  void run() {
    camera.run();
    last_positions = current_positions;
    current_positions = project2Camera();
    updateRenderList();
    motion_field.update(last_positions, current_positions, render_list);
  }

  void render(boolean show_motion_field) {
    // build mesh
    camera.open();
    noStroke();
    beginShape(TRIANGLES);
    for (int i = 0; i < height - 1; i++)
      for (int j = 0; j < width - 1; j++) {
        PVector pos0 = point_cloud.getPosition(i * width + j);
        PVector pos1 = point_cloud.getPosition(i * width + j + 1);
        PVector pos2 = point_cloud.getPosition((i + 1) * width + j + 1);
        PVector pos3 = point_cloud.getPosition((i + 1) * width + j);
        fill(point_cloud.getColor(i * width + j));
        vertex(pos0.x, pos0.y, pos0.z);
        fill(point_cloud.getColor(i * width + j + 1));
        vertex(pos1.x, pos1.y, pos1.z);
        fill(point_cloud.getColor((i + 1) * width + j + 1));
        vertex(pos2.x, pos2.y, pos2.z);

        fill(point_cloud.getColor((i + 1) * width + j + 1));
        vertex(pos2.x, pos2.y, pos2.z);
        fill(point_cloud.getColor((i + 1) * width + j + 1));
        vertex(pos3.x, pos3.y, pos3.z);
        fill(point_cloud.getColor(i * width + j));
        vertex(pos0.x, pos0.y, pos0.z);
      }
    endShape();
    if (show_motion_field) {
      camera.close();
      motion_field.render();
    }
  }
}
