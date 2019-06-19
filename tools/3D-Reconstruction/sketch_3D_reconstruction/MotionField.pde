class MotionField {
  Camera camera;
  int block_size;
  PointCloud point_cloud;
  ArrayList<PVector> last_positions;
  ArrayList<PVector> current_positions;
  ArrayList<PVector> motion_field;

  MotionField(Camera camera, PointCloud point_cloud, int block_size) {
    this.camera = camera;
    this.point_cloud = point_cloud;
    this.block_size = block_size;
    this.last_positions = new ArrayList<PVector>();
    this.current_positions = new ArrayList<PVector>();
    this.motion_field = new ArrayList<PVector>();
    for (int i = 0; i < point_cloud.points.size(); i++) {
      PVector pos = point_cloud.points.get(i);
      PVector proj_pos = getProjectedPos(pos);
      last_positions.add(proj_pos);
      current_positions.add(proj_pos);
    }
  }

  PVector getProjectedPos(PVector pos) {
    float[] cam_mat = camera.getCameraMat();
    PVector trans_pos =
        PVector.sub(pos, camera.pos);  // translate based on camera's position
    PVector rot_pos =
        MatxVec3(cam_mat, trans_pos);  // rotate based on camera angle
    PVector proj_pos = new PVector(0, 0, 0);
    proj_pos.x =
        height / 2.0f * rot_pos.x / (rot_pos.z) / tan(camera.fov / 2.0f);
    proj_pos.y =
        height / 2.0f * rot_pos.y / (rot_pos.z) / tan(camera.fov / 2.0f);
    proj_pos.z = trans_pos.z;

    return proj_pos;
  }

  void run() {
    last_positions = current_positions;
    // take care
    current_positions = new ArrayList<PVector>();
    for (int i = 0; i < point_cloud.points.size(); i++) {
      PVector pos = point_cloud.points.get(i);
      PVector proj_pos = getProjectedPos(pos);
      current_positions.add(proj_pos);
    }
  }

  ArrayList<PVector> getMotionField() {
    ArrayList<PVector> depth = new ArrayList<PVector>();
    IntList pixel_idx = new IntList();
    for (int i = 0; i < width * height; i++)
      depth.add(new PVector(Float.POSITIVE_INFINITY, -1));
    for (int i = 0; i < current_positions.size(); i++) {
      PVector pos = current_positions.get(i);
      int row = int(pos.y + height / 2);
      int col = int(pos.x + width / 2);
      int idx = row * width + col;
      if (row >= 0 && row < height && col >= 0 && col < width) {
        PVector depth_info = depth.get(idx);
        if (depth_info.y == -1) {
          depth.set(idx, new PVector(pos.z, pixel_idx.size()));
          pixel_idx.append(i);
        } else if (pos.z < depth_info.x) {
          depth.set(idx, new PVector(pos.z, depth_info.y));
          pixel_idx.set(int(depth_info.y), i);
        }
      }
    }
    motion_field = new ArrayList<PVector>();
    int r_num = height / block_size, c_num = width / block_size;
    for (int i = 0; i < r_num * c_num; i++)
      motion_field.add(new PVector(0, 0, 0));
    for (int i = 0; i < pixel_idx.size(); i++) {
      PVector cur_pos = current_positions.get(pixel_idx.get(i));
      PVector last_pos = last_positions.get(pixel_idx.get(i));
      int row = int(cur_pos.y + height / 2);
      int col = int(cur_pos.x + width / 2);
      int idx = row / block_size * c_num + col / block_size;
      PVector mv = PVector.sub(last_pos, cur_pos);
      PVector acc_mv = motion_field.get(idx);
      motion_field.set(
          idx, new PVector(acc_mv.x + mv.x, acc_mv.y + mv.y, acc_mv.z + 1));
    }
    for (int i = 0; i < r_num * c_num; i++) {
      PVector mv = motion_field.get(i);
      if (mv.z > 0) {
        motion_field.set(i, new PVector(mv.x / mv.z, mv.y / mv.z, 0));
      }
    }
    return motion_field;
  }

  void showMotionField() {
    ortho(-width, 0, -height, 0);
    camera(0, 0, 0, 0, 0, 1, 0, 1, 0);
    getMotionField();
    int r_num = height / block_size, c_num = width / block_size;
    for (int i = 0; i < r_num; i++)
      for (int j = 0; j < c_num; j++) {
        PVector mv = motion_field.get(i * c_num + j);
        float ox = j * block_size + 0.5f * block_size;
        float oy = i * block_size + 0.5f * block_size;
        stroke(255, 0, 0);
        line(ox, oy, ox + mv.x, oy + mv.y);
      }
  }
}
