#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.

  Vector3D intersect_test(0.);
  double test_min_t_x, test_min_t_y, test_min_t_z, min_t;
  double test_max_t_x, test_max_t_y, test_max_t_z, max_t;
  double test_mid;
  //x测试
  intersect_test.x = 1.;
  test_min_t_x = dot((min-r.o), intersect_test)/dot(r.d, intersect_test);
  test_max_t_x = dot((max-r.o), intersect_test)/dot(r.d, intersect_test);
  if (test_max_t_x < test_min_t_x)
  {
    test_mid = test_max_t_x;
    test_max_t_x = test_min_t_x;
    test_min_t_x = test_mid;
  }
  intersect_test.x = 0.;
  //y测试
  intersect_test.y = 1.;
  test_min_t_y = dot((min-r.o), intersect_test)/dot(r.d, intersect_test);
  test_max_t_y = dot((max-r.o), intersect_test)/dot(r.d, intersect_test);
  if (test_max_t_y < test_min_t_y)
  {
    test_mid = test_max_t_y;
    test_max_t_y = test_min_t_y;
    test_min_t_y = test_mid;
  }
  intersect_test.y = 0.;
  //z测试
  intersect_test.z = 1.;
  test_min_t_z = dot((min-r.o), intersect_test)/dot(r.d, intersect_test);
  test_max_t_z = dot((max-r.o), intersect_test)/dot(r.d, intersect_test);
  if (test_max_t_z < test_min_t_z)
  {
    test_mid = test_max_t_z;
    test_max_t_z = test_min_t_z;
    test_min_t_z = test_mid;
  }
  intersect_test.z = 0.;

  min_t = std::max(test_min_t_x, test_min_t_y);
  max_t = std::min(test_max_t_x, test_max_t_y);
  if (min_t > max_t)
  {
    return false;
  }
  if (test_min_t_z > max_t || test_max_t_z < min_t)
  {
    return false;
  }
  min_t = std::max(test_min_t_z, min_t);
  max_t = std::max(test_max_t_z, max_t);
  if (min_t > max_t)
  {
    return false;
  }
  bool final_overlap_test = t1>=min_t && t0<=max_t;
  return final_overlap_test;

}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
