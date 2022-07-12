#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

  // TODO Assignment 7: Part 4
  // compute position and direction of ray=from the input sensor sample coordinate.
  // Note: use rndR and rndTheta to uniformly sample a unit disk.
  x = 2 * (x - 0.5) * tan(radians(hFov) * 0.5);
  y = 2 * (y - 0.5) * tan(radians(vFov) * 0.5);
  double pLens_x, pLens_y;
  pLens_x =lensRadius * pow(rndR, 0.5) * cos(rndTheta);
  pLens_x =lensRadius * pow(rndR, 0.5) * sin(rndTheta);
  Vector3D pLens(pLens_x, pLens_y, 0.0);
  Vector3D pFocus(x*focalDistance, y*focalDistance, -focalDistance);
  Vector3D ray_dir;
  ray_dir = pFocus - pLens;
  ray_dir /= ray_dir.norm();
  
  Ray r(c2w * pLens + pos, c2w * ray_dir);
  r.min_t = nClip;
  r.max_t = fClip;
  return r;
}


} // namespace CGL
