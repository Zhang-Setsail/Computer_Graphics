#include "sphere.h"

#include <cmath>

#include "pathtracer/bsdf.h"
#include "util/sphere_drawing.h"

namespace CGL {
namespace SceneObjects {

bool Sphere::test(const Ray &r, double &t1, double &t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.


  return true;

}

bool Sphere::has_intersection(const Ray &r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double a = (r.d).norm2();
  double b = 2 * dot((r.o - o), r.d);
  double c = (r.o - o).norm2() - r2;

  double b_sqar_4ac = b*b-4*a*c;
  if (b_sqar_4ac >= 0)
  {
    double t_large = (-b+pow(b_sqar_4ac, 0.5))/(2*a);
    double t_small = (-b-pow(b_sqar_4ac, 0.5))/(2*a);
    if ((t_large>=r.min_t && t_large <= r.max_t) || (t_small >= r.min_t && t_small <= r.max_t))
    {
      return true;
    }
    return false;
  }
  return false;
}

bool Sphere::intersect(const Ray &r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  double a = (r.d).norm2();
  double b = 2 * dot((r.o - o), r.d);
  double c = (r.o - o).norm2() - r2;
  double b_sqar_4ac = b*b-4*a*c;

  if (b_sqar_4ac >= 0)
  {
    double t_large = (-b+pow(b_sqar_4ac, 0.5))/(2*a);
    double t_small = (-b-pow(b_sqar_4ac, 0.5))/(2*a);
    if (t_small >= r.min_t && t_small <= r.max_t)
    {
      r.max_t = t_small;
      i->t = t_small;
      i->primitive = this;
      i->n = r.o - o + t_small * r.d;
      i->n.normalize();
      i->bsdf = get_bsdf();
      return true;
    }
    else if (t_large >= r.min_t && t_large <= r.max_t)
    {
      r.max_t = t_large;
      i->t = t_large;
      i->primitive = this;
      i->n = r.o - o + t_large * r.d;
      i->n.normalize();
      i->bsdf = get_bsdf();
      return true;
    }
  }
  return false;
}

void Sphere::draw(const Color &c, float alpha) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color &c, float alpha) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

} // namespace SceneObjects
} // namespace CGL
