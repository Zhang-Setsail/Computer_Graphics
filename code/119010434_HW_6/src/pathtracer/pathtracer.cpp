#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"


using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
  // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Vector3D L_out;

  // TODO (Part 3): Write your sampling loop here
  // TODO BEFORE YOU BEGIN
  // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading 
  Vector3D wi, nextDir, nextOri;
  // Ray next_ray;
  double pdf;
  double wi_dot;
  Vector3D normal_vector_obj_spa(0, 0, 1);
  Intersection next_insec;
  Vector3D light_next_insec;
  for (int i = 0; i < num_samples; i++)
  {
    Vector3D f;
    f = isect.bsdf->sample_f(w_out, &wi, &pdf);
    wi_dot = dot(wi, normal_vector_obj_spa);
    nextDir = o2w * wi;
    nextOri = hit_p;
    Ray next_ray = Ray(nextOri, nextDir, INF_D, 1);
    next_ray.min_t = EPS_F;
    // an epsilon constant that can be used to avoid numerical precision issues.
    if (wi_dot>0 && bvh->intersect(next_ray, &next_insec))
    {
      light_next_insec = next_insec.bsdf->get_emission();
      L_out = L_out + light_next_insec*f*wi_dot/pdf;
    }
  }
  L_out = L_out/num_samples;
  return L_out;
  
  return Vector3D(1.0);

}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out;

  // scene->lights is a vector storing pointers of all the lights in the scene
  //遍历所有光线
  Vector3D wi, nextDir, nextOri;
  // Ray next_ray;
  double pdf;
  double wi_dot;
  double distance;
  Intersection next_insec;
  Vector3D emitted_radiance;
  Vector3D out_total_light;
  Vector3D out_single_light;
  int sample_times = 0;
  for (auto light = scene->lights.begin(); light != scene->lights.end(); light++)
  {
    for (int i = 0; i < ns_area_light; i++)
    {
      sample_times = sample_times + 1;
      emitted_radiance = (*light)->sample_L(hit_p, &wi, &distance, &pdf);
      wi_dot = dot(wi, isect.n);
      nextDir = wi;
      nextOri = hit_p;
      Ray next_ray = Ray(nextOri, nextDir, distance-EPS_F, 1);
      next_ray.min_t = EPS_F;
      if (wi_dot>0 && !(bvh->intersect(next_ray, &next_insec)))
      {
        Vector3D f;
        f = isect.bsdf->f(w_out, wi);
        out_single_light = emitted_radiance*f*wi_dot/pdf;
        out_total_light = out_total_light + out_single_light;
      }
      //当碰到光源时停止is_delta_light()
      if ((*light)->is_delta_light())
      {
        break;
      }
    }
    L_out = L_out + out_total_light/sample_times;
  }
  return L_out;

  return Vector3D(1.0);

}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light

  return isect.bsdf->get_emission();
  //只取光源的光
}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`

  if (direct_hemisphere_sample)
  {
    return estimate_direct_lighting_hemisphere(r, isect);
  }
  return estimate_direct_lighting_importance(r, isect);
  
  return Vector3D(1.0);
}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D L_out(0, 0, 0);

  // TODO: Part 4, Task 2
  // Returns the one bounce radiance + radiance from extra bounces at this point.
  // Should be called recursively to simulate extra bounces.

  int light_num = scene->lights.size();
  int sample_num = light_num * ns_area_light;
  L_out = L_out + one_bounce_radiance(r, isect);

  // In theory, the probability of terminating can be arbitrary -- but we suggest using a termination probability between 0.3 or 0.4.
  double terminate_prob = 0.35;
  Vector3D wi;
  double pdf;
  Vector3D f;
  double wi_dot;
  Vector3D nextDir, nextOri;
  Vector3D normal_vector_obj_spa(0, 0, 1);
  Intersection next_insec;
  Vector3D light_single;
  // 在终止之前
  if (coin_flip(1.-terminate_prob))
  {
    f = isect.bsdf->sample_f(w_out, &wi, &pdf);
    wi_dot = dot(wi, normal_vector_obj_spa);
    nextDir = o2w * wi;
    nextOri = hit_p;
    Ray next_ray = Ray(nextOri, nextDir, INF_D, r.depth - 1);
    next_ray.min_t = EPS_F;
    if (next_ray.depth > 1)
    {
      if (wi_dot)
      {
        if (bvh->intersect(next_ray, &next_insec))
        {
          light_single = at_least_one_bounce_radiance(next_ray, next_insec);
          L_out = L_out + light_single*f* wi_dot/pdf/(1-terminate_prob);
        } 
      }
    }
  }
  return L_out;
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.
  
  if (!bvh->intersect(r, &isect)){
    // std::cout<<!bvh->intersect(r, &isect)<<endl;
    return envLight ? envLight->sample_dir(r) : L_out;
  }


  // L_out = (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);

  // TODO (Part 3): Return the direct illumination.
  L_out = zero_bounce_radiance(r, isect);
  // L_out += one_bounce_radiance(r, isect);//一次弹射

  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct
  L_out += at_least_one_bounce_radiance(r, isect);//多次弹射

  return L_out;
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"

  int num_cam_ray_per_pixel = ns_aa;
  Vector2D origin = Vector2D(x, y);
  Vector3D radiance_total(0.);
  Vector3D scene_radiance;
  int radiance_num = 0;
  float s1 = 0, s2 = 0;
  float sigma_2;
  float u;
  float I;
  int already_sample_num = 0;
  for (int i = 0; i < num_cam_ray_per_pixel; i++)
  {
    bool need_adapt = (already_sample_num == samplesPerBatch);
    if (need_adapt)
    {
      already_sample_num = 0;
      u = s1 / (float)radiance_num;
      sigma_2 = (1/(float)radiance_num-1) * (s2-(s1*s1)/(float)radiance_num);
      I = 1.96*pow(sigma_2, 0.5)/pow((double)radiance_num, 0.5);
      bool pixel_convergence = (I <= maxTolerance*u);
      if (pixel_convergence)
      {
        break;
      } 
    }
    Vector2D pixel_sample = origin + gridSampler->get_sample();
    double pos_x, pos_y;
    pos_x = pixel_sample.x / sampleBuffer.w;
    pos_y = pixel_sample.y / sampleBuffer.h;
    Ray world_ray = camera->generate_ray(pos_x, pos_y);
    world_ray.depth = max_ray_depth;
    scene_radiance = est_radiance_global_illumination(world_ray);

    float sample_illuminance = scene_radiance.illum();
    s1 = s1 + sample_illuminance;
    s2 = s2 + pow(sample_illuminance, 2);
    radiance_total = radiance_total + scene_radiance;
    already_sample_num = already_sample_num + 1;
    radiance_num = radiance_num + 1;
  }
  radiance_total = radiance_total / radiance_num;
  sampleBuffer.update_pixel(radiance_total, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = radiance_num;
  
}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
