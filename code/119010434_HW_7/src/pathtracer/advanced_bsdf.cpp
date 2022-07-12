#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Assignment 7: Part 1
  // Implement MirrorBSDF
  return Vector3D();
}

void MirrorBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Mirror BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    ImGui::TreePop();
  }
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D h) {
  // TODO Assignment 7: Part 2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
  
  
  
  double beckmann_distribution;
  // double roughness = 0.005;
  // roughness = 0.05;
  // roughness = 0.1;
  // roughness = 0.3;
  // alpha is roughness number.
  double d_up = exp(-pow(tan(acos(dot(Vector3D(0., 0., 1.), h)))/ alpha, 2));
  double d_down = PI * pow(alpha, 2) * pow(cos(acos(dot(Vector3D(0., 0., 1.), h))), 4);
  beckmann_distribution = d_up / d_down;

  return beckmann_distribution;
}

Vector3D MicrofacetBSDF::F(const Vector3D wi) {
  // TODO Assignment 7: Part 2
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Vector3D.
  
  
  Vector3D rs, rp;
  double cos_i = cos_theta(wi);
  rs = (eta*eta + k*k - 2.*eta*cos_i + cos_i*cos_i) / (eta*eta + k*k + 2.*eta*cos_i + cos_i*cos_i);
  rp = ((eta*eta + k*k)*cos_i*cos_i - 2.*eta*cos_i + 1.) / ((eta*eta + k*k)*cos_i*cos_i + 2.*eta*cos_i + 1);
  return (rs+rp) / 2.;
  return Vector3D();
}

Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
  // TODO Assignment 7: Part 2
  // Implement microfacet model here.
  if (wo.z <= 0 || wi.z <= 0)
  {
    return Vector3D(0.);
  }
  
  Vector3D f_result, h;
  h = wo + wi;
  h.normalize();
  f_result = (F(wi) * G(wo, wi) * D(h)) / (4 * dot(Vector3D(0, 0, 1), wo) * dot(Vector3D(0,0,1), wi));
  return f_result;

  return Vector3D();
}

Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Assignment 7: Part 2
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.
  
  
  
  Vector2D uni_random_num = sampler.get_sample();
  double r1, r2;
  r1 = uni_random_num.x;
  r2 = uni_random_num.y;
  double theta_h = atan(pow(-alpha*alpha*log(1-r1), 0.5));
  double phi_h = 2.0 * PI * r2;
  double p_theta = ((2*sin(theta_h)) / (alpha*alpha*pow(cos(theta_h), 3))) * exp(-pow(tan(theta_h)/alpha, 2));
  double p_phi = 1/(2.0*PI);

  Vector3D h;
  h.x = sin(theta_h) * cos(phi_h);
  h.y = sin(theta_h) * sin(phi_h);
  h.z = cos(theta_h);

  // h = wo + wi;
  // *wi = h.norm() - wo.norm();
  // wi->normalize();
  *wi = 2 * (dot(h, wo)) * h - wo;
  // *wi = cosineHemisphereSampler.get_sample(pdf);
  if ((*wi).z <= 0)
  {
    *pdf = 0.;
    return Vector3D(0.);
  }
  double p_w_h = p_theta * p_phi / sin(theta_h);
  *pdf = p_w_h/(4.0 * dot(*wi, h)); 
  
  return MicrofacetBSDF::f(wo, *wi);

}

void MicrofacetBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Micofacet BSDF"))
  {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Assignment 7: Part 1
  // Implement RefractionBSDF
  return Vector3D();
}

void RefractionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Assignment 7: Part 1
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305
  return Vector3D();
}

void GlassBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

  // TODO Assignment 7: Part 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.


}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

  // TODO Assignment 7: Part 1
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  return true;

}

} // namespace CGL
