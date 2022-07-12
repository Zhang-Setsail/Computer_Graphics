#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {

  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.


  BBox bbox;

  for (auto p = start; p != end; p++) {
    BBox bb = (*p)->get_bbox();
    bbox.expand(bb);
  }

  BVHNode *node = new BVHNode(bbox);

  int number_primitives = (end - start);

  if (number_primitives <= max_leaf_size)
  {
    node->start = start;
    node->end = end;
    node->l = NULL;
    node->r = NULL;
  }
  else {
    Vector3D avg(0.);
    for (auto p = start; p != end; p++) {
      avg = avg + (*p)->get_bbox().centroid();
    }
    avg = avg/(float)number_primitives;

    int split_axis = 0; //使用哪条轴切分
    float min_BBox_face = MAXFLOAT; //最小加权（乘物体数量）面积切分

    for (int axis = 0; axis < 3; axis++)
    {
      std::vector<Primitive*> new_left_primitives;
      int new_left_num = 0;
      std::vector<Primitive*> new_right_primitives;
      int new_right_num = 0;

      for (auto p = start; p != end; p++) {
        bool primitive_to_left = (*p)->get_bbox().centroid()[axis] <= avg[axis];
        if (primitive_to_left) {
          new_left_primitives.push_back(*p);
          new_left_num = new_left_num + 1;
        }
        else {
          new_right_primitives.push_back(*p);
          new_right_num = new_right_num + 1;
        }
      }
      //把所有物体根据轴与平均距离分开
      BBox new_left_bbox;
      BBox new_right_bbox;
      for (int i = 0; i < new_left_num; i++)
      {
        BBox bb = new_left_primitives[i]->get_bbox();
        new_left_bbox.expand(bb);
      }
      for (int i = 0; i < new_right_num; i++)
      {
        BBox bb = new_right_primitives[i]->get_bbox();
        new_right_bbox.expand(bb);
      }
      float left_and_right_face;
      float new_left_x = new_left_bbox.extent.x;
      float new_left_y = new_left_bbox.extent.y;
      float new_left_z = new_left_bbox.extent.z;
      float new_right_x = new_right_bbox.extent.x;
      float new_right_y = new_right_bbox.extent.y;
      float new_right_z = new_right_bbox.extent.z;

      left_and_right_face = new_left_num * (new_left_x*new_left_y + new_left_x*new_left_z + new_left_y*new_left_z);
      left_and_right_face = left_and_right_face + new_right_num * (new_right_x*new_right_y + new_right_x*new_right_z + new_right_y*new_right_z);

      if (left_and_right_face < min_BBox_face)
      {
        min_BBox_face = left_and_right_face;
        split_axis = axis;
      }
    }
    std::vector<Primitive*> real_left_primitives;
    int left_primitives_num = 0;
    std::vector<Primitive*> real_right_primitives;
    int right_primitives_num = 0;
    //已经获得要切分的轴，上面是实际切分的轴分开的物体
    for (auto p = start; p != end; p++) {
      bool primitive_to_left = (*p)->get_bbox().centroid()[split_axis] <= avg[split_axis];
      if (primitive_to_left) {
        real_left_primitives.push_back(*p);
        left_primitives_num = left_primitives_num + 1;
      }
      else {
        real_right_primitives.push_back(*p);
        right_primitives_num = right_primitives_num + 1;
      }
    }
    std::vector<Primitive *>::iterator center, center_right;
    center = start;
    center_right = start;
    for (int i = 0; i < number_primitives; i++)
    {
      if (i < left_primitives_num)
      {
        *center = real_left_primitives[i];
        center++;
        center_right++;
      }
      else if (i >= left_primitives_num)
      {
        *center_right = real_right_primitives[i-left_primitives_num];
        center_right++;
      }
    }
    node->l = construct_bvh(start, center, max_leaf_size);
    node->r = construct_bvh(center, end, max_leaf_size);
  }
  return node;

  
}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.
  
  
  total_isects++;

  if (!node->bb.intersect(ray, ray.min_t, ray.max_t))
  {
    return false;
  }
  if (!node->isLeaf())
  {
    return has_intersection(ray, node->l) || has_intersection (ray, node->r);
    //如果是root节点需要查看两个叶节点
  }
  // for (auto p : primitives) {
  for (auto leaf_primitives = node->start; leaf_primitives != node->end; leaf_primitives++)
  {
    // if (p->has_intersection(ray))
    if ((*leaf_primitives)->has_intersection(ray))
    {
      return true;
    }
  }
  return false;

  // for (auto p : primitives) {
  //   total_isects++;
  //   if (p->has_intersection(ray))
  //     return true;
  // }
  // return false;

}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  total_isects++;

  if (!(node->bb.intersect(ray, ray.min_t, ray.max_t)))
  {
    return false;
  }
  if (!node->isLeaf())
  {
    bool left_son_return = intersect(ray, i, node->l);
    bool right_son_return = intersect (ray, i, node->r);
    //zhe liu yin wei intersect de yuan yin ,yi ding yao dan du yin yong mei ge han shu;
    bool sum_insecter_return = left_son_return || right_son_return;
    return sum_insecter_return;
    //如果是root节点需要查看两个叶节点
  }
  bool hit = false;
  // for (auto p : primitives) {
  for (auto leaf_primitives = node->start; leaf_primitives != node->end; leaf_primitives++)
  {
    // if (p->has_intersection(ray))
    hit = hit || (*leaf_primitives)->intersect(ray, i);
  }
  return hit;
    

  // bool hit = false;
  // for (auto p : primitives) {
  //   total_isects++;
  //   hit = p->intersect(ray, i) || hit;
  // }
  // return hit;


}

} // namespace SceneObjects
} // namespace CGL
