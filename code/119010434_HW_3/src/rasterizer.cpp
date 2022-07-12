#include "rasterizer.h"
#include <vector>


using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)


    sample_buffer[y * width + x] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  bool pointInSideTriangle(float i, float j, float x0, float y0, float x1, float y1, float x2, float y2)
  {
    Vector3D v0, v1, v2;
    v0.x = x0; v1.x = x1; v2.x = x2;
    v0.y = y0; v1.y = y1; v2.y = y2;
    v0.z = 0;  v1.z = 0;  v2.z = 0;
    Vector3D v_draw;
    v_draw.x = i + 0.5; v_draw.y = j + 0.5; v_draw.z = 0;
    Vector3D ans_1 = CGL::cross((v1-v0),(v_draw-v0));
    Vector3D ans_2 = CGL::cross((v2-v1),(v_draw-v1));
    Vector3D ans_3 = CGL::cross((v0-v2),(v_draw-v2));
    float ans1 = ans_1[2];
    float ans2 = ans_2[2];
    float ans3 = ans_3[2];
    if (ans1 >= 0 && ans2 >= 0 && ans3 >= 0)
    {
      return true;
    }
    else if (ans1 <= 0 && ans2 <= 0 && ans3 <= 0)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  float find_max(float x0, float x1, float x2)
  {
    if (x0 >= x1 && x0 >= x2)
    {
      return x0;
    }
    else if (x1 >= x0 && x1 >=x2)
    {
      return x1;
    }
    else 
    {
      return x2;
    }
  }

  float find_min(float x0, float x1, float x2)
  {
    if (x0 <= x1 && x0 <= x2)
    {
      return x0;
    }
    else if (x1 <= x0 && x1 <=x2)
    {
      return x1;
    }
    else 
    {
      return x2;
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    cout<<"ok"<<endl;
    cout<<"x0:"<<x0<<"y0:"<<y0<<"x1:"<<x1<<"y1:"<<y1<<"x2:"<<x2<<"y2:"<<y2<<endl;
    float max_x = max(x0, x1);
    // max_x = ceil(max(max_x, x2));
    max_x = find_max(x0, x1, x2);
    max_x = ceil(max_x);
    float min_x = min(x0, x1);
    // min_x = floor(min(max_x, x2));
    min_x = find_min(x0, x1, x2);
    min_x = floor(min_x);

    float max_y = max(y0, y1);
    // max_y = ceil(max(max_y, y2));
    max_y = find_max(y0, y1, y2);
    max_y = ceil(max_y);
    float min_y = min(y0, y1);
    // min_y = floor(min(min_y, y2));
    min_y = find_min(y0, y1, y2);
    min_y = floor(min_y);

    //min_x = 0;
    
    // cout<<"max.x:"<<max_x<<"max.y"<<max_y<<endl;
    // cout<<"min.x:"<<min_x<<"min.y"<<min_y<<endl;


    bool inSideTriangle = false;
    // this->sample_rate = 9;
    // cout<<"sr:"<<this->sample_rate<<endl;
    if (this->sample_rate == 1)
    {
      for (float i = min_x; i <= max_x; i++)
      {
        for (float j = min_y; j <= max_y; j++)
        {
          inSideTriangle = pointInSideTriangle(i, j, x0, y0, x1, y1, x2, y2);
          
          if (inSideTriangle)
          {
            fill_pixel((int)i, (int)j, color);
          }   
        }
      }
    }
    else
    {
      // cout<<"nihao"<<endl;
      // TODO: Task 2: Update to implement super-sampled rasterization
      // ssaa/msaa
      int sample_rate_now = this->sample_rate;
      Color new_color;
      float cover_point_num;

      Color color_list[sample_rate_now];
      Color color_white;
      color_white.r = 1;
      color_white.g = 1;
      color_white.b = 1;
      Color color_black;
      color_black.r = 0;
      color_black.g = 0;
      color_black.b = 0;
      
      

      for (float i = min_x; i <= max_x; i++)
      {
        for (float j = min_y; j < max_y; j++)
        {
          cover_point_num = 0.0;
          for (int k = 0; k < sample_rate_now; k++)
          {
            float x_frag, y_frag;
            float frag_num;
            frag_num = 1 + sqrt((float)sample_rate_now);
            float unit_frag = 1/frag_num;

            // cout<<"unit_frag: "<<unit_frag<<endl;
            
            int space_num;
            space_num = (int)sqrt((float)sample_rate_now);
            x_frag = k/space_num + 1;
            y_frag = k%space_num + 1;
            // cout<<"k:"<<k<<" space_num:"<<space_num<<" xfrag:"<<x_frag<<" yfrag:"<<y_frag<<endl;
            // cout<<"point_pos: ("<<x_frag*unit_frag<<","<<y_frag*unit_frag<<")"<<endl;
            inSideTriangle = pointInSideTriangle(i+(x_frag*unit_frag), j+(y_frag*unit_frag), x0, y0, x1, y1, x2, y2);

            if (inSideTriangle)
            {
              color_list[k] = color;
            }
            else
            {
              color_list[k] = color_white;
            }
          }
          for (int list_num = 0; list_num < sample_rate_now; list_num++)
          {
            /* code */
            // cout<<"color part:"<<color_list[list_num]<<endl;
            // cout<<new_color<<endl;
            new_color = new_color + color_list[list_num] * (1.0/sample_rate_now);
          }
          
          // cout<<"covernum: "<<cover_point_num<<endl;
          // cout<<"old color: "<<color<<endl;
          // cout<<"scale: "<<cover_point_num/sample_rate_now<<endl;
          // cout<<new_color<<endl;
          // if (new_color.r+new_color.g+new_color.b > 0)
          
          /* code */
          fill_pixel((int)i, (int)j, new_color);
          // cout<<"covernum: "<<cover_point_num<<endl;
          
          new_color = color_black;
        }
      }
    }
  }

  float Get_Cood_0(float i, float j, float x0, float y0, float x1, float y1, float x2, float y2);
  float Get_Cood_1(float i, float j, float x0, float y0, float x1, float y1, float x2, float y2);
  float Get_Cood_2(float i, float j, float x0, float y0, float x1, float y1, float x2, float y2);


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle

    float max_x = max(x0, x1);
    // max_x = ceil(max(max_x, x2));
    max_x = find_max(x0, x1, x2);
    max_x = ceil(max_x);
    float min_x = min(x0, x1);
    // min_x = floor(min(max_x, x2));
    min_x = find_min(x0, x1, x2);
    min_x = floor(min_x);

    float max_y = max(y0, y1);
    // max_y = ceil(max(max_y, y2));
    max_y = find_max(y0, y1, y2);
    max_y = ceil(max_y);
    float min_y = min(y0, y1);
    // min_y = floor(min(min_y, y2));
    min_y = find_min(y0, y1, y2);
    min_y = floor(min_y);

    bool inSideTriangle = false;

    float cood_0, cood_1, cood_2;
    cood_0 = 0.;
    cood_1 = 0.;
    cood_2 = 0.;

    for (float i = min_x; i <= max_x; i++)
    {
      for (float j = min_y; j <= max_y; j++)
      {
        inSideTriangle = pointInSideTriangle(i+0.5, j+0.5, x0, y0, x1, y1, x2, y2);
        
        if (inSideTriangle)
        {
          //P = iA(0) + jB(1) + kC(2)
          //i = -(Px-Bx)(Cy-By)+(Py-By)(Cy-By)/-ABCBABCB
          //j = -(Px-Cx)(Ay-Cy)+(Py-Cy)(Ax-Cx)/-BCACBCAC
          // cood_0 = (-(i+0.5-x1)*(y2-y1)+(j+0.5-y1)*(x2-x1))/(-(x0-x1)*(y2-y1)+(y0-y1)*(y2-y1));
          // cood_1 = (-(i+0.5-x2)*(y0-y2)+(j+0.5-y2)*(x0-x2))/(-(x1-x2)*(y0-y2)+(y1-y2)*(x0-x2));
          // cood_2 = (-(i+0.5-x0)*(y1-y0)+(j+0.5-y0)*(x1-x0))/(-(x2-x0)*(y1-y0)+(y2-y0)*(x1-x0));
          cood_0 = Get_Cood_0(i, j, x0, y0, x1, y1, x2, y2);
          cood_1 = Get_Cood_1(i, j, x0, y0, x1, y1, x2, y2);
          cood_2 = Get_Cood_2(i, j, x0, y0, x1, y1, x2, y2);
          // cood_2 = 1-cood_0-cood_1;
          // Color color = cood_0 * c0;
          // Color color = cood_1 * c1;
          // Color color = cood_2 * c2;

          // Color color = cood_0 * c0 + cood_1 * c1 + cood_2 * c2;
          Color color = cood_1 * c1 + cood_2 * c2;
          // Color color = cood_0 * c0 + cood_1 * c1;
          // Color color = c2;
          fill_pixel((int)i, (int)j, color);
        }   
      }
    }
  }

  float Get_Cood_0(float i, float j, float x0, float y0, float x1, float y1, float x2, float y2)
  {
    return (-(i-x1)*(y2-y1)+(j-y1)*(x2-x1))/(-(x0-x1)*(y2-y1)+(y0-y1)*(y2-y1));
  }

  float Get_Cood_1(float i, float j, float x0, float y0, float x1, float y1, float x2, float y2)
  {
    return (-(i-x2)*(y0-y2)+(j-y2)*(x0-x2))/(-(x1-x2)*(y0-y2)+(y1-y2)*(x0-x2));
  }

  float Get_Cood_2(float i, float j, float x0, float y0, float x1, float y1, float x2, float y2)
  {
    return (-(i-x0)*(y1-y0)+(j-y0)*(x1-x0))/(-(x2-x0)*(y1-y0)+(y2-y0)*(x1-x0));
  }

  Vector3D barycentric(Vector2D A, Vector2D B, Vector2D C, Vector2D P) {
    Vector3D s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    Vector3D u = cross(s[0], s[1]);
    if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return Vector3D(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return Vector3D(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
  }

  void test_cood(float i=1, float j=1, float x0=0, float y0=0, float x1=2, float y1=0, float x2=0, float y2=2)
  {
    float cood_0, cood_1, cood_2;
    //cood_0 = (-(i+0.5-x1)*(y2-y1)+(j+0.5-y1)*(x2-x1))/(-(x0-x1)*(y2-y1)+(y0-y1)*(y2-y1));
    cood_0 = Get_Cood_0(i, j, x0, y0, x1, y1, x2, y2);
    //cood_1 = (-(i+0.5-x2)*(y0-y2)+(j+0.5-y2)*(x0-x2))/(-(x1-x2)*(y0-y2)+(y1-y2)*(x0-x2));
    cood_1 = Get_Cood_1(i, j, x0, y0, x1, y1, x2, y2);
    // cood_2 = 1-cood_0-cood_1;
    cood_2 = Get_Cood_2(i, j, x0, y0, x1, y1, x2, y2);
    cout<<"test: cood0:"<<cood_0<<" cood1:"<<cood_1<<" cood2:"<<cood_2<<endl;
  }

  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
    cout<<"test1"<<endl;
    // test_cood();

    // u0 = u0 * tex.width;
    // v0 = v0 * tex.height;
    // u1 = u1 * tex.width;
    // v1 = v1 * tex.height;
    // u1 = u2 * tex.width;
    // v1 = v2 * tex.height;
    cout<<" u0:"<<u0<<" v0:"<<v0<<" u1:"<<u1<<" v1:"<<v1<<" u2:"<<u2<<" v2:"<<v2<<endl;
    cout<<"width:"<<tex.width<<" heigth:"<<tex.width<<endl;
    cout<<"test2"<<endl;

    float max_x = max(x0, x1);
    // max_x = ceil(max(max_x, x2));
    max_x = find_max(x0, x1, x2);
    max_x = ceil(max_x);
    float min_x = min(x0, x1);
    // min_x = floor(min(max_x, x2));
    min_x = find_min(x0, x1, x2);
    min_x = floor(min_x);

    float max_y = max(y0, y1);
    // max_y = ceil(max(max_y, y2));
    max_y = find_max(y0, y1, y2);
    max_y = ceil(max_y);
    float min_y = min(y0, y1);
    // min_y = floor(min(min_y, y2));
    min_y = find_min(y0, y1, y2);
    min_y = floor(min_y);
    cout<<"test3"<<endl;

    bool inSideTriangle = false;

    float cood_0, cood_1, cood_2;
    cood_0 = 0.;
    cood_1 = 0.;
    cood_2 = 0.;
    float cood_0_dx, cood_1_dx, cood_2_dx;
    float cood_0_dy, cood_1_dy, cood_2_dy;

    float u_texture, v_texture;
    Color color_get;
    cout<<"test4"<<endl;

    for (float i = min_x; i <= max_x; i++)
    {
      for (float j = min_y; j <= max_y; j++)
      {
        inSideTriangle = pointInSideTriangle(i+0.5, j+0.5, x0, y0, x1, y1, x2, y2);
        
        if (inSideTriangle)
        {
          //cood_0 = (-(i+0.5-x1)*(y2-y1)+(j+0.5-y1)*(x2-x1))/(-(x0-x1)*(y2-y1)+(y0-y1)*(y2-y1));
          cood_0 = Get_Cood_0(i, j, x0, y0, x1, y1, x2, y2);
          //cood_1 = (-(i+0.5-x2)*(y0-y2)+(j+0.5-y2)*(x0-x2))/(-(x1-x2)*(y0-y2)+(y1-y2)*(x0-x2));
          cood_1 = Get_Cood_1(i, j, x0, y0, x1, y1, x2, y2);
          //cood_2 = 1-cood_0-cood_1;
          cood_2 = Get_Cood_2(i, j, x0, y0, x1, y1, x2, y2);
          cout<<"cood0:"<<cood_0<<" cood1:"<<cood_1<<" cood2:"<<cood_2<<endl;
          cout<<"puls:"<<cood_0 + cood_1 + cood_2<<endl;
          Vector2D A, B, C, P;
          A.x = x0; A.y = y0; B.x = x1; B.y = y1; C.x = x2; C.y = y2;
          P.x = i; P.y = j;
          cood_0 = barycentric(A,B,C,P).x;
          cood_1 = barycentric(A,B,C,P).y;
          cood_2 = barycentric(A,B,C,P).z;

          if (this->lsm == L_NEAREST || this->lsm == L_LINEAR)
          {
            cout<<"test???"<<endl;
            Vector2D PX, PY;
            PX.x = i+1;
            PX.y = j;
            // cood_0_dx = Get_Cood_0(i+1, j, x0, y2, x1, y1, x2, y2);
            // cood_1_dx = Get_Cood_1(i+1, j, x0, y2, x1, y1, x2, y2);
            // cood_2_dx = Get_Cood_2(i+1, j, x0, y2, x1, y1, x2, y2);
            cood_0_dx = barycentric(A,B,C,PX).x;
            cood_1_dx = barycentric(A,B,C,PX).y;
            cood_2_dx = barycentric(A,B,C,PX).z;

            PY.x = i;
            PY.y = j+1;
            // cood_0_dy = Get_Cood_0(i, j+1, x0, y2, x1, y1, x2, y2);
            // cood_1_dy = Get_Cood_1(i, j+1, x0, y2, x1, y1, x2, y2);
            // cood_2_dy = Get_Cood_2(i, j+1, x0, y2, x1, y1, x2, y2);
            cood_0_dy = barycentric(A,B,C,PY).x;
            cood_1_dy = barycentric(A,B,C,PY).y;
            cood_2_dy = barycentric(A,B,C,PY).z;

            SampleParams Params_send;
            Params_send.p_uv.x = cood_0 * u0 + cood_1 * u1 + cood_2 * u2;
            Params_send.p_uv.y = cood_0 * v0 + cood_1 * v1 + cood_2 * v2;

            Params_send.p_dx_uv.x = cood_0_dx*u0 + cood_1_dx*u1 + cood_2_dx*u2;
            Params_send.p_dx_uv.y = cood_0_dx*v0 + cood_1_dx*v1 + cood_2_dx*v2;

            Params_send.p_dy_uv.x = cood_0_dy*u0 + cood_1_dy*u1 + cood_2_dy*u2;
            Params_send.p_dy_uv.y = cood_0_dy*v0 + cood_1_dy*v1 + cood_2_dy*v2;

            Params_send.lsm = this->lsm;
            Params_send.psm = this->psm;

            //float level_num = tex.get_level(Params_send);

            //int level_num_int = round(level_num);

            color_get = tex.sample(Params_send);
            fill_pixel((int)i, (int)j, color_get);
          }
          else
          {
            u_texture = cood_0 * u0 + cood_1 * u1 + cood_2 * u2;
            cout<<"cood0:"<<cood_0<<" cood1:"<<cood_1<<" cood2:"<<cood_2<<endl;
            cout<<"puls:"<<cood_0 + cood_1 + cood_2<<endl;
            v_texture = cood_0 * v0 + cood_1 * v1 + cood_2 * v2;
            Vector2D uv_get;
            uv_get.x = u_texture;
            uv_get.y = v_texture;
            cout<<"test???"<<endl;
            if (this->psm == P_NEAREST)
            {
              /* code */
              cout<<"test1"<<endl;
              color_get = tex.sample_nearest(uv_get, 0);
            }
            else if (this->psm == P_LINEAR)
            {
              /* code */
              cout<<"test2"<<endl;
              color_get = tex.sample_bilinear(uv_get, 0);
            }
            else {
              cout<<"something wrong!!"<<endl;
              color_get = tex.sample_nearest(uv_get, 0);
            }
            cout<<"test!!!"<<endl;
            fill_pixel((int)i, (int)j, color_get);
          }
        }   
      }
    }
  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support


    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = sample_buffer[y * width + x];

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
