#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

  Color Texture::sample(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
    float level_num = get_level(sp);

    int level_num_int = round(level_num);
    if (sp.lsm == L_NEAREST)
    {
      if (sp.psm == P_NEAREST)
      {
        return sample_nearest(sp.p_uv, level_num_int);
      }
      else if (sp.psm == P_LINEAR)
      {
        return sample_bilinear(sp.p_uv, level_num_int);
      }
    }
    else if (sp.lsm == L_LINEAR)
    {
      if (sp.psm == P_NEAREST)
      {
        float low_level_weight = level_num - floor(level_num);
        float high_level_weight = ceil(level_num) - level_num;
        cout<<"low_level_wei:"<<low_level_weight<<endl;
        cout<<"low_level:"<<floor(level_num)<<endl;
        cout<<"high_level:"<<ceil(level_num)<<endl;
        cout<<"high_level_wei:"<<high_level_weight<<endl;
        cout<<"float level"<<level_num<<endl;
        Color low_level_color, high_level_color;
        low_level_color = sample_nearest(sp.p_uv, floor(level_num));
        high_level_color = sample_nearest(sp.p_uv, ceil(level_num));

        return low_level_color*low_level_weight + high_level_color*high_level_weight;
      }
      else if (sp.psm == P_LINEAR)
      {
        float low_level_weight = level_num - floor(level_num);
        float high_level_weight = ceil(level_num) - level_num;
        Color low_level_color, high_level_color;
        cout<<"LOW!!!!!"<<endl<<endl<<endl;
        low_level_color = sample_bilinear(sp.p_uv, floor(level_num));
        cout<<"High!!!!!"<<endl;
        high_level_color = sample_bilinear(sp.p_uv, ceil(level_num));

        return low_level_color*low_level_weight + high_level_color*high_level_weight;
      }
    }
// return magenta for invalid level
    return Color(1, 0, 1);
  }

  float Texture::get_level(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
    auto& mip = mipmap[0];

    Vector2D sp_uv = sp.p_uv;
    Vector2D sp_uv_dx = sp.p_dx_uv;
    Vector2D sp_uv_dy = sp.p_dy_uv;
    float dx_ans, dy_ans;
    dx_ans = pow(pow(sp_uv_dx.x*mip.width - sp_uv.x*mip.width, 2) + pow(sp_uv_dx.y*mip.height - sp_uv.y*mip.height, 2), 0.5);
    dy_ans = pow(pow(sp_uv_dy.x*mip.width - sp_uv.x*mip.width, 2) + pow(sp_uv_dy.y*mip.height - sp_uv.y*mip.height, 2), 0.5);
    float level_ans;
    level_ans = max(dx_ans, dy_ans);
    level_ans = log(level_ans)/log(2);

    return level_ans;
  }

  Color MipLevel::get_texel(int tx, int ty) {
    return Color(&texels[tx * 3 + ty * width * 3]);
  }

  Color Texture::sample_nearest(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    cout<<level<<endl;
    if (level == -1)
    {
      level = 0;
    }
    auto& mip = mipmap[level];
    float uv_x = uv.x;
    float uv_y = uv.y;
    int x_int = round(uv_x*mip.width);
    cout<<"uv.x:"<<uv_x<<" uv.y:"<<uv_y<<endl;
    int y_int = round(uv_y*mip.height);
    cout<<"uv.x.int:"<<x_int<<" uv.y.int:"<<y_int<<endl;
    Color Texture_color = mip.get_texel(x_int,y_int);
    // return magenta for invalid level
    return Texture_color;
  }

  Color Texture::sample_bilinear(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    if (level == -1)
    {
      level = 0;
    }
    auto& mip = mipmap[level];
    cout<<"level"<<level;
    
    float uv_x = uv.x*mip.width;
    float uv_y = uv.y*mip.height;
    int x_int, y_int;
    x_int = floor(uv_x);
    y_int = floor(uv_y);
    float xlength_0_0 = uv_x - x_int; float ylength_0_0 = uv_y - y_int;
    float xlength_0_1 = 1 - (uv_x - x_int); float ylength_0_1 = uv_y - y_int;
    float xlength_1_0 = uv_x - x_int; float ylength_1_0 = 1 - (uv_y - y_int);
    float xlength_1_1 = 1 - (uv_x - x_int); float ylength_1_1 = 1 - (uv_y - y_int);

    float color_cover_0_0 = xlength_0_0 * ylength_0_0;
    cout<<"uvx:"<<uv_x<<" uvx int:"<<x_int<<endl;
    cout<<"00x:"<<xlength_0_0<<" 00y:"<<ylength_0_0<<endl;
    cout<<"01x:"<<xlength_0_1<<" 01y:"<<ylength_0_1<<endl;
    cout<<"10x:"<<xlength_1_0<<" 10y:"<<ylength_1_0<<endl;
    cout<<"11x:"<<xlength_1_1<<" 11y:"<<ylength_1_1<<endl;
    float color_cover_0_1 = xlength_0_1 * ylength_0_1;
    float color_cover_1_0 = xlength_1_0 * ylength_1_0;
    float color_cover_1_1 = xlength_1_1 * ylength_1_1;
    cout<<"color plus:"<<color_cover_0_0+color_cover_0_1+color_cover_1_0+color_cover_1_1<<endl;
    cout<<"color cover 00:"<<color_cover_0_0<<endl;
    cout<<"color cover 01:"<<color_cover_0_1<<endl;
    cout<<"color cover 10:"<<color_cover_1_0<<endl;
    cout<<"color cover 11:"<<color_cover_1_1<<endl;
    Color Texture_color = mip.get_texel(x_int, y_int) * color_cover_0_0
                        + mip.get_texel(x_int, y_int+1) * color_cover_0_1
                        + mip.get_texel(x_int+1, y_int) * color_cover_1_0
                        + mip.get_texel(x_int+1, y_int+1) * color_cover_1_1;

    // return magenta for invalid level
    return Texture_color;
  }



  /****************************************************************************/

  // Helpers

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
