#include "fastLUT3D.h"

fastLUT3D::fastLUT3D(int array_size)
{
  key_array_size = pow(array_size, 3);
  one_dim_size = array_size;
  
  init();
}

fastLUT3D::~fastLUT3D()
{
  delete [] key_array;
  delete [] xyz_array;
}

void fastLUT3D::init()
{
  key_array = new float[key_array_size];
  xyz_array = new float[key_array_size];
}

void fastLUT3D::setElement(uint32_t x, uint32_t y, uint32_t z, float value, bool morton)
{
  if (morton)
    {
      uint32_t key = EncodeMorton(x, y, z);
      key_array[key] = value;
    }
  else
    {
      uint32_t key = index(x, y, z);
      xyz_array[key] = value;
    }
}

inline float fastLUT3D::getElement(uint32_t x, uint32_t y, uint32_t z, bool morton)
{
  if (morton)
    {
      uint64_t key = EncodeMorton(x, y, z);
      return key_array[key];
    }
  else
    {
      uint32_t key = index(x, y, z);
      return xyz_array[key];
    }
}

void fastLUT3D::setBounds(float xlo, float xhi, float ylo, float yhi, float zlo, float zhi)
{
  xbound[0] = xlo; xbound[1] = xhi;
  ybound[0] = ylo; ybound[1] = yhi;
  zbound[0] = zlo; zbound[1] = zhi;

  range[0] = xhi - xlo;
  range[1] = yhi - ylo;
  range[2] = zhi - zlo;

  cell_size[0] = range[0] / (one_dim_size);
  cell_size[1] = range[1] / (one_dim_size);
  cell_size[2] = range[2] / (one_dim_size);

  inv_cell_size[0] = 1.0f / cell_size[0];
  inv_cell_size[1] = 1.0f / cell_size[1];
  inv_cell_size[2] = 1.0f / cell_size[2];
}

inline uint32_t fastLUT3D::index(uint32_t x, uint32_t y, uint32_t z)
{
  return x + one_dim_size * (y + one_dim_size * z);
}

inline uint64_t fastLUT3D::EncodeMorton(uint32_t x, uint32_t y, uint32_t z)
{
  uint64_t answer = 0;
  answer = morton256_z[(z >> 16) & 0xFF ] | // we start by shifting the third byte, since we only look at the first 21 bits
    morton256_y[(y >> 16) & 0xFF ] |
    morton256_x[(x >> 16) & 0xFF ];
  answer = answer << 48 | morton256_z[(z >> 8) & 0xFF ] | // shifting second byte
    morton256_y[(y >> 8) & 0xFF ] |
    morton256_x[(x >> 8) & 0xFF ];
  answer = answer << 24 |
    morton256_z[(z) & 0xFF ] | // first byte
    morton256_y[(y) & 0xFF ] |
    morton256_x[(x) & 0xFF ];
  return answer;
  
  return answer;//(Part1By2(z) << 2) + (Part1By2(y) << 1) + Part1By2(x);
}

inline uint32_t fastLUT3D::Part1By2(uint32_t x)
{
  x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
  x = (x ^ (x << 16)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
  x = (x ^ (x <<  8)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
  x = (x ^ (x <<  4)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
  x = (x ^ (x <<  2)) & 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
  return x;
}

void fastLUT3D::DecodeMorton(uint32_t code, uint32_t &x, uint32_t &y, uint32_t &z)
{
  x = Compact1By2(code >> 0);
  y = Compact1By2(code >> 1);
  z = Compact1By2(code >> 2);
}

inline uint32_t fastLUT3D::Compact1By2(uint32_t x)
{
  x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
  x = (x ^ (x >>  2)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
  x = (x ^ (x >>  4)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
  x = (x ^ (x >>  8)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
  x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
  return x;
}

float fastLUT3D::Interpolate(float x, float y, float z, bool morton)
{
  float xd = (x - xbound[0]) * inv_cell_size[0];
  float yd = (y - ybound[0]) * inv_cell_size[1];
  float zd = (z - zbound[0]) * inv_cell_size[2];
    
  uint32_t ubx = (uint32_t)xd;
  uint32_t uby = (uint32_t)yd;
  uint32_t ubz = (uint32_t)zd;

  uint32_t obx = ubx + 1;
  uint32_t oby = uby + 1;
  uint32_t obz = ubz + 1;

  float v[] = { getElement( ubx, uby, ubz, morton ), getElement( ubx, uby, obz, morton ),
		getElement( ubx, oby, ubz, morton ), getElement( ubx, oby, obz, morton ),
		getElement( obx, uby, ubz, morton ), getElement( obx, uby, obz, morton ),
		getElement( obx, oby, ubz, morton ), getElement( obx, oby, obz, morton ) };

  xd -= (float)ubx;
  yd -= (float)uby;
  zd -= (float)ubz;
  
  float i1 = v[0] * (1 - zd) + v[1] * zd;
  float i2 = v[2] * (1 - zd) + v[3] * zd;
  float j1 = v[4] * (1 - zd) + v[5] * zd;
  float j2 = v[6] * (1 - zd) + v[7] * zd;

  
  float w1 = i1 * (1 - yd) + i2 * yd;
  float w2 = j1 * (1 - yd) + j2 * yd;

  //_mm_fmadd_ps()
  float result = w1 * (1 - xd) + w2 * xd;

  //std::cout << w1 << " " << w2 << " " << v[0] << std::endl;
  return result;
}

void fastLUT3D::makeMortonPlot(TH3I &h)
{
  for (int i = 0; i < h.GetNbinsX(); ++i)
    for (int j = 0; j < h.GetNbinsY(); ++j)
      for (int k = 0; k < h.GetNbinsZ(); ++k)
	h.SetBinContent(i+1, j+1, k+1, getElement(i,j,k));//EncodeMorton(i,j,k));
}
