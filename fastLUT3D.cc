#include "fastLUT3D.h"

fastLUT3D::fastLUT3D(int array_size)
{
  key_array_size = pow(array_size, 3);
  one_dim_size = array_size;

  GPU_mode = false; // off by default
  
  init();
}

fastLUT3D::~fastLUT3D()
{
  delete [] xyz_array;

#ifdef __CUDACC__
  cudaFree(gpu_params);
  cudaFreeArray(array);
  cudaDestroyTextureObject(tex);
#endif
}

void fastLUT3D::init()
{
  xyz_array = new float[key_array_size];
}

void fastLUT3D::setElement(uint32_t x, uint32_t y, uint32_t z, float value, bool morton)
{
  uint32_t key;
  if (morton)
    {
      key = EncodeMorton(x, y, z);
    }
  else
    {
      key = index(x, y, z);
    }
  xyz_array[key] = value;
}

inline float fastLUT3D::getElement(uint32_t x, uint32_t y, uint32_t z, bool morton)
{
  uint64_t key;
  if (morton)
    {
      key = EncodeMorton(x, y, z);
    }
  else
    {
      key = index(x, y, z);
    }
  return xyz_array[key];
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

#ifdef __CUDACC__
  float cpu_params[4] = {xbound[0], ybound[0], zbound[0], inv_cell_size[0]};
  cudaMalloc((void**)&gpu_params, 4 * sizeof(float));
  cudaMemcpy(gpu_params, cpu_params, 4 * sizeof(float), cudaMemcpyHostToDevice);
#endif
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
  //return answer;
  
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

HEMI_DEV_CALLABLE_MEMBER float fastLUT3D::Interpolate(float x, float y, float z, bool morton)
{
  #ifdef HEMI_DEV_CODE
  float invcellsize = gpu_params[4];
  float xd = (x - gpu_params[0]) * invcellsize + 0.5f;
  float yd = (y - gpu_params[1]) * invcellsize + 0.5f;
  float zd = (z - gpu_params[2]) * invcellsize + 0.5f;
  
  return tex3D<float>(tex, xd, yd, zd);
  #else

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
  #endif
}

/*void fastLUT3D::makeMortonPlot(TH3I &h)
{
  for (int i = 0; i < h.GetNbinsX(); ++i)
    for (int j = 0; j < h.GetNbinsY(); ++j)
      for (int k = 0; k < h.GetNbinsZ(); ++k)
	h.SetBinContent(i+1, j+1, k+1, getElement(i,j,k));//EncodeMorton(i,j,k));
	}*/

void fastLUT3D::GPUmode(bool on)
{
  if (on)
    {
      #ifdef __CUDACC__
      GPU_mode = true;
      createCudaArray();
      createTexture();
      #else
      printf("ERROR: GPU mode selected, but not compiled with GPU support.\n");
      exit(-1);
      #endif
    }
  else
    {
      GPU_mode = false;
    }
}
  
// gpu only functions
#ifdef __CUDACC__

void fastLUT3D::createCudaArray()
{
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  cudaMemcpy3DParms copyParams = {0};
  
  cudaExtent volumeSize = make_cudaExtent(one_dim_size, one_dim_size, one_dim_size);
  cudaMalloc3DArray(&array, &channelDesc, volumeSize);
  copyParams.srcPtr = make_cudaPitchedPtr((void*)xyz_array, volumeSize.width * sizeof(float), volumeSize.width, volumeSize.height);
  copyParams.dstArray = array;
  copyParams.extent   = volumeSize;
  copyParams.kind     = cudaMemcpyHostToDevice;
  cudaMemcpy3D(&copyParams);
}

void fastLUT3D::createTexture()
{
  // created array, now to bind it to a texture object
  struct cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType =  cudaResourceTypeArray;

  resDesc.res.linear.devPtr = array;
  resDesc.res.linear.sizeInBytes = one_dim_size * one_dim_size * one_dim_size * sizeof(float);

  resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
  resDesc.res.linear.desc.x = 32;
  struct cudaTextureDesc texDesc;
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.normalizedCoords = false;
  //texDesc.normalizedCoords = true;
  texDesc.readMode = cudaReadModeElementType;
  //texDesc.filterMode = cudaFilterModePoint;
  texDesc.filterMode = cudaFilterModeLinear;

  texDesc.addressMode[0] = cudaAddressModeClamp;
  texDesc.addressMode[1] = cudaAddressModeClamp;
  texDesc.addressMode[2] = cudaAddressModeClamp;

  // make the texture object
  cudaCreateTextureObject(&tex, &resDesc, &texDesc, NULL);
}

#endif
