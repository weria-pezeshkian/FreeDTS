/*
  Compute the potential of a point in a scalar field

  U_{em}(R) = \sum_j w_j \cdot V_{em}(r_j)
  V_{em}(r) = \begin{cases}
    \m_Xi \left[1 - \frac{\theta(r) - \theta_{thr}}{\theta_{max} -
  \theta_{thr}}\right] & \text{if } \theta(r) \geq \theta_{thr} \ \m_Xi &
  \text{if } \theta(r) < \theta_{thr} \end{cases}

  Where:

  - $\theta(r)$ is the scalar field containing EM densities.
  - $r$ is a position vector in 3D space.
  - $\m_Xi$ is a force scaling factor.
  - $\theta_{max}$ is the mam_Ximum value of the densities.
  - $\theta_{thr}$ is a threshold value to exclude solvent contribution.
  - $w_j$ is a per-atom weighting factor.

  See the MDFF formalism (https://www.ks.uiuc.edu/Research/mdff/method.html)

  L.G. Trabuco, E. Villa, K. Mitra, J. Frank, and K. Schulten. Structure, 16,
  673-683 (2008)

  Copyright (c) 2024-2025 European Molecular Biology Laboratory

  Author: Valentin Maurer <valentin.maurer@embl-hamburg.de>
*/
#include <cmath>
#include <numeric>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>

inline std::vector<float> gaussDerivative(
    const unsigned int dim, const unsigned int kernel_size, const double sigma = 1.0
  ){
    if (dim >= 3) {
        throw std::invalid_argument("dim must be < 3");
    }

    const int m_Offset = kernel_size / 2;
    std::vector<float> diff(kernel_size);
    std::vector<float> smooth(kernel_size);
    std::vector<float> ret(kernel_size * kernel_size * kernel_size);

    const float sigma2 = sigma * sigma;
    for (int i = 0; i < kernel_size; ++i) {
        float x = -(static_cast<float>(i) - m_Offset);
        smooth[i] = std::exp(-(x * x)/(2.0 * sigma2));
        diff[i] = smooth[i] * (-x/(sigma2));
    }

    for (int k = 0; k < kernel_size; ++k) {
        float z_val = (dim == 2) ? diff[k] : smooth[k];
        for (int j = 0; j < kernel_size; ++j) {
            float y_val = (dim == 1) ? diff[j] : smooth[j];
            for (int i = 0; i < kernel_size; ++i) {
                float x_val = (dim == 0) ? diff[i] : smooth[i];
                ret[i + kernel_size * (j + kernel_size * k)] = x_val * y_val * z_val;
            }
        }
    }

    return ret;
}

class HarmonicPotentialCalculator {
private:
  const std::vector<float> m_Data; // Volume (x,y,z) with length N = x * y * z.
  const std::vector<unsigned int> m_Shape;  // Volume shape x, y, z.
  const std::vector<double> m_Sampling; // Distance/Voxel conversion factor.
  const std::vector<double> m_Offsets;       // Grid offsets from centering.
  const double m_DataScaling; // Scaling factor for data to invert contrast.
  const double m_Xi;       // Force scaling factor.
  const double m_ThetaThr; // Threshold for density to be foreground.
  double m_ThetaMax;       // Maximum value in m_Data.
  const double m_StepSize; // Step size along gradient.
  const unsigned int m_KernelSize; // Gradient kernel size.
  const size_t m_StrideY; // y stride of m_Data.
  const size_t m_StrideZ; // z stride of m_Data.
  std::vector<float> m_gradient_X; // Gradient kernel in x.
  std::vector<float> m_gradient_Y; // Gradient kernel in y.
  std::vector<float> m_gradient_Z; // Gradient kernel in z.
  std::vector<float> m_boundaryValues; // Per-plane averages used for padding.
  const float m_Padding = 0.4; // Scaling factor for padding values.

  inline bool isInBounds(const unsigned int x, const unsigned int y, const unsigned int z) const noexcept {
    return x < m_Shape[0] && y < m_Shape[1] && z < m_Shape[2];
  }

  inline unsigned int safeDoubleToUInt(const double val) const noexcept {
    return val < 0.0 ? std::numeric_limits<unsigned int>::max() : static_cast<unsigned int>(val);
  }

  __attribute__((always_inline))
  inline float getValue(const int x, const int y, const int z) const noexcept {
    return m_Data[z * m_StrideZ + y * m_StrideY + x];
  }

  inline float interpolateTheta(const double x, const double y, const double z) const noexcept {
    // Wrapping will cause negative values to fail in bounds check
    const unsigned int x0 = safeDoubleToUInt(x);
    const unsigned int y0 = safeDoubleToUInt(y);
    const unsigned int z0 = safeDoubleToUInt(z);

    if (!isInBounds(x0, y0, z0) || !isInBounds(x0 + 1, y0 + 1, z0 + 1)){

      int offset = z0;
      if (z0 >= m_Shape[2]) {
          offset = m_Shape[2] - 1;
          if (z < 0){
            offset = 0;
          }
      }

      return m_Padding * m_boundaryValues[offset] * m_ThetaMax;
    }

    const float xd = static_cast<float>(x - x0);
    const float yd = static_cast<float>(y - y0);
    const float zd = static_cast<float>(z - z0);

    const float xd_comp = 1.0f - xd;
    const float yd_comp = 1.0f - yd;
    const float zd_comp = 1.0f - zd;

    const float v000 = getValue(x0, y0, z0);
    const float v100 = getValue(x0 + 1, y0, z0);
    const float v010 = getValue(x0, y0 + 1, z0);
    const float v110 = getValue(x0 + 1, y0 + 1, z0);
    const float v001 = getValue(x0, y0, z0 + 1);
    const float v101 = getValue(x0 + 1, y0, z0 + 1);
    const float v011 = getValue(x0, y0 + 1, z0 + 1);
    const float v111 = getValue(x0 + 1, y0 + 1, z0 + 1);

    const float c00 = v000 * xd_comp + v100 * xd;
    const float c10 = v010 * xd_comp + v110 * xd;
    const float c01 = v001 * xd_comp + v101 * xd;
    const float c11 = v011 * xd_comp + v111 * xd;

    const float c0 = c00 * yd_comp + c10 * yd;
    const float c1 = c01 * yd_comp + c11 * yd;

    return m_DataScaling * (c0 * zd_comp + c1 * zd);
  }

  std::array<double, 3> interpolateThetaGrad(const double x, const double y, const double z) const {
      std::array<double, 3> gradient = {0.0, 0.0, 0.0};

      if (m_StepSize == 0.0f){
        return gradient;
      }

      const unsigned int x0 = safeDoubleToUInt(x);
      const unsigned int y0 = safeDoubleToUInt(y);
      const unsigned int z0 = safeDoubleToUInt(z);

      if (!isInBounds(x0, y0, z0) || !isInBounds(x0 + 1, y0 + 1, z0 + 1)){
        return gradient;
      }

      // For larger number of vertices, consider pre-computing the volume gradient
      const int n = m_KernelSize / 2;
      for (int k = 0; k < m_KernelSize; ++k) {
          const int zkn = z + k - n;
          const int km = m_KernelSize * k;

          for (int j = 0; j < m_KernelSize; ++j) {
              const int ykn = y + j - n;
              const int jkm = m_KernelSize * (j + km);

              for (int i = 0; i < m_KernelSize; ++i) {
                  const int idx = i + jkm;
                  const double val = static_cast<double>(interpolateTheta(x + i - n, ykn, zkn));

                  gradient[0] += val * m_gradient_X[idx];
                  gradient[1] += val * m_gradient_Y[idx];
                  gradient[2] += val * m_gradient_Z[idx];
              }
          }
      }

      const double norm = std::sqrt(gradient[0] * gradient[0] +
                            gradient[1] * gradient[1] +
                            gradient[2] * gradient[2]);
      if (norm > 1e-10) {
        const double scaleFactor = m_StepSize / norm;
        gradient[0] *= scaleFactor;
        gradient[1] *= scaleFactor;
        gradient[2] *= scaleFactor;
      }
      return gradient;
  }


void printDebugInfo() const {
    const int width = 16; // Increased for longer parameter names

    std::cout << "\n===== HMFF Configuration =====\n";

    // Volume information
    std::cout << std::left << std::setw(width) << "Reference" <<
      ": (X, Y, Z)" << "\n";
    std::cout << std::left << std::setw(width) << "Volume Shape" <<
      ": (" << m_Shape[0] << ", " << m_Shape[1] << ", " << m_Shape[2] << ")\n";
    std::cout << std::left << std::setw(width) << "Mesh Scaling" <<
      ": (" << std::fixed << std::setprecision(3)
            << m_Sampling[0] << ", "
            << m_Sampling[1] << ", "
            << m_Sampling[2] << ")\n";
    std::cout << std::left << std::setw(width) << "Mesh Offset" <<
      ": (" << std::fixed << std::setprecision(3)
            << m_Offsets[0] << ", "
            << m_Offsets[1] << ", "
            << m_Offsets[2] << ")\n";

    std::cout << std::left << std::setw(width) << "Xi"
      << ": " << m_Xi << "\n";
    std::cout << std::left << std::setw(width) << "ThetaThr"
      << ": " << m_ThetaThr << "\n";
    std::cout << std::left << std::setw(width) << "ThetaMax"
      << ": " << m_ThetaMax << "\n";
    std::cout << std::left << std::setw(width) << "StepSize"
      << ": " << m_StepSize << "\n";
    std::cout << std::left << std::setw(width) << "KernelSize"
      << ": " << m_KernelSize << "\n";
    std::cout << std::left << std::setw(width) << "DataScaling"
      << ": " << m_DataScaling << "\n";
    std::cout << std::left << std::setw(width) << "PaddingFactor"
      << ": " << m_Padding << "\n\n";
}


public:
  HarmonicPotentialCalculator(const std::vector<float> &data,
                              const std::vector<unsigned int> &shape,
                              const std::vector<double> &sampling,
                              const double m_Xi,
                              const double m_ThetaThr,
                              const std::vector<double> &offsets,
                              const double m_DataScaling,
                              const double step_size,
                              const unsigned int kernel_size = 3,
                              float data_percentile = 0.999,
                              const double padding_factor = 0.4)
      : m_Data(data),
        m_Shape(shape),
        m_Sampling(sampling), m_Xi(m_Xi),
        m_ThetaThr(m_ThetaThr),
        m_Offsets(offsets),
        m_DataScaling(m_DataScaling),
        m_StepSize(step_size),
        m_KernelSize(kernel_size),
        m_StrideY(shape[0]),
        m_StrideZ(shape[0] * shape[1]),
        m_Padding(padding_factor) {

    if (shape.size() != sampling.size()) {
      std::cerr << "Shape and sampling rate dimensions must match.\n";
      exit(1);
    }

    if (offsets.size() != 3) {
      std::cerr << "Offsets vector must contain exactly 3 values (x,y,z).\n";
      exit(1);
    }

    size_t totalSize = shape[0] * shape[1] * shape[2];
    if (totalSize != data.size()) {
      std::cerr << "Data size doesn't match the specified shape.\n";
      exit(1);
    }

    std::vector<float> dataCopy = data;
    float rank = (m_DataScaling > 0 ? data_percentile : 1.0 - data_percentile);
    size_t k = std::floor(rank * (totalSize - 1));
    std::nth_element(dataCopy.begin(), dataCopy.begin() + k, dataCopy.end());
    float theta_max = m_DataScaling * dataCopy[k];

    m_ThetaMax = static_cast<double>(theta_max);

    m_gradient_X = gaussDerivative(0, m_KernelSize, 1.0);
    m_gradient_Y = gaussDerivative(1, m_KernelSize, 1.0);
    m_gradient_Z = gaussDerivative(2, m_KernelSize, 1.0);

    // Neighbor weighted per slice z-value
    m_boundaryValues.resize(shape[2]);
    for (size_t z = 0; z < shape[2]; ++z) {
        m_boundaryValues[z] = computePlaneAverage(data, z);
    }

    std::vector<float> temp = m_boundaryValues;
    for (size_t z = 0; z < shape[2]; ++z) {
        if (z == 0) {
            m_boundaryValues[z] = (temp[z] + temp[z + 1]) / 2.0f;
        }
        else if (z == shape[2] - 1) {
            m_boundaryValues[z] = (temp[z - 1] + temp[z]) / 2.0f;
        }
        else {
            m_boundaryValues[z] = (temp[z - 1] + temp[z] + temp[z + 1]) / 3.0f;
        }
    }

    float minVal = *std::min_element(m_boundaryValues.begin(), m_boundaryValues.end());
    float maxVal = *std::max_element(m_boundaryValues.begin(), m_boundaryValues.end());
    float range = maxVal - minVal;

    if (range > 0) {
        for (size_t i = 0; i < m_boundaryValues.size(); ++i) {
            m_boundaryValues[i] = (m_boundaryValues[i] - minVal) / range;
        }
    }
    else {
        std::fill(m_boundaryValues.begin(), m_boundaryValues.end(), 0.5f);
    }

    printDebugInfo();
  }


  float computePlaneAverage(const std::vector<float>& inputData, size_t zIndex) {
      size_t planeSize = m_Shape[0] * m_Shape[1];
      size_t startIdx = zIndex * planeSize;
      size_t endIdx = startIdx + planeSize;

      float sum = 0.0f;
      for (size_t i = startIdx; i < endIdx; ++i) {
          sum += std::abs(inputData[i]);
      }
      return sum / planeSize;
  }

  inline double computeEnergy(const double x, const double y, const double z) const noexcept {

    const double scaled_x = (x - m_Offsets[0]) * m_Sampling[0];
    const double scaled_y = (y - m_Offsets[1]) * m_Sampling[1];
    const double scaled_z = (z - m_Offsets[2]) * m_Sampling[2];

    const double theta = std::min(
      static_cast<double>(interpolateTheta(scaled_x, scaled_y, scaled_z)), m_ThetaMax
    );

    if (theta > m_ThetaThr) {
      return m_Xi * (1 - (theta - m_ThetaThr) / (m_ThetaMax - m_ThetaThr));
    }
    return m_Xi;

  }

  inline std::array<double, 3> computeEnergyGradient(const double x, const double y,
                              const double z) const {
    const double scaled_x = (x - m_Offsets[0]) * m_Sampling[0];
    const double scaled_y = (y - m_Offsets[1]) * m_Sampling[1];
    const double scaled_z = (z - m_Offsets[2]) * m_Sampling[2];

    return interpolateThetaGrad(scaled_x, scaled_y, scaled_z);
  }

};
