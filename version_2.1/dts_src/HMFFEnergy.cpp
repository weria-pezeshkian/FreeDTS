/*
  Compute the energy of a whole system in a scalar field using HMFF.

  Copyright (c) 2024-2025 European Molecular Biology Laboratory

  Author: Valentin Maurer <valentin.maurer@embl-hamburg.de>
*/
#include "HMFFEnergy.h"
#include "Nfunction.h"
#include "State.h"
#include <time.h>
#include <sstream>

HMFFEnergy::HMFFEnergy(State *pState, std::string data) : Energy(pState){
  std::vector<std::string> ndata = Nfunction::Split(data);

  if(ndata.size() < 5 || ndata.size() > 8){
    std::cerr << "MDFF Energy Parameters:\n"
              << "  Required: path xi theta_thr scaling offset\n"
              << "  Optional: [inversion] [step_size] [padding_factor] [kernel_size] [data_percentile]\n\n"
              << "  path       : Path to MRC file containing EM density map\n"
              << "  xi         : Force scaling factor (typical: 1.0-50.0)\n"
              << "  theta_thr  : Threshold value to exclude solvent contribution (typical 0.0) \n"
              << "  scaling    : Mesh scaling factor for distance/voxel conversion\n"
              << "  offset     : Grid offset (single value or comma-separated x,y,z)\n"
              << "  inversion  : Invert density contrast (0=false, 1=true, default: 1)\n"
              << "  step_size  : Step size along gradient (default: 0.0)\n"
              << "  padding    : Scaling factor for boundary padding (default: 0.4)\n"
              << "  kernel     : Gradient kernel size, must be odd (default: 3)\n"
              << "  percentile : theta_max determined as percentile of data (0.0-1.0, default: 0.999)\n\n"
              << "Examples:\n"
              << "  Basic: path.mrc 5.0 0 0.005 24.0\n"
              << "  Full:  path.mrc 5.0 0 0.005 24.0 1 0.0 0.4\n";
    exit(1);
  }

  // Required parameters
  const double xi = (float) Nfunction::String_to_Double(ndata[1]);
  const double theta_thr = (float) Nfunction::String_to_Double(ndata[2]);
  const double mesh_scaling = Nfunction::String_to_Double(ndata[3]);

  // Parse mesh offset from (0,0,0).
  std::vector<double> offsets;
  std::istringstream ss(ndata[4]);
  std::string value;
  std::vector<std::string> offset_parts;
  while (std::getline(ss, value, ',')) {
    offset_parts.push_back(value);
  }

  if (offset_parts.size() == 1) {
    double single_offset = Nfunction::String_to_Double(offset_parts[0]);
    offsets = {single_offset, single_offset, single_offset};
  } else if (offset_parts.size() == 3) {
    offsets = {
      Nfunction::String_to_Double(offset_parts[0]),
      Nfunction::String_to_Double(offset_parts[1]),
      Nfunction::String_to_Double(offset_parts[2])
    };
  } else {
    std::cerr << "Error: offset must be either a single value or three comma-separated values (x,y,z)\n";
    exit(1);
  }

  double data_scaling = 1.0;
  if (ndata.size() >= 6) {
    if ((bool) Nfunction::String_to_Double(ndata[5])) {
      data_scaling = -1.0;
    }
  }

  const double step_size = (ndata.size() >= 7) ?
    (double) Nfunction::String_to_Double(ndata[6]) : 0.0;

  const float padding_factor = (ndata.size() >= 8) ?
    (float) Nfunction::String_to_Double(ndata[7]) : 0.4;

  unsigned int kernel_size = 3;
  if (ndata.size() >= 9) {
    kernel_size = (unsigned int) Nfunction::String_to_Double(ndata[8]);
    if (kernel_size % 2 == 0) {
      std::cerr << "Error: kernel_size must be odd (got " << kernel_size << ")\n";
      exit(1);
    }
    if (kernel_size < 3) {
      std::cerr << "Error: kernel_size must be >= 3 (got " << kernel_size << ")\n";
      exit(1);
    }
  }

  float data_percentile = 0.999f; // Default
  if (ndata.size() >= 10) {
    data_percentile = (float) Nfunction::String_to_Double(ndata[9]);
    if (data_percentile <= 0.0f || data_percentile >= 1.0f) {
      std::cerr << "Error: data_percentile must be between 0.0 and 1.0 (got "
                << data_percentile << ")\n";
      exit(1);
    }
  }

  MRCParser parser(ndata[0]);
  const std::vector<float> &mrcData = parser.getData();
  std::vector<int32_t> int32t_shape = parser.getShape();
  const std::vector<unsigned int> shape(int32t_shape.begin(), int32t_shape.end());
  std::vector<float> float_sampling = parser.getSampling();
  std::vector<double> sampling(float_sampling.begin(), float_sampling.end());

  // Transform mesh to 1U/Voxel and then bin to SamplingU/Voxel
  for (size_t i = 0; i < sampling.size(); i++) {
      sampling[i] *= mesh_scaling;
      sampling[i] = 1 / sampling[i];
  }

  m_pCalculator = new HarmonicPotentialCalculator(
    mrcData, shape, sampling, xi, theta_thr, offsets, data_scaling, step_size,
    kernel_size, data_percentile, padding_factor
  );
  if (m_pCalculator == nullptr) {
    std::cerr << "Error: Failed to allocate memory for "
                 "Inclusion_Interaction_Map \n";
    exit(1);
  }
}

Vec3D HMFFEnergy::CalculateGradient(vertex *p_vertex) {
    const std::array<double, 3> grad = m_pCalculator->computeEnergyGradient(
        p_vertex->GetXPos(), p_vertex->GetYPos(), p_vertex->GetZPos()
    );
    return Vec3D(grad[0], grad[1], grad[2]);
}


double HMFFEnergy::SingleVertexEnergy(vertex *p_vertex) {
    // Get the standard energy from the base Energy class
    double energy = Energy::SingleVertexEnergy(p_vertex);

    // Add HMFF contribution
    energy += m_pCalculator->computeEnergy(
        p_vertex->GetXPos(), p_vertex->GetYPos(), p_vertex->GetZPos()
    );

    // Update the vertex with the new total energy
    // Note: we are updating the vertex energy twice. This could be avoided
    // by making the corresponding members in Energy non-private or re-implementing
    // them here, risking synchronization issues in the future.
    p_vertex->UpdateEnergy(energy);
    return energy;
}
