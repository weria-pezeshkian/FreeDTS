/*
  Parser for the CCP4/MRC format.

  Copyright (c) 202-2025 European Molecular Biology Laboratory

  Author: Valentin Maurer <valentin.maurer@embl-hamburg.de>
*/
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <cstdint>

template <typename T> float to_float(T c) { return static_cast<float>(c); }

class MRCParser {
private:
  struct Header {
    int32_t nx, ny, nz;                // 12
    int32_t mode;                      // 16
    int32_t nxstart, nystart, nzstart; // 28
    float mx, my, mz;                  // 40
    float cella[3], cellb[3];          // 64
    int32_t mapc, mapr, maps;          // 64
    float dmin, dmax, dmean;           // 76
    int32_t ispg;                      // 88
    int32_t nsymbt;                    // 96
    char extra[100];                   // 196
    int32_t originx, originy, originz; // 208
    char map[4];                       // 212
    int32_t machst;                    // 216
    float rms;                         // 220
    int32_t nlabl;                     // 224
    char label[800];                   // 1024
  };

  Header m_Header;
  std::vector<float> m_Data;

  template <typename T>
  void readData(std::ifstream &file, size_t dataSize) {
      m_Data.resize(dataSize);

      const size_t chunkSize = std::min(dataSize, size_t(1024 * 1024 * 64));
      std::vector<T> buffer(chunkSize);

      for (size_t offset = 0; offset < dataSize; offset += chunkSize) {
          size_t currentChunkSize = std::min(chunkSize, dataSize - offset);
          buffer.resize(currentChunkSize);
          file.read(reinterpret_cast<char *>(buffer.data()), currentChunkSize * sizeof(T));
          std::transform(buffer.begin(), buffer.end(), m_Data.begin() + offset, to_float<T>);
      }
  }

  void readFile(const std::string &filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
      std::cerr << "Cannot open file: " + filename + "\n";
      exit(1);
    }

    file.read(reinterpret_cast<char *>(&m_Header), sizeof(m_Header));
    if (std::string(m_Header.map, 4) != "MAP ") {
      std::cerr << "Invalid MRC file: MAP m_Header not found\n";
      exit(1);
    }

    // Factor in extended header
    if (m_Header.nsymbt > 0) {
      file.seekg(m_Header.nsymbt, std::ios::cur);
    }

    size_t dataSize = m_Header.nx * m_Header.ny * m_Header.nz;

    switch (m_Header.mode) {
    case 0: // 8-bit signed integer
      readData<int8_t>(file, dataSize);
      break;
    case 1: // 16-bit signed integer
      readData<int16_t>(file, dataSize);
      break;
    case 2: // 32-bit float
      readData<float>(file, dataSize);
      break;
    case 6: // 16-bit unsigned integer
      readData<uint16_t>(file, dataSize);
      break;
    default:
      std::cerr << "Unsupported MODE: " + std::to_string(m_Header.mode);
      exit(1);
    }
  }

public:
  MRCParser(const std::string &filename) { readFile(filename); }

  const std::vector<float> &getData() const { return m_Data; }
  const std::vector<int32_t> getShape() const {
    const std::vector<int32_t> shape = {m_Header.nx, m_Header.ny, m_Header.nz};
    return shape;
  }
  const std::vector<float> getSampling() const {
    const std::vector<float> sampling = {m_Header.cella[0] / m_Header.nx,
                                      m_Header.cella[1] / m_Header.ny,
                                      m_Header.cella[2] / m_Header.nz};
    return sampling;
  }
};
