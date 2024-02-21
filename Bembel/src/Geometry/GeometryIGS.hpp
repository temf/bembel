// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_SRC_GEOMETRY_GEOMETRYIGS_HPP_
#define BEMBEL_SRC_GEOMETRY_GEOMETRYIGS_HPP_

namespace Bembel {

/**
 * \ingroup Geometry
 * \brief loads  geometry from IGES file. Note that the direction
 *of the normals must be consistent.
 *
 * \param file_name path/filename pointing to the geometry file
 * \return std::vector of NURBS::Patch describing geometry
 */
std::vector<Patch> LoadGeometryFileIGS(const std::string& file_name) noexcept {
  std::ifstream file;
  std::vector<int> patch_lines;
  file.open(file_name);
  if (!file) {
    std::cerr << "File " << file_name << " doesn't exist!";
    exit(1);
  }
  std::string row;
  getline(file, row);

  if (row[72] != 'S') {
    std::cerr << "Format of the file not supported! I assume that character 72 "
                 "denotes the section";
    exit(1);
  }
  // character 72 denotes the section and the information starts in Section D
  while (row[72] != 'D') {
    getline(file, row);
    assert(!file.eof() && "End of file should not be found here!");
  }

  // collect two lines info of each patch from Directory section
  while (row[72] != 'P') {
    std::stringstream iss;
    std::string word;
    std::vector<int> info1, info2;
    iss.str(row);
    while (iss >> word) {
      info1.push_back(std::stoi(word));
    }
    assert(info1[0] == 128 && "Entry type must be NURBS surface!");
    // this info1 in not going to be used further

    getline(file, row);
    iss.clear();
    iss.str(row);
    while (iss >> word) {
      info2.push_back(std::stoi(word));
    }
    getline(file, row);
    // This entry denotes how many lines correspond to a patch
    patch_lines.push_back(info2[3]);
    assert(!file.eof() && "End of file should not be found here!");
  }

  // main loop over patches
  std::vector<Patch> out;
  out.reserve(patch_lines.size());
  const int number_of_patches = patch_lines.size();
  for (auto i = 0; i < number_of_patches; ++i) {
    std::stringstream iss;
    std::string word;
    std::vector<double> data;
    for (auto j = 0; j < patch_lines[i]; ++j) {
      // characters 0 to 63 contain data
      std::string raw_data = row.substr(0, 64);
      raw_data.erase(std::remove(raw_data.begin(), raw_data.end(), ' '),
                     raw_data.end());

      iss.str(raw_data);
      while (std::getline(iss, word, ',')) {
        data.push_back(atof(word.c_str()));
      }
      iss.clear();

      getline(file, row);
    }

    const int K1 = data[1];
    const int K2 = data[2];
    const int M1 = data[3];
    const int M2 = data[4];

    const int N1 = 1 + K1 - M1;
    const int N2 = 1 + K2 - M2;
    const int A = N1 + 2 * M1;
    const int B = N2 + 2 * M2;
    const int C = (1 + K1) * (1 + K2);

    // The data of the first entries is read and not needed any more
    data.erase(data.begin(), data.begin() + 10);

    // the + 1 is necessary because the construct excludes the last iterator
    std::vector<double> tempknt1(data.begin(), data.begin() + A + 1);
    std::vector<double> tempknt2(data.begin() + A + 1,
                                 data.begin() + A + B + 2);

    const int M = K2 + 1;
    const int N = K1 + 1;

    Eigen::MatrixXd weights(M, N);
    auto it_weights = data.begin() + A + B + 2;
    for (auto entry = it_weights; entry != it_weights + C; ++entry) {
      int index = entry - it_weights;
      weights((int)index / N, index % N) = *entry;
    }

    Eigen::MatrixXd x_coordinates(M, N);
    Eigen::MatrixXd y_coordinates(M, N);
    Eigen::MatrixXd z_coordinates(M, N);
    auto it_points = data.begin() + A + B + C + 2;
    for (auto entry = it_points; entry != it_points + 3 * C; entry += 3) {
      int index = (int)(entry - it_points) / 3;
      x_coordinates((int)index / N, index % N) = *entry;
      y_coordinates((int)index / N, index % N) = *(entry + 1);
      z_coordinates((int)index / N, index % N) = *(entry + 2);
    }

    // we need to transfer to homogeneous coordinates
    x_coordinates = x_coordinates.cwiseProduct(weights);
    y_coordinates = y_coordinates.cwiseProduct(weights);
    z_coordinates = z_coordinates.cwiseProduct(weights);

    std::vector<Eigen::MatrixXd> tmp = {x_coordinates, y_coordinates,
                                        z_coordinates, weights};

    Bembel::Patch tempPatch;
    tempPatch.init_Patch(tmp, tempknt1, tempknt2);
    out.push_back(tempPatch);
  }
  file.close();
  return out;
}

/**
 * \ingroup Geometry
 * \brief This function writes a given section to the specified file.
 *
 * \param file_name path/filename pointing to the geometry file
 * \param section vector containing the lines to be written to the file
 * \param section_char char to describe the section. Use 'S', 'G' or 'D'
 */
void writeSection(std::string file_name,
                  const std::vector<std::string>& section,
                  const char section_char) {
  std::ofstream file(file_name, std::ios::app);
  for (auto it = section.begin(); it != section.end(); ++it) {
    const int index = std::distance(section.begin(), it);
    file << std::left << std::setw(72) << *it << section_char << std::right
         << std::setw(7) << index + 1 << "\n";
  }
  file.close();
  return;
}

/**
 * \ingroup Geometry
 * \brief This function writes a the parameter section of the IGES file.
 *
 * \param file_name path/filename pointing to the geometry file
 * \param section vector containing the lines to be written to the file
 * \param patch_idx index counting in uneven number
 * \param start_count Start counting the lines with this index
 * \return Last index of the line 
 */
int writeParameterSection(std::string file_name,
                          const std::vector<std::string>& section,
                          const int patch_idx, const int start_count) {
  std::ofstream file(file_name, std::ios::app);
  int last_line = start_count;
  for (auto it = section.begin(); it != section.end(); ++it) {
    file << std::left << std::setw(64) << *it
         << std::right << std::setw(8) << patch_idx << "P"
         << std::right << std::setw(7) << last_line << "\n";
    ++last_line;
  }
  file.close();
  return last_line;
}

/**
 * \ingroup Geometry
 * \brief Transform a vector of Data in a vector of lines with a specific
 *length.
 *
 * \param data vector containing the entries written to the section
 * \param linewidth Maximum length of the lines
 * \return Vector with strings of each line.
 */
std::vector<std::string> makeSection(std::vector<std::string> data,
                                     const int linewidth = 72) {
  std::vector<std::string> out;
  out.reserve(data.size());

  std::string line = "";
  const int number_of_entries = data.size();
  for (auto i = 0; i < number_of_entries; ++i) {
    std::string item = data[i];
    if (i < number_of_entries - 1)
      item += ",";
    else
      item += ";";

    const int line_length = line.length() + item.length();
    if (line_length > linewidth) {
      out.push_back(line);
      line = "";
    }
    line += item;
  }
  out.push_back(line);
  return out;
}

/**
 * \ingroup Geometry
 * \brief Write the IGES header into the given file.
 * 
 * \param file_name File name to write the header in.
 */
void writeIGSHeader(std::string file_name) {
  std::vector<std::string> section(4);
  section[0] = "";
  section[1] = "IGES obtained from Bembel.";
  section[2] = "See <http://www.bembel.eu>";
  section[3] = "";

  writeSection(file_name, section, 'S');
  return;
}

/**
 * \brief Convert a double to a string with given precision.
 * 
 * \param d double to be converted.
 * \param precision Precision of the conversion.
 */
std::string double_to_string(double d, const int precision) {
  std::ostringstream stm;
  stm << std::setprecision(precision) << d;
  return stm.str();
}

/** * \ingroup Geometry
 * \brief Write the IGES Global section into the given file.
 *
 * \param file_name File name to write the section in.
 * \return Number of lines of this section.
 */
int writeGlobalSection(std::string file_name) {
  std::vector<std::string> out(24);
  std::time_t now = std::time(nullptr);

  std::tm* localTime = std::localtime(&now);
  char dateString[11];  // "YYYY-MM-DD\0"
  std::strftime(dateString, sizeof(dateString), "%Y-%m-%d", localTime);
  // Parameter Deliminator Character
  out[0] = "1H,";
  // Record Delimiter Character
  out[1] = "1H;";
  // Product ID from Sender
  out[2] = "Bembel";
  // File Name
  out[3] = file_name;
  // System ID
  out[4] = "Bembel";
  // Pre-processor Version
  out[5] = "writeIGSFile";
  // Number of Bits for Integers
  out[6] = "32";
  // Single Precision Magnitude
  out[7] = "75";
  // Single Precision Significance
  out[8] = "6";
  // Double Precision Magnitude
  out[9] = "75";
  // Double Precision Significance
  out[10] = "15";
  // Product ID for Receiver
  out[11] = "Nurbs from Bembel";
  // Model Space Scale
  out[12] = "1.0";
  // Unit Flag (6 = metres)
  out[13] = "6";
  // Units  (metres = "M")
  out[14] = "M";
  // Maximum Number of Line Weights
  out[15] = "1000";
  // Size of Maximum Line Width
  out[16] = "1.0";
  // Date and Time of file generation
  out[17] = dateString;
  // Minimum User-intended Resolution
  out[18] = "0.000001";
  // Approximate Maximum Coordinate
  out[19] = "10000.0";
  // Name of Author
  out[20] = "Maximilian Nolte";
  // Author's Organization
  out[21] = "CEM - TU Darmstadt";
  // IGES Version Number (3 = IGES version 2.0)
  out[22] = "3";
  // Drafting Standard Code (0 = no standard)
  out[23] = "0";

  std::vector<std::string> section = makeSection(out);

  writeSection(file_name, section, 'G');
  return section.size();
}

/**
 * \ingroup Geometry
 * \brief Write the IGES Directory section into the given file.
 *
 * \param file_name File name to write the section in.
 * \param start_idx Vector containing the start lines of the sections.
 * \param number_of_lines Vector containing the lengths of the sections.
 * \return Number of lines of this section.
 */
int writeDirectory(std::string file_name, std::vector<int> start_idx,
                    std::vector<int> number_of_lines) {
  assert(start_idx.size() == number_of_lines.size());

  std::ofstream fid(file_name, std::ios::app);
  int tot_number_of_lines = 0;
  const int number_of_patches = start_idx.size();
  for (auto i = 0; i < number_of_patches; ++i) {
    fid << std::setw(8) << "128"
        << std::setw(8) << start_idx[i]
        << std::setw(8) << 0
        << std::setw(8) << 0
        << std::setw(8) << 0
        << std::setw(8) << 0
        << std::setw(8) << 0
        << std::setw(8) << 0
        << std::setw(8) << 0
        << "D"
        << std::setw(7) << (i + 1) * 2 - 1<< "\n";

    fid << std::setw(8) << "128"
        << std::setw(8) << 0
        << std::setw(8) << 0
        << std::setw(8) << number_of_lines[i]
        << std::setw(8) << 0
        << std::setw(8) << " "
        << std::setw(8) << " "
        << std::setw(8) << " "
        << std::setw(8) << 0
        << "D"
        << std::setw(7) << (i + 1) * 2 << "\n";

    tot_number_of_lines += 2;
  }
  return tot_number_of_lines;
}

/**
 * \ingroup Geometry
 * \brief This Function writes the patch from Bembel into a vector which can be
 *written into an IGES file.
 *
 * This function assumes that the patch knot vectors do not contains internal
 * nodes.
 * 
 * \param patch Bembel patch containing all information
 * \param precision Number of significant digits to save floats
 * \return Vector with patch data in IGES format.
 */
std::vector<std::string> writePatchData(const Patch& patch,
                                        const int precision) {
  const int degree_x = patch.polynomial_degree_x_;
  const int degree_y = patch.polynomial_degree_y_;

  assert(patch.unique_knots_x_.size() == 2 &&
         "We can not handle internal knots");
  assert(patch.unique_knots_y_.size() == 2 &&
         "We can not handle internal knots");

  std::vector<double> knots_x(2 * degree_x, 0);
  for (auto i = 0; i <= degree_x; ++i) knots_x[degree_x + i] = 1;
  std::vector<double> knots_y(2 * degree_y, 0);
  for (auto i = 0; i <= degree_y; ++i) knots_y[degree_y + i] = 1;

  const std::array<double, 2> span_x = {0, 1};
  const std::array<double, 2> span_y = {0, 1};

  const int number_of_points_x = knots_x.size() - degree_x;
  const int number_of_points_y = knots_y.size() - degree_y;

  const int size = patch.data_.size() / 4;
  std::vector<double> weights(size);
  std::vector<double> coordinates_x(size);
  std::vector<double> coordinates_y(size);
  std::vector<double> coordinates_z(size);
  // transform from wx, wy, wz to x, y, z
  // further more switch row major and column major
  for (auto i = 0; i < size; ++i) {
    const int rowIndex = i % number_of_points_y;
    const int colIndex = i / number_of_points_y;

    // Calculate the index for the corresponding element in row-major order
    const int j = rowIndex * number_of_points_x + colIndex;

    double weight = patch.data_[i * 4 + 3];
    weights[j] = weight;
    coordinates_x[j] = patch.data_[i * 4] / weight;
    coordinates_y[j] = patch.data_[i * 4 + 1] / weight;
    coordinates_z[j] = patch.data_[i * 4 + 2] / weight;
  }

  const int number_of_data_entries = 16;
  const int size_data =
      number_of_data_entries + 4 * size + knots_x.size() + knots_y.size();
  std::vector<std::string> patch_data(size_data);

  patch_data[0] = "128";
  patch_data[1] = std::to_string(number_of_points_x - 1);
  patch_data[2] = std::to_string(number_of_points_y - 1);
  patch_data[3] = std::to_string(degree_x - 1);
  patch_data[4] = std::to_string(degree_y - 1);
  patch_data[5] = "0";
  patch_data[6] = "0";
  patch_data[7] = "0";
  patch_data[8] = "0";
  patch_data[9] = "0";

  for (auto i = 0; i < knots_x.size(); ++i) {
    patch_data[10 + i] = double_to_string(knots_x[i], precision);
  }

  const int idx_knots_y = 10 + knots_x.size();
  for (auto i = 0; i < knots_y.size(); ++i) {
    patch_data[idx_knots_y + i] = double_to_string(knots_y[i], precision);
  }

  const int idx_weights = 10 + knots_x.size() + knots_y.size();
  for (auto i = 0; i < size; ++i) {
    patch_data[idx_weights + i] = double_to_string(weights[i], precision);

    patch_data[idx_weights + size + i * 3] =
        double_to_string(coordinates_x[i], precision);
    patch_data[idx_weights + size + i * 3 + 1] =
        double_to_string(coordinates_y[i], precision);
    patch_data[idx_weights + size + i * 3 + 2] =
        double_to_string(coordinates_z[i], precision);
  }
  patch_data[idx_weights + 4 * size] = std::to_string(span_x[0]);
  patch_data[idx_weights + 4 * size + 1] = std::to_string(span_x[1]);
  patch_data[idx_weights + 4 * size + 2] = std::to_string(span_y[0]);
  patch_data[idx_weights + 4 * size + 3] = std::to_string(span_y[1]);
  patch_data[idx_weights + 4 * size + 4] = "0";
  patch_data[idx_weights + 4 * size + 5] = "0";

  return patch_data;
}

/**
 * \ingroup Geometry
 * \brief Writes Geometry into an IGES file format.
 *
 * \param geometry PatchVector which is written to the file.
 * \param file_name File name to write to.
 * \param precision Significant number of digits for writing floats.
 */
void writeIGSFile(const std::vector<Patch>& geometry, std::string file_name,
                  const int precision = 6) {
  std::ofstream out(file_name);
  out.close();
  writeIGSHeader(file_name);
  const int size_global = writeGlobalSection(file_name);

  const int number_of_patches = geometry.size();
  std::vector<int> start_idx(number_of_patches);
  std::vector<int> end_idx(number_of_patches);
  std::vector<std::vector<std::string>> patch_sections(number_of_patches);
  int start_line = 1;
  for (auto i = 0; i < number_of_patches; ++i) {
    const std::vector<std::string> patch_data =
        writePatchData(geometry[i], precision);
    patch_sections[i] = makeSection(patch_data, 64);

    start_idx[i] = start_line;
    start_line += patch_sections[i].size();
    end_idx[i] = patch_sections[i].size();
  }
  const int size_directory = writeDirectory(file_name, start_idx, end_idx);

  int first_line = 1;
  int size_parameter = 0;
  for (auto it = patch_sections.begin(); it != patch_sections.end(); ++it) {
    const int i = std::distance(patch_sections.begin(), it);
    const int last_idx =
        writeParameterSection(file_name, *it, 1 + 2 * i, first_line);
    first_line = last_idx;
    size_parameter += (*it).size();
  }

  // Terminate Section
  std::stringstream last_section;
  last_section << std::setw(7) << 4 << "S"
               << std::setw(7) << size_global << "G"
               << std::setw(7) << size_directory << "D"
               << std::setw(7) << size_parameter << "P";


  std::ofstream file(file_name, std::ios::app);
  file << std::left << std::setw(72) << last_section.str()
       << "T"
       << std::right << std::setw(7) << 1 << "\n";
  file.close();

  return;
}
}  // namespace Bembel
#endif  // BEMBEL_SRC_GEOMETRY_GEOMETRYIGS_HPP_
