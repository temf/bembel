// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _C_BEMBEL_LOGGER_
#define _C_BEMBEL_LOGGER_


namespace Bembel {
namespace IO {

/**
 * \ingroup IO
 * \brief A simple class for producing output.
 *
 *  A simple class for producing output. Its not meant as a logger in
 * the classical sense, but rather to write measurement data in a machine &
 * human readable output. For example, it can be used to generate files to be
 * read by pgfplots.
 *
 */

template <int N = 10>
class Logger {
  std::string _sep = " ";
  std::fstream _file;
  std::ios_base::openmode _mode = std::fstream::trunc;
  std::string name_ = "default.log";

 public:
  Logger() {}
  Logger(std::string name) { name_ = name; }
  Logger(std::string name, std::string sep) {
    name_ = name;
    _sep = sep;
  }
  Logger(std::string name, std::string sep, std::ios_base::openmode mode) {
    name_ = name;
    _sep = sep;
    _mode = mode;
  }
  Logger(std::string name, std::ios_base::openmode mode) {
    name_ = name;
    _mode = mode;
  }
  Logger(std::ios_base::openmode mode) { _mode = mode; }
  template <typename T, typename... Args>
  void term(T t, Args... arg) {
    std::cout << std::setprecision(N - 4) << std::setw(N + _sep.size())
              << std::left << t << _sep;
    term(arg...);
  }
  template <typename T>
  void term(T t) {
    std::cout << std::setprecision(N - 4) << std::setw(N + _sep.size()) << t
              << std::endl;
  }
  template <typename T, typename... Args>
  void file(T t, Args... arg) {
    if (!(_file.is_open())) {
      _file.open(name_, std::fstream::out | std::fstream::trunc);
    }
    _file << std::setprecision(N - 4) << std::setw(N + _sep.size()) << std::left
          << t << _sep;
    file(arg...);
  }
  template <typename T>
  void file(T t) {
    if (!(_file.is_open())) {
      _file.open(name_, std::fstream::out | std::fstream::trunc);
    }
    std::setprecision(N);
    _file << std::setprecision(N - 4) << std::setw(N + _sep.size()) << std::left
          << t << std::endl;
  }
  template <typename... Args>
  void both(Args... arg) {
    file(arg...);
    term(arg...);
  }
  ~Logger() { _file.close(); }
};
}  // namespace Util
}  // namespace Bembel
#endif
