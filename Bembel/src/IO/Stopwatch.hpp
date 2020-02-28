// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_IO_STOPWATCH_H_
#define BEMBEL_IO_STOPWATCH_H_

namespace Bembel {

namespace IO {
/**
 * \ingroup IO
 * \brief A simple class for benchmarking.
 *
 *	A simple stopwatch class for benchmarking. start() starts the timer,
 *	lap() gives the time from either start or the last lap, stop() gives the
 *  overall time since start, and get_data() returns a std::vector<double> with
 *  all lap-times. It is written such that its own timekeeping is not measured.
 */

class Stopwatch {
  std::chrono::time_point<std::chrono::high_resolution_clock> p_;
  double sum_;
  std::vector<double> durations_;
  bool has_been_started_;
  bool has_lapped_;

 public:
  inline Stopwatch() {
    has_been_started_ = false;
    has_lapped_ = false;
  };

  void tic(void) { p_ = std::chrono::high_resolution_clock::now(); }
  double toc(void) {
    double dtime = std::chrono::duration<double, std::ratio<1, 1>>(
                       std::chrono::high_resolution_clock::now() - p_)
                       .count();
    return dtime;
  }
  inline void start() {
    durations_ = {};
    assert(has_been_started_ == false && "Stopwatch already running!");
    has_been_started_ = true;
    p_ = std::chrono::high_resolution_clock::now();
  }
  inline double lap() {
    double out = std::chrono::duration<double, std::ratio<1, 1>>(
                     std::chrono::high_resolution_clock::now() - p_)
                     .count();
    durations_.push_back(out);
    assert(has_been_started_ == true && "Stopwatch must be started first!");
    has_lapped_ = true;
    p_ = std::chrono::high_resolution_clock::now();
    return out;
  }
  inline double stop() {
    double out = std::chrono::duration<double, std::ratio<1, 1>>(
                     std::chrono::high_resolution_clock::now() - p_)
                     .count();
    assert(has_been_started_ == true && "Stopwatch must be started first!");
    if (has_lapped_) {
      has_lapped_ = false;
      has_been_started_ = false;
      return std::accumulate(durations_.begin(), durations_.end(), 0.0);
    }
    return out;
  }
  inline std::vector<double> get_data() { return durations_; }
};
}  // namespace IO
}  // namespace Bembel
#endif
