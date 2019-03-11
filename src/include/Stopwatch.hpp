// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _UTLI_C_STOPWATCH_
#define _UTLI_C_STOPWATCH_

#include <assert.h>
#include <chrono>
#include <numeric>
#include <utility>
#include <vector>

namespace Bembel {

namespace Util {
/**
 * @brief A simple class for benchmarking.
 *
 *	A simple stopwatch class for benchmarking. start() starts the timer,
 *	lap() gives the time from either start or the last lap, stop() gives the
 *  overall time since start, and get_data() returns a std::vector<double> with
 *  all lap-times. It is written such that its own timekeeping is not measured.
 */

class Stopwatch {
  std::chrono::time_point<std::chrono::high_resolution_clock> _p;
  double _sum;
  std::vector<double> _durations;
  bool _has_been_started;
  bool _has_lapped;

 public:
  inline Stopwatch() {
    _has_been_started = false;
    _has_lapped = false;
  };
  inline void start() {
    _durations = {};
    assert(_has_been_started == false && "Stopwatch already running!");
    _has_been_started = true;
    _p = std::chrono::high_resolution_clock::now();
  }
  inline double lap() {
    double out = std::chrono::duration<double, std::ratio<1, 1>>(
                     std::chrono::high_resolution_clock::now() - _p)
                     .count();
    _durations.push_back(out);
    assert(_has_been_started == true && "Stopwatch must be started first!");
    _has_lapped = true;
    _p = std::chrono::high_resolution_clock::now();
    return out;
  }
  inline double stop() {
    double out = std::chrono::duration<double, std::ratio<1, 1>>(
                     std::chrono::high_resolution_clock::now() - _p)
                     .count();
    assert(_has_been_started == true && "Stopwatch must be started first!");
    if (_has_lapped) {
      _has_lapped = false;
      _has_been_started = false;
      return std::accumulate(_durations.begin(), _durations.end(), 0.0);
    }
    return out;
  }
  inline std::vector<double> get_data() { return _durations; }
};
}  // namespace Util
}  // namespace Bembel
#endif
