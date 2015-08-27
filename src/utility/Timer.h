#ifndef TurboStructured_utility_Timer
#define TurboStructured_utility_Timer

#include <cassert>
#include <chrono>
#include <mutex>
#include <memory>
#include <thread>

using namespace std;
using namespace std::chrono;

//Simple timer object
class Timer {
  using clock_t = high_resolution_clock;
  clock_t::time_point _pauseTime;
public:
	bool IsActive;
  clock_t::time_point StartTime;
  clock_t::time_point EndTime;
  clock_t::duration ElapsedTime;

	Timer() {
		IsActive = false;
		ElapsedTime = clock_t::duration(0);
	};

	double ElapsedTimeMilliseconds() {
		return duration_cast<duration<double, std::milli>>(ElapsedTime).count();
	};

	void Start() {
		IsActive = true;
		StartTime = clock_t::now();
		_pauseTime = StartTime;
		ElapsedTime = clock_t::duration(0);
	};

	void Resume() {
		if (!IsActive) {
			_pauseTime = clock_t::now();
			IsActive = true;
		};
	};

	void Pause() {
		if (IsActive) {
			ElapsedTime += clock_t::now() - _pauseTime;
			IsActive = false;
		};
	};

	void Reset() {
		if (!IsActive) {
			ElapsedTime = clock_t::duration(0);
		};
	};

	void Stop() {
		Pause();
		EndTime = clock_t::now();
	};
};

#endif