#ifndef ATOMICENERGY_TIMER_H
#define ATOMICENERGY_TIMER_H
#include <chrono>
#include <iostream>
#include <sstream>

class Timer {
    using timestamp = std::chrono::time_point<std::chrono::high_resolution_clock>;
    using clock = std::chrono::high_resolution_clock;
private:
    timestamp startTime;
    timestamp endTime;
    bool isRunning = false;
    int timerIdx;

    static int timersCount;

public:
    Timer();

    void Start();
    void Stop();
    std::string ElapsedMilliseconds();
};

#endif //ATOMICENERGY_TIMER_H
