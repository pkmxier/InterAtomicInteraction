#include "Timer.h"

int Timer::timersCount = 0;

Timer::Timer() {
    timerIdx = timersCount;
    ++timersCount;
}

void Timer::Start() {
    startTime = clock::now();
    isRunning = true;
}

void Timer::Stop() {
    endTime = clock::now();
    isRunning = false;
}

std::string Timer::ElapsedMilliseconds(){
    timestamp endTime;

    if (isRunning) {
        endTime = clock::now();
    }
    else {
        endTime = this->endTime;
    }

    std::stringstream ss;

    double cnt = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    ss << "Timer{" << timerIdx << "}: elapsed " << cnt << " s";

    return ss.str();
}