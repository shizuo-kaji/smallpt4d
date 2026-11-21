#pragma once

#include <iostream>
#include <mutex>

class ProgressBar {
public:
    ProgressBar(unsigned size, unsigned max_ticks)
        : size_(size), max_ticks_(max_ticks), current_(0) {}

    void end() const {
        std::cout << std::endl;
    }

    void print() {
        std::lock_guard<std::mutex> guard(tick_mutex_);
        const unsigned safe_ticks = (max_ticks_ == 0) ? 1u : max_ticks_;
        const unsigned pos = size_ * current_ / safe_ticks;
        std::cout << "\r[";
        for (unsigned i = 0; i < size_; ++i) {
            if (i < pos) {
                std::cout << "=";
            } else if (i == pos) {
                std::cout << ">";
            } else {
                std::cout << ".";
            }
        }
        const unsigned percent = current_ * 100u / safe_ticks;
        std::cout << "] " << percent << " %";
        std::cout.flush();
    }

    void tick() {
        std::lock_guard<std::mutex> guard(tick_mutex_);
        if (current_ < max_ticks_) {
            ++current_;
        }
    }

private:
    unsigned size_;
    unsigned max_ticks_;
    unsigned current_;
    std::mutex tick_mutex_;
};

