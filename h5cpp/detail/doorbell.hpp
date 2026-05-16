#pragma once
// Doorbell: a sequenced wake primitive. C++20 uses atomic wait/notify
// (lock-free, maps to OS futex). C++17 falls back to mutex+condvar.
// API is the same either way.
#include <cstdint>
#ifdef __cpp_lib_atomic_wait
#  include <atomic>
namespace h5::detail {
    struct doorbell_t {
        void ring() noexcept {
            seq_.fetch_add(1, std::memory_order_release);
            seq_.notify_one();
        }
        void ring_all() noexcept {
            seq_.fetch_add(1, std::memory_order_release);
            seq_.notify_all();
        }
        void wait(std::uint32_t last) const noexcept {
            seq_.wait(last, std::memory_order_acquire);
        }
        std::uint32_t load() const noexcept {
            return seq_.load(std::memory_order_acquire);
        }
    private:
        std::atomic<std::uint32_t> seq_{0};
    };
} // namespace h5::detail
#else
#  include <condition_variable>
#  include <mutex>
namespace h5::detail {
    struct doorbell_t {
        void ring() noexcept {
            { std::lock_guard<std::mutex> lk(m_); ++seq_; }
            cv_.notify_one();
        }
        void ring_all() noexcept {
            { std::lock_guard<std::mutex> lk(m_); ++seq_; }
            cv_.notify_all();
        }
        void wait(std::uint32_t last) {
            std::unique_lock<std::mutex> lk(m_);
            cv_.wait(lk, [&]{ return seq_ != last; });
        }
        std::uint32_t load() const noexcept {
            std::lock_guard<std::mutex> lk(const_cast<std::mutex&>(m_));
            return seq_;
        }
    private:
        mutable std::mutex m_;
        std::condition_variable cv_;
        std::uint32_t seq_{0};
    };
} // namespace h5::detail
#endif
