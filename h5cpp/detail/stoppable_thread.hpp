#pragma once
// stoppable_thread_t: a thread that supports cooperative cancellation.
// C++20: thin alias over std::jthread (stop_token injected automatically).
// C++17: std::thread wrapper with atomic stop flag; user lambda receives
//        h5::detail::stop_token_t as first argument.
#include <thread>
#include <atomic>
#include <functional>

#ifdef __cpp_lib_jthread
#  include <stop_token>
namespace h5::detail {
    using stop_token_t       = std::stop_token;
    using stoppable_thread_t = std::jthread;
} // namespace h5::detail

#else
namespace h5::detail {

    class stop_token_t {
        std::shared_ptr<std::atomic<bool>> flag_;
        friend class stoppable_thread_t;
        explicit stop_token_t(std::shared_ptr<std::atomic<bool>> f) noexcept
            : flag_(std::move(f)) {}
    public:
        stop_token_t() = default;
        [[nodiscard]] bool stop_requested() const noexcept {
            return flag_ && flag_->load(std::memory_order_acquire);
        }
    };

    class stoppable_thread_t {
        std::shared_ptr<std::atomic<bool>> flag_ =
            std::make_shared<std::atomic<bool>>(false);
        std::thread t_;
    public:
        stoppable_thread_t() = default;

        template<class Fn, class... Args>
        explicit stoppable_thread_t(Fn&& fn, Args&&... args) {
            stop_token_t st{flag_};
            t_ = std::thread(std::forward<Fn>(fn), std::move(st),
                             std::forward<Args>(args)...);
        }

        ~stoppable_thread_t() {
            if (t_.joinable()) {
                flag_->store(true, std::memory_order_release);
                t_.join();
            }
        }

        stoppable_thread_t(const stoppable_thread_t&) = delete;
        stoppable_thread_t& operator=(const stoppable_thread_t&) = delete;
        stoppable_thread_t(stoppable_thread_t&&) noexcept = default;
        stoppable_thread_t& operator=(stoppable_thread_t&&) noexcept = default;

        void request_stop() noexcept {
            flag_->store(true, std::memory_order_release);
        }
        [[nodiscard]] bool joinable() const noexcept { return t_.joinable(); }
        void join() { if (t_.joinable()) t_.join(); }
    };

} // namespace h5::detail
#endif
