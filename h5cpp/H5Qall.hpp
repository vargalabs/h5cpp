/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 *
 * Bounded lock-free queue and ring-buffer arena — amalgamated for h5cpp.
 * Ported from sigma/include/sigma/{queue,ring}.hpp under MIT relicense.
 *
 * Requires C++20: std::bit_ceil, std::atomic::wait/notify, std::stop_token.
 *
 * Namespaces provided:
 *   bounded::spsc::queue_t<T,N>  — 1 producer, 1 consumer
 *   bounded::mpsc::queue_t<T,N>  — M producers, 1 consumer  (Vyukov)
 *   bounded::spmc::queue_t<T,N>  — 1 producer, M consumers  (Vyukov)
 *   bounded::mpmc::queue_t<T,N>  — M producers, M consumers (Vyukov)
 *   bounded::ring::adaptor_t     — byte-arena ring on top of any queue_t
 *   bounded::ring::spsc_t<N,K>   — spsc ring alias (N ctrl slots, K bytes)
 *   bounded::ring::mpsc_t<N,K>   — mpsc ring alias
 *   sigma::doorbell_t            — standalone notify primitive
 *
 * h5::impl aliases (h5-scoped names for the I/O layer):
 *   h5::impl::spsc_queue_t<T,N>
 *   h5::impl::mpsc_queue_t<T,N>
 *   h5::impl::spmc_queue_t<T,N>
 *   h5::impl::mpmc_queue_t<T,N>
 *   h5::impl::ring_spsc_t<N,K>
 *   h5::impl::ring_mpsc_t<N,K>
 *   h5::impl::staging_ring_t<NP,NC,N,K>  — selects spsc or mpsc by concurrency
 */
#pragma once

#if __cplusplus < 202002L
#  error "H5Qall.hpp requires C++20 or later"
#endif

#include <array>
#include <atomic>
#include <bit>          // std::bit_ceil — C++20
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <new>
#include <optional>
#include <stop_token>   // std::stop_token — C++20
#include <type_traits>
#include <utility>

// ============================================================================
// bounded::impl — internal slot / cell / alignment helpers
// ============================================================================
namespace bounded::impl {

    // hardware_destructive_interference_size in libstdc++ is gated on
    // __GCC_DESTRUCTIVE_SIZE, which GCC always defines but Clang only sets from
    // version 19 onward.  Apple Clang (libc++) derives the constant independently
    // and is unaffected.  Fall back to 64 — correct for x86-64 and Apple Silicon.
    //
    // GCC emits -Winterference-size when the constant is used in a header because
    // the value is technically target-tuning-dependent.  The queue internals are
    // never part of a cross-TU ABI boundary, so the warning is not actionable here.
#ifdef __cpp_lib_hardware_interference_size
#  if defined(__GNUC__) && !defined(__clang__)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Winterference-size"
#  endif
    inline constexpr std::size_t cache_line_size_v =
        std::hardware_destructive_interference_size;
#  if defined(__GNUC__) && !defined(__clang__)
#    pragma GCC diagnostic pop
#  endif
#else
    inline constexpr std::size_t cache_line_size_v = 64;
#endif

    constexpr std::size_t ceil_pow2(std::size_t x) noexcept {
        if (x <= 1) return 1;
        return std::bit_ceil(x);
    }

    template<typename T>
    struct slot_t {
        T* ptr() noexcept {
            return std::launder(reinterpret_cast<T*>(&storage[0]));
        }
        const T* ptr() const noexcept {
            return std::launder(reinterpret_cast<const T*>(&storage[0]));
        }
        template<typename... args_t>
        void emplace(args_t&&... args)
            noexcept(std::is_nothrow_constructible_v<T, args_t...>)
        {
            ::new (static_cast<void*>(&storage[0])) T(std::forward<args_t>(args)...);
        }
        void destroy() noexcept(std::is_nothrow_destructible_v<T>) {
            std::destroy_at(ptr());
        }
        alignas(T) std::byte storage[sizeof(T)];
    };

    constexpr std::size_t align_up(std::size_t n, std::size_t a) noexcept {
        return (n + (a - 1)) & ~(a - 1);
    }

    template<typename T>
    struct alignas(cache_line_size_v) cell_t {
        std::atomic<std::size_t> seq{0};
        slot_t<T> slot;
        static constexpr std::size_t payload_size_v =
            sizeof(std::atomic<std::size_t>) + sizeof(slot_t<T>);
        static constexpr std::size_t padded_size_v =
            align_up(payload_size_v, cache_line_size_v);
        static constexpr std::size_t pad_size_v =
            padded_size_v - payload_size_v;
        [[no_unique_address]] std::byte pad[pad_size_v ? pad_size_v : 1];
    };

    static_assert(alignof(cell_t<int>) >= cache_line_size_v);

} // namespace bounded::impl

// ============================================================================
// sigma::doorbell_t — standalone atomic notify primitive
// ============================================================================
namespace sigma {

    struct alignas(bounded::impl::cache_line_size_v) doorbell_t {
        void ring() noexcept {
            seq.fetch_add(1, std::memory_order_release);
            seq.notify_one();
        }
        void ring_all() noexcept {
            seq.fetch_add(1, std::memory_order_release);
            seq.notify_all();
        }
        std::atomic<std::uint32_t> seq{0};
    };

} // namespace sigma

// ============================================================================
// bounded::spsc — single-producer, single-consumer queue
// ============================================================================
namespace bounded::spsc {

    template<typename T, std::size_t capacity>
    struct queue_t {
        using cache_line_t =
            std::integral_constant<std::size_t,
                bounded::impl::cache_line_size_v>;

        queue_t() = default;
        queue_t(const queue_t&) = delete;
        queue_t& operator=(const queue_t&) = delete;

        ~queue_t() {
            auto t = tail.load(std::memory_order_relaxed);
            const auto h = head.load(std::memory_order_relaxed);
            while (t != h) {
                buffer[t].destroy();
                t = (t + 1) & mask;
            }
        }

        static constexpr std::size_t raw_capacity()    noexcept { return capacity_pow2; }
        static constexpr std::size_t usable_capacity() noexcept { return capacity_pow2 - 1; }

        bool empty() const noexcept {
            const auto t = tail.load(std::memory_order_relaxed);
            const auto h = head.load(std::memory_order_acquire);
            return t == h;
        }
        bool full() const noexcept {
            const auto h    = head.load(std::memory_order_relaxed);
            const auto next = (h + 1) & mask;
            const auto t    = tail.load(std::memory_order_acquire);
            return next == t;
        }
        std::size_t approx_size() const noexcept {
            const auto t = tail.load(std::memory_order_acquire);
            const auto h = head.load(std::memory_order_acquire);
            return (h - t) & mask;
        }

        bool push(const T& v) noexcept(std::is_nothrow_copy_constructible_v<T>) {
            return emplace(v);
        }
        bool push(T&& v) noexcept(std::is_nothrow_move_constructible_v<T>) {
            return emplace(std::move(v));
        }
        template<typename... args_t>
        bool emplace(args_t&&... args)
            noexcept(std::is_nothrow_constructible_v<T, args_t...>)
        {
            const auto h    = head.load(std::memory_order_relaxed);
            const auto next = (h + 1) & mask;
            if (next == tail.load(std::memory_order_acquire)) return false;
            buffer[h].emplace(std::forward<args_t>(args)...);
            head.store(next, std::memory_order_release);
            seq.fetch_add(1, std::memory_order_release);
            seq.notify_one();
            return true;
        }
        bool pop(T& out)
            noexcept(std::is_nothrow_move_assignable_v<T> &&
                     std::is_nothrow_move_constructible_v<T>)
        {
            const auto t = tail.load(std::memory_order_relaxed);
            if (t == head.load(std::memory_order_acquire)) return false;
            out = std::move(*buffer[t].ptr());
            buffer[t].destroy();
            tail.store((t + 1) & mask, std::memory_order_release);
            return true;
        }
        bool wait_pop(T& out, std::stop_token st)
            noexcept(std::is_nothrow_move_assignable_v<T> &&
                     std::is_nothrow_move_constructible_v<T>)
        {
            std::stop_callback on_stop(st, [this]() noexcept { notify_all(); });
            auto last_seq = sequence();
            while (!st.stop_requested()) {
                if (pop(out)) return true;
                if (!wait_for_data(st, last_seq)) return false;
            }
            return false;
        }
        std::uint32_t sequence() const noexcept {
            return seq.load(std::memory_order_acquire);
        }
        void notify_all() noexcept {
            seq.fetch_add(1, std::memory_order_release);
            seq.notify_all();
        }
        bool wait_for_data(std::stop_token st, std::uint32_t& last_seq) const noexcept {
            if (st.stop_requested()) return false;
            if (!empty()) return true;
            seq.wait(last_seq, std::memory_order_acquire);
            last_seq = seq.load(std::memory_order_acquire);
            return !st.stop_requested();
        }

        static constexpr std::size_t capacity_pow2 = impl::ceil_pow2(capacity);
        static_assert(capacity_pow2 >= 2, "capacity must be >= 2");
        static constexpr std::size_t mask = capacity_pow2 - 1;

        alignas(cache_line_t::value) std::atomic<std::size_t>  head{0};
        alignas(cache_line_t::value) std::atomic<std::size_t>  tail{0};
        alignas(cache_line_t::value) std::atomic<std::uint32_t> seq{0};
        std::array<impl::slot_t<T>, capacity_pow2> buffer;
    };

} // namespace bounded::spsc

// ============================================================================
// bounded::mpsc — many-producer, single-consumer (Vyukov)
// ============================================================================
namespace bounded::mpsc {

    template<typename T, std::size_t capacity>
    struct queue_t {
        using cache_line_t =
            std::integral_constant<std::size_t,
                bounded::impl::cache_line_size_v>;

        queue_t() noexcept {
            for (std::size_t i = 0; i < capacity_pow2; ++i)
                buffer[i].seq.store(i, std::memory_order_relaxed);
        }
        queue_t(const queue_t&) = delete;
        queue_t& operator=(const queue_t&) = delete;
        ~queue_t() {
            T tmp;
            while (pop(tmp)) {}
        }

        static constexpr std::size_t raw_capacity()    noexcept { return capacity_pow2; }
        static constexpr std::size_t usable_capacity() noexcept { return capacity_pow2; }

        bool empty_approx() const noexcept {
            return tail.load(std::memory_order_acquire) ==
                   head.load(std::memory_order_acquire);
        }
        std::size_t approx_size() const noexcept {
            const auto t = tail.load(std::memory_order_acquire);
            const auto h = head.load(std::memory_order_acquire);
            return (h >= t) ? (h - t) : 0;
        }

        bool push(const T& v) noexcept(std::is_nothrow_copy_constructible_v<T>) {
            return emplace(v);
        }
        bool push(T&& v) noexcept(std::is_nothrow_move_constructible_v<T>) {
            return emplace(std::move(v));
        }
        template<typename... args_t>
        bool emplace(args_t&&... args)
            noexcept(std::is_nothrow_constructible_v<T, args_t...>)
        {
            std::size_t pos = head.load(std::memory_order_relaxed);
            while (true) {
                auto& cell      = buffer[pos & mask];
                const auto cseq = cell.seq.load(std::memory_order_acquire);
                const auto dif  = static_cast<std::ptrdiff_t>(cseq) -
                                  static_cast<std::ptrdiff_t>(pos);
                if (dif == 0) {
                    if (head.compare_exchange_weak(pos, pos + 1,
                            std::memory_order_relaxed, std::memory_order_relaxed)) {
                        cell.slot.emplace(std::forward<args_t>(args)...);
                        cell.seq.store(pos + 1, std::memory_order_release);
                        doorbell.fetch_add(1, std::memory_order_release);
                        doorbell.notify_one();
                        return true;
                    }
                } else if (dif < 0) {
                    return false; // full
                } else {
                    pos = head.load(std::memory_order_relaxed);
                }
            }
        }
        bool pop(T& out)
            noexcept(std::is_nothrow_move_assignable_v<T> &&
                     std::is_nothrow_move_constructible_v<T>)
        {
            const auto pos  = tail.load(std::memory_order_relaxed);
            auto& cell      = buffer[pos & mask];
            const auto cseq = cell.seq.load(std::memory_order_acquire);
            const auto dif  = static_cast<std::ptrdiff_t>(cseq) -
                              static_cast<std::ptrdiff_t>(pos + 1);
            if (dif < 0) return false;
            out = std::move(*cell.slot.ptr());
            cell.slot.destroy();
            cell.seq.store(pos + capacity_pow2, std::memory_order_release);
            tail.store(pos + 1, std::memory_order_release);
            return true;
        }
        bool wait_pop(T& out, std::stop_token st)
            noexcept(std::is_nothrow_move_assignable_v<T> &&
                     std::is_nothrow_move_constructible_v<T>)
        {
            std::stop_callback on_stop(st, [this]() noexcept { notify_all(); });
            auto last = sequence();
            while (!st.stop_requested()) {
                if (pop(out)) return true;
                if (!wait_for_data(st, last)) return false;
            }
            return false;
        }
        std::uint32_t sequence() const noexcept {
            return doorbell.load(std::memory_order_acquire);
        }
        void notify_all() noexcept {
            doorbell.fetch_add(1, std::memory_order_release);
            doorbell.notify_all();
        }
        bool wait_for_data(std::stop_token st, std::uint32_t& last) const noexcept {
            if (st.stop_requested()) return false;
            if (doorbell.load(std::memory_order_acquire) != last) return true;
            doorbell.wait(last, std::memory_order_acquire);
            last = doorbell.load(std::memory_order_acquire);
            return !st.stop_requested();
        }

        static constexpr std::size_t capacity_pow2 = impl::ceil_pow2(capacity);
        static_assert(capacity_pow2 >= 2, "capacity must be >= 2");
        static constexpr std::size_t mask = capacity_pow2 - 1;

        alignas(cache_line_t::value) std::atomic<std::size_t>  head{0};
        alignas(cache_line_t::value) std::atomic<std::size_t>  tail{0};
        alignas(cache_line_t::value) std::atomic<std::uint32_t> doorbell{0};
        std::array<impl::cell_t<T>, capacity_pow2> buffer;
    };

} // namespace bounded::mpsc

// ============================================================================
// bounded::spmc — single-producer, many-consumer (Vyukov)
// ============================================================================
namespace bounded::spmc {

    template<typename T, std::size_t capacity>
    struct queue_t {
        using cache_line_t =
            std::integral_constant<std::size_t,
                bounded::impl::cache_line_size_v>;

        queue_t() noexcept {
            for (std::size_t i = 0; i < capacity_pow2; ++i)
                buffer[i].seq.store(i, std::memory_order_relaxed);
        }
        queue_t(const queue_t&) = delete;
        queue_t& operator=(const queue_t&) = delete;
        ~queue_t() {
            T tmp;
            while (pop(tmp)) {}
        }

        static constexpr std::size_t raw_capacity()    noexcept { return capacity_pow2; }
        static constexpr std::size_t usable_capacity() noexcept { return capacity_pow2; }

        bool push(const T& v) noexcept(std::is_nothrow_copy_constructible_v<T>) {
            return emplace(v);
        }
        bool push(T&& v) noexcept(std::is_nothrow_move_constructible_v<T>) {
            return emplace(std::move(v));
        }
        template<typename... args_t>
        bool emplace(args_t&&... args)
            noexcept(std::is_nothrow_constructible_v<T, args_t...>)
        {
            const auto pos  = head.load(std::memory_order_relaxed);
            auto& cell      = buffer[pos & mask];
            const auto cseq = cell.seq.load(std::memory_order_acquire);
            const auto dif  = static_cast<std::ptrdiff_t>(cseq) -
                              static_cast<std::ptrdiff_t>(pos);
            if (dif < 0) return false;
            cell.slot.emplace(std::forward<args_t>(args)...);
            cell.seq.store(pos + 1, std::memory_order_release);
            head.store(pos + 1, std::memory_order_release);
            doorbell.fetch_add(1, std::memory_order_release);
            doorbell.notify_one();
            return true;
        }
        bool pop(T& out)
            noexcept(std::is_nothrow_move_assignable_v<T> &&
                     std::is_nothrow_move_constructible_v<T>)
        {
            std::size_t pos = tail.load(std::memory_order_relaxed);
            while (true) {
                auto& cell      = buffer[pos & mask];
                const auto cseq = cell.seq.load(std::memory_order_acquire);
                const auto dif  = static_cast<std::ptrdiff_t>(cseq) -
                                  static_cast<std::ptrdiff_t>(pos + 1);
                if (dif < 0) return false;
                if (dif == 0) {
                    if (tail.compare_exchange_weak(pos, pos + 1,
                            std::memory_order_relaxed, std::memory_order_relaxed)) {
                        out = std::move(*cell.slot.ptr());
                        cell.slot.destroy();
                        cell.seq.store(pos + capacity_pow2, std::memory_order_release);
                        return true;
                    }
                } else {
                    pos = tail.load(std::memory_order_relaxed);
                }
            }
        }
        bool wait_pop(T& out, std::stop_token st)
            noexcept(std::is_nothrow_move_assignable_v<T> &&
                     std::is_nothrow_move_constructible_v<T>)
        {
            std::stop_callback on_stop(st, [this]() noexcept { notify_all(); });
            auto last_seq = sequence();
            while (!st.stop_requested()) {
                if (pop(out)) return true;
                if (!wait_for_data(st, last_seq)) return false;
            }
            return false;
        }
        std::uint32_t sequence() const noexcept {
            return doorbell.load(std::memory_order_acquire);
        }
        void notify_all() noexcept {
            doorbell.fetch_add(1, std::memory_order_release);
            doorbell.notify_all();
        }
        bool wait_for_data(std::stop_token st, std::uint32_t& last_seq) const noexcept {
            if (st.stop_requested()) return false;
            if (doorbell.load(std::memory_order_acquire) != last_seq) return true;
            doorbell.wait(last_seq, std::memory_order_acquire);
            last_seq = doorbell.load(std::memory_order_acquire);
            return !st.stop_requested();
        }

        static constexpr std::size_t capacity_pow2 = impl::ceil_pow2(capacity);
        static_assert(capacity_pow2 >= 2, "capacity must be >= 2");
        static constexpr std::size_t mask = capacity_pow2 - 1;

        alignas(cache_line_t::value) std::atomic<std::size_t>  head{0};
        alignas(cache_line_t::value) std::atomic<std::size_t>  tail{0};
        alignas(cache_line_t::value) std::atomic<std::uint32_t> doorbell{0};
        std::array<impl::cell_t<T>, capacity_pow2> buffer;
    };

} // namespace bounded::spmc

// ============================================================================
// bounded::mpmc — many-producer, many-consumer (Vyukov)
// ============================================================================
namespace bounded::mpmc {

    template<typename T, std::size_t capacity>
    struct queue_t {
        using cache_line_t =
            std::integral_constant<std::size_t,
                bounded::impl::cache_line_size_v>;

        queue_t() noexcept {
            for (std::size_t i = 0; i < capacity_pow2; ++i)
                buffer[i].seq.store(i, std::memory_order_relaxed);
        }
        queue_t(const queue_t&) = delete;
        queue_t& operator=(const queue_t&) = delete;
        ~queue_t() {
            T tmp;
            while (try_pop(tmp)) {}
        }

        static constexpr std::size_t raw_capacity()    noexcept { return capacity_pow2; }
        static constexpr std::size_t usable_capacity() noexcept { return capacity_pow2; }

        bool try_push(const T& v) noexcept(std::is_nothrow_copy_constructible_v<T>) {
            return try_emplace(v);
        }
        bool try_push(T&& v) noexcept(std::is_nothrow_move_constructible_v<T>) {
            return try_emplace(std::move(v));
        }
        template<typename... args_t>
        bool try_emplace(args_t&&... args)
            noexcept(std::is_nothrow_constructible_v<T, args_t...>)
        {
            std::size_t pos = head.load(std::memory_order_relaxed);
            while (true) {
                auto& cell      = buffer[pos & mask];
                const auto cseq = cell.seq.load(std::memory_order_acquire);
                const auto dif  = static_cast<std::ptrdiff_t>(cseq) -
                                  static_cast<std::ptrdiff_t>(pos);
                if (dif == 0) {
                    if (head.compare_exchange_weak(pos, pos + 1,
                            std::memory_order_relaxed, std::memory_order_relaxed)) {
                        cell.slot.emplace(std::forward<args_t>(args)...);
                        cell.seq.store(pos + 1, std::memory_order_release);
                        doorbell.fetch_add(1, std::memory_order_release);
                        doorbell.notify_one();
                        return true;
                    }
                } else if (dif < 0) {
                    return false;
                } else {
                    pos = head.load(std::memory_order_relaxed);
                }
            }
        }
        bool try_pop(T& out)
            noexcept(std::is_nothrow_move_assignable_v<T> &&
                     std::is_nothrow_move_constructible_v<T>)
        {
            std::size_t pos = tail.load(std::memory_order_relaxed);
            while (true) {
                auto& cell      = buffer[pos & mask];
                const auto cseq = cell.seq.load(std::memory_order_acquire);
                const auto dif  = static_cast<std::ptrdiff_t>(cseq) -
                                  static_cast<std::ptrdiff_t>(pos + 1);
                if (dif == 0) {
                    if (tail.compare_exchange_weak(pos, pos + 1,
                            std::memory_order_relaxed, std::memory_order_relaxed)) {
                        out = std::move(*cell.slot.ptr());
                        cell.slot.destroy();
                        cell.seq.store(pos + capacity_pow2, std::memory_order_release);
                        return true;
                    }
                } else if (dif < 0) {
                    return false;
                } else {
                    pos = tail.load(std::memory_order_relaxed);
                }
            }
        }
        bool wait_pop(T& out, std::stop_token st)
            noexcept(std::is_nothrow_move_assignable_v<T> &&
                     std::is_nothrow_move_constructible_v<T>)
        {
            std::stop_callback on_stop(st, [this]() noexcept { notify_all(); });
            auto last = sequence();
            while (!st.stop_requested()) {
                if (try_pop(out)) return true;
                if (!wait_for_data(st, last)) return false;
            }
            return false;
        }
        std::uint32_t sequence() const noexcept {
            return doorbell.load(std::memory_order_acquire);
        }
        void notify_all() noexcept {
            doorbell.fetch_add(1, std::memory_order_release);
            doorbell.notify_all();
        }
        bool wait_for_data(std::stop_token st, std::uint32_t& last) const noexcept {
            if (st.stop_requested()) return false;
            if (doorbell.load(std::memory_order_acquire) != last) return true;
            doorbell.wait(last, std::memory_order_acquire);
            last = doorbell.load(std::memory_order_acquire);
            return !st.stop_requested();
        }

        static constexpr std::size_t capacity_pow2 = impl::ceil_pow2(capacity);
        static_assert(capacity_pow2 >= 2, "capacity must be >= 2");
        static constexpr std::size_t mask = capacity_pow2 - 1;

        alignas(cache_line_t::value) std::atomic<std::size_t>  head{0};
        alignas(cache_line_t::value) std::atomic<std::size_t>  tail{0};
        alignas(cache_line_t::value) std::atomic<std::uint32_t> doorbell{0};
        std::array<impl::cell_t<T>, capacity_pow2> buffer;
    };

} // namespace bounded::mpmc

// ============================================================================
// bounded::ring — byte-arena ring buffer on top of any queue_t
//
// Producer: push(data, length, type) — copies bytes into arena, posts ctrl_t
// Consumer: pop() → ctrl_t;  data(ctrl_t) → const byte*;  release(ctrl_t)
//
// Single-producer assumption: head_ is not protected by a CAS. The inner
// queue_t may be SPSC or MPSC depending on the ctrl slot concurrency needed.
// ============================================================================
namespace bounded::ring {

    struct ctrl_t {
        uint16_t type;
        uint16_t length;
        uint32_t position;
    };
    static_assert(sizeof(ctrl_t) == 8);

    constexpr uint32_t align_up(uint32_t x, uint32_t a) noexcept {
        return (x + (a - 1u)) & ~(a - 1u);
    }

    template<typename queue_t, std::size_t n_bytes_v, uint32_t align_v = 16>
    class adaptor_t {
    public:
        static_assert((n_bytes_v & (n_bytes_v - 1)) == 0,
            "n_bytes_v must be a power of two");

        adaptor_t() = default;
        explicit adaptor_t(queue_t q) : queue_(std::move(q)) {}

        bool push(const void* data, uint32_t length, uint16_t type) noexcept {
            if (!data || length == 0 || length > 0xFFFFu) return false;
            uint32_t pos   = 0;
            std::byte* dst = nullptr;
            if (!try_reserve(length, pos, dst)) return false;
            std::memcpy(dst, data, length);
            return queue_.push(ctrl_t{type, static_cast<uint16_t>(length), pos});
        }

        std::optional<ctrl_t> pop() noexcept {
            ctrl_t token;
            if (!queue_.pop(token)) return std::nullopt;
            return token;
        }

        const std::byte* data(const ctrl_t& c) const noexcept {
            return buffer_.data() + c.position;
        }

        void release(const ctrl_t& c) noexcept {
            tail_.store(tail_.load(std::memory_order_relaxed) + c.length,
                        std::memory_order_release);
        }

        uint32_t    capacity()    const noexcept { return static_cast<uint32_t>(n_bytes_v); }
        std::size_t approx_size() const noexcept { return queue_.approx_size(); }

    private:
        bool try_reserve(uint32_t length, uint32_t& out_pos, std::byte*& out_ptr) noexcept {
            const uint32_t head     = head_.load(std::memory_order_relaxed);
            const uint32_t tail     = tail_.load(std::memory_order_acquire);
            const uint32_t used     = head - tail;
            const uint32_t free_    = static_cast<uint32_t>(n_bytes_v) - used;
            const uint32_t ahead    = align_up(head, align_v);
            const uint32_t pad_aln  = ahead - head;
            if (free_ < pad_aln + length) return false;

            const uint32_t mask = static_cast<uint32_t>(n_bytes_v - 1);
            uint32_t pos        = ahead & mask;
            if (pos + length > n_bytes_v) {
                const uint32_t tail_room = static_cast<uint32_t>(n_bytes_v) - pos;
                if (free_ < pad_aln + tail_room + length) return false;
                if (!queue_.push(ctrl_t{0u, static_cast<uint16_t>(tail_room), pos}))
                    return false;
                head_.store(ahead + tail_room, std::memory_order_release);
                return try_reserve(length, out_pos, out_ptr);
            }
            out_pos = pos;
            out_ptr = buffer_.data() + pos;
            head_.store(ahead + length, std::memory_order_release);
            return true;
        }

        queue_t                                    queue_{};
        std::atomic<uint32_t>                      head_{0}, tail_{0};
        alignas(64) std::array<std::byte, n_bytes_v> buffer_{};
    };

    // Convenience aliases
    template<std::size_t ctrl_cap, std::size_t ring_bytes>
    using spsc_t = adaptor_t<bounded::spsc::queue_t<ctrl_t, ctrl_cap>, ring_bytes>;

    template<std::size_t ctrl_cap, std::size_t ring_bytes>
    using mpsc_t = adaptor_t<bounded::mpsc::queue_t<ctrl_t, ctrl_cap>, ring_bytes>;

} // namespace bounded::ring

// ============================================================================
// h5::impl — h5cpp-scoped aliases and thread-topology selector
// ============================================================================
namespace h5::impl {

    template<typename T, std::size_t N>
    using spsc_queue_t = bounded::spsc::queue_t<T, N>;

    template<typename T, std::size_t N>
    using mpsc_queue_t = bounded::mpsc::queue_t<T, N>;

    template<typename T, std::size_t N>
    using spmc_queue_t = bounded::spmc::queue_t<T, N>;

    template<typename T, std::size_t N>
    using mpmc_queue_t = bounded::mpmc::queue_t<T, N>;

    template<std::size_t ctrl_cap, std::size_t ring_bytes>
    using ring_spsc_t = bounded::ring::spsc_t<ctrl_cap, ring_bytes>;

    template<std::size_t ctrl_cap, std::size_t ring_bytes>
    using ring_mpsc_t = bounded::ring::mpsc_t<ctrl_cap, ring_bytes>;

    // Select queue topology from producer/consumer counts.
    // NP == 1 && NC == 1 → spsc (no CAS on either side)
    // NP  > 1 && NC == 1 → mpsc (CAS on producer head)
    // All other combinations fall back to mpsc (consumer always single in h5cpp).
    template<std::size_t NP, std::size_t NC,
             std::size_t ctrl_cap, std::size_t ring_bytes>
    using staging_ring_t = std::conditional_t<
        (NP == 1 && NC == 1),
        bounded::ring::spsc_t<ctrl_cap, ring_bytes>,
        bounded::ring::mpsc_t<ctrl_cap, ring_bytes>>;

} // namespace h5::impl
