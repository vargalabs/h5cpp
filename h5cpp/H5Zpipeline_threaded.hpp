/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once
#include <atomic>
#include <cstring>
#include <memory>
#include <thread>
#include <vector>
#include "H5Zpipeline.hpp"
#include "H5Qall.hpp"

// Workers handle compression only; H5Dwrite_chunk is always called from the
// main thread.  Work items use raw ::hid_t integers (not h5cpp RAII wrappers)
// so no HDF5 API calls occur on worker threads, preserving compatibility with
// non-thread-safe HDF5 builds.

namespace h5::impl {

namespace detail {

    // Uncompressed chunk dispatched to a worker for filter application.
    struct raw_work_t {
        std::unique_ptr<std::byte[]> data;
        std::array<hsize_t, H5CPP_MAX_RANK> offset{};
        std::size_t nbytes{0};
        ::hid_t   ds_id{H5I_UNINIT};
        ::hid_t   dxpl_id{H5I_UNINIT};

        raw_work_t() = default;
        raw_work_t(raw_work_t&&) = default;
        raw_work_t& operator=(raw_work_t&&) = default;
        raw_work_t(const raw_work_t&) = delete;
        raw_work_t& operator=(const raw_work_t&) = delete;
    };

    // Filtered chunk ready for H5Dwrite_chunk, returned to main thread.
    struct done_work_t {
        aligned_ptr data;       // compressed/filtered bytes (aligned allocation)
        std::array<hsize_t, H5CPP_MAX_RANK> offset{};
        std::size_t nbytes{0};  // compressed size
        std::uint32_t mask{0};
        ::hid_t   ds_id{H5I_UNINIT};
        ::hid_t   dxpl_id{H5I_UNINIT};

        done_work_t() = default;
        done_work_t(done_work_t&&) = default;
        done_work_t& operator=(done_work_t&&) = default;
        done_work_t(const done_work_t&) = delete;
        done_work_t& operator=(const done_work_t&) = delete;
    };

} // namespace detail

struct threaded_pipeline_t : public pipeline_t<threaded_pipeline_t> {
    threaded_pipeline_t() = default;

    ~threaded_pipeline_t() {
        // Drain all in-flight chunks (workers compress, main thread writes).
        flush();
        // std::jthread dtors call request_stop() then join() automatically.
    }

    threaded_pipeline_t(const threaded_pipeline_t&) = delete;
    threaded_pipeline_t& operator=(const threaded_pipeline_t&) = delete;
    threaded_pipeline_t(threaded_pipeline_t&&) = delete;
    threaded_pipeline_t& operator=(threaded_pipeline_t&&) = delete;

    // Drains all in-flight work and calls H5Dwrite_chunk from the calling
    // (main) thread.  Must be called before reading back written data.
    void flush() {
        while (in_flight_.load(std::memory_order_acquire) > 0) {
            drain_done();
            std::this_thread::yield();
        }
        drain_done(); // final sweep after counter hits zero
    }

    void write_chunk_impl(const hsize_t* offset, std::size_t nbytes, const void* ptr) {
        ensure_workers();
        drain_done(); // opportunistic: write any already-compressed chunks

        detail::raw_work_t work;
        work.data = std::make_unique<std::byte[]>(nbytes);
        std::memcpy(work.data.get(), ptr, nbytes);
        std::copy(offset, offset + this->rank, work.offset.data());
        work.nbytes  = nbytes;
        work.ds_id   = static_cast<::hid_t>(this->ds);
        work.dxpl_id = static_cast<::hid_t>(this->dxpl);

        in_flight_.fetch_add(1, std::memory_order_release);
        while (!compress_queue_.push(std::move(work))) {
            drain_done(); // relieve back-pressure
            std::this_thread::yield();
        }
    }

    void read_chunk_impl(const hsize_t* offset, std::size_t nbytes, void* /*ptr*/) {
        uint32_t filter_mask;
        std::size_t length = nbytes;
        if (tail == 0) {
#if H5_VERSION_GE(2,0,0)
            std::size_t buf_size = nbytes;
            H5Dread_chunk2(ds, dxpl, offset, &filter_mask, chunk0, &buf_size);
#else
            H5Dread_chunk(ds, dxpl, offset, &filter_mask, chunk0);
#endif
            return;
        }
        void* read_target = (tail % 2 == 1) ? chunk1 : chunk0;
#if H5_VERSION_GE(2,0,0)
        std::size_t buf_size = nbytes;
        H5Dread_chunk2(ds, dxpl, offset, &filter_mask, read_target, &buf_size);
#else
        H5Dread_chunk(ds, dxpl, offset, &filter_mask, read_target);
#endif
        void* src = read_target;
        void* dst = (read_target == chunk0)
            ? static_cast<void*>(chunk1) : static_cast<void*>(chunk0);
        for (hsize_t j = tail; j > 0; --j) {
            const hsize_t fi = j - 1;
            length = filter[fi](dst, src, length,
                flags[fi] | H5Z_FLAG_REVERSE, cd_size[fi], cd_values[fi]);
            void* tmp = src; src = dst; dst = tmp;
        }
    }

private:
    // Called from main thread only: pop finished compressed chunks and write.
    void drain_done() {
        detail::done_work_t result;
        while (done_queue_.pop(result)) {
            H5Dwrite_chunk(result.ds_id, result.dxpl_id, result.mask,
                           result.offset.data(), result.nbytes, result.data.get());
            in_flight_.fetch_sub(1, std::memory_order_release);
        }
    }

    void ensure_workers() {
        std::call_once(init_flag_, [this] {
            constexpr unsigned cfg = static_cast<unsigned>(H5CPP_PIPELINE_WORKERS);
            const unsigned n = (cfg > 0) ? cfg
                                         : std::max(1u, std::thread::hardware_concurrency());
            workers_.reserve(n);
            for (unsigned i = 0; i < n; ++i)
                workers_.emplace_back([this](std::stop_token st) { worker_loop(st); });
        });
    }

    void worker_loop(std::stop_token st) {
        const std::size_t scratch = filter::filter_scratch_bound(block_size);
        aligned_ptr wbuf0 = make_aligned(H5CPP_MEM_ALIGNMENT, scratch);
        aligned_ptr wbuf1 = make_aligned(H5CPP_MEM_ALIGNMENT, scratch);

        detail::raw_work_t work;
        while (compress_queue_.wait_pop(work, st)) {
            detail::done_work_t result = compress(work, wbuf0, wbuf1);
            while (!done_queue_.push(std::move(result)))
                std::this_thread::yield();
        }
    }

    detail::done_work_t compress(const detail::raw_work_t& work,
                                 aligned_ptr& wbuf0, aligned_ptr& wbuf1)
    {
        detail::done_work_t out;
        out.offset  = work.offset;
        out.ds_id   = work.ds_id;
        out.dxpl_id = work.dxpl_id;

        std::size_t length = work.nbytes;

        if (tail == 0) {
            // No filters: copy raw bytes into an aligned buffer.
            const std::size_t sz = filter::filter_scratch_bound(length);
            out.data  = make_aligned(H5CPP_MEM_ALIGNMENT, sz);
            std::memcpy(out.data.get(), work.data.get(), length);
            out.nbytes = length;
            out.mask   = 0;
            return out;
        }

        // First filter: work.data → wbuf0.
        length = filter[0](wbuf0.get(), work.data.get(), length,
                           flags[0], cd_size[0], cd_values[0]);
        if (!length) out.mask |= 1u;

        // Subsequent filters: ping-pong between wbuf0 and wbuf1.
        void* src = wbuf0.get();
        void* dst = wbuf1.get();
        for (hsize_t j = 1; j < tail; ++j) {
            length = filter[j](dst, src, length, flags[j], cd_size[j], cd_values[j]);
            if (!length) out.mask |= (1u << j);
            std::swap(src, dst);
        }

        // Copy filtered result into fresh aligned buffer so scratch is reusable.
        const std::size_t sz = filter::filter_scratch_bound(length);
        out.data = make_aligned(H5CPP_MEM_ALIGNMENT, sz);
        std::memcpy(out.data.get(), src, length);
        out.nbytes = length;
        return out;
    }

    // compress_queue_: main thread pushes raw chunks; workers pop and compress.
    bounded::spmc::queue_t<detail::raw_work_t,  64> compress_queue_;
    // done_queue_:     workers push compressed chunks; main thread pops and writes.
    bounded::mpsc::queue_t<detail::done_work_t, 64> done_queue_;

    std::atomic<int>          in_flight_{0};
    std::vector<std::jthread> workers_;
    std::once_flag            init_flag_;
};

} // namespace h5::impl
