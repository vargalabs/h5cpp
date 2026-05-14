// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Temporary diagnostic: localise the SegFault inside h5::write(fd, path, ref)
// on Windows MSVC by manually unrolling every step the high-level overload
// performs internally.  The last successful probe before the crash pins the
// offending call.  Delete this file once root-caused.

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <vector>
#include <string>
#include <hdf5.h>

#include <h5cpp/core>
#include <h5cpp/io>

namespace {
    std::uintptr_t g_sp_top = 0;

    inline std::ptrdiff_t stack_used() {
        volatile char marker;
        std::uintptr_t sp = reinterpret_cast<std::uintptr_t>(const_cast<char*>(&marker));
        if (g_sp_top == 0 || sp > g_sp_top) g_sp_top = sp;
        return static_cast<std::ptrdiff_t>(g_sp_top - sp);
    }

    void probe(const char* label) {
        std::fprintf(stderr, "[STACKPROBE] %-50s used=%lld\n",
                     label, static_cast<long long>(stack_used()));
        std::fflush(stderr);
    }
}

int main() {
    probe("main entry");

    const char* fname = "win_stack_probe.h5";
    std::error_code ec;
    std::filesystem::remove(fname, ec);
    probe("after std::filesystem::remove");

    probe("before h5::create (fd)");
    h5::fd_t fd = h5::create(fname, H5F_ACC_TRUNC);
    probe("after  h5::create (fd)");
    std::fprintf(stderr, "      fd = %lld\n", (long long)static_cast<hid_t>(fd));

    std::vector<double> data{1.0, 2.0, 3.0};

    // ---- Manually unroll h5::write(fd, "vec", data) one step at a time ----
    //
    // 1) error stack mute / H5Lexists / unmute  (the "does it already exist?" probe)
    probe("before h5::mute");
    h5::mute();
    probe("after  h5::mute");

    probe("before H5Lexists");
    htri_t exists = H5Lexists(static_cast<hid_t>(fd), "vec", H5P_DEFAULT);
    probe("after  H5Lexists");
    std::fprintf(stderr, "      H5Lexists -> %d\n", (int)exists);

    probe("before h5::unmute");
    h5::unmute();
    probe("after  h5::unmute");

    // 2) size of payload  (impl::size(ref) -> std::array<size_t,1> -> count_t)
    probe("before impl::size(vec)");
    auto sz = h5::impl::size(data);
    probe("after  impl::size(vec)");
    std::fprintf(stderr, "      sz[0] = %zu\n", sz[0]);

    probe("before h5::count_t count{sz}");
    h5::count_t count = sz;
    probe("after  h5::count_t count{sz}");
    std::fprintf(stderr, "      count.rank=%d count[0]=%llu\n",
                 count.rank, (unsigned long long)count[0]);

    // 3) current_dims derivation
    probe("before h5::impl::get_current_dims(count)");
    h5::current_dims_t current_dims = h5::impl::get_current_dims(count);
    probe("after  h5::impl::get_current_dims(count)");

    // 4) dataset creation (same path basics_io uses and which is proven to work)
    probe("before h5::create<double>(fd,path,current_dims)");
    h5::ds_t ds = h5::create<double>(fd, "vec", current_dims);
    probe("after  h5::create<double>(fd,path,current_dims)");
    std::fprintf(stderr, "      ds = %lld\n", (long long)static_cast<hid_t>(ds));

    // 5) raw pointer to data
    probe("before impl::data(vec)");
    const double* ptr = h5::impl::data(data);
    probe("after  impl::data(vec)");
    std::fprintf(stderr, "      ptr = %p  ptr[0]=%g\n", (void*)ptr, ptr[0]);

    // 6) file_space + rank
    probe("before H5Dget_space(ds)");
    hid_t file_space_id = H5Dget_space(static_cast<hid_t>(ds));
    probe("after  H5Dget_space(ds)");
    std::fprintf(stderr, "      file_space = %lld\n", (long long)file_space_id);

    probe("before H5Sget_simple_extent_ndims");
    int rank = H5Sget_simple_extent_ndims(file_space_id);
    probe("after  H5Sget_simple_extent_ndims");
    std::fprintf(stderr, "      rank = %d\n", rank);

    // 7) mem_space + select_all
    probe("before h5::create_simple(n_elements)");
    h5::sp_t mem_space = h5::create_simple(static_cast<hsize_t>(count));
    probe("after  h5::create_simple(n_elements)");
    std::fprintf(stderr, "      mem_space = %lld\n",
                 (long long)static_cast<hid_t>(mem_space));

    probe("before h5::select_all(mem_space)");
    h5::select_all(mem_space);
    probe("after  h5::select_all(mem_space)");

    probe("before H5Sselect_all(file_space)");
    H5Sselect_all(file_space_id);
    probe("after  H5Sselect_all(file_space)");

    // 8) THE PRIME SUSPECT: default-constructed h5::dt_t<double>
    probe("before h5::dt_t<double>{} default ctor");
    h5::dt_t<double> type;
    probe("after  h5::dt_t<double>{} default ctor");
    std::fprintf(stderr, "      h5::dt_t<double> -> %lld\n",
                 (long long)static_cast<hid_t>(type));

    // 9) raw HDF5 C call — bypass h5cpp wrappers
    probe("before H5Dwrite (raw C call)");
    herr_t werr = H5Dwrite(static_cast<hid_t>(ds),
                           static_cast<hid_t>(type),
                           static_cast<hid_t>(mem_space),
                           file_space_id,
                           H5P_DEFAULT,
                           ptr);
    probe("after  H5Dwrite (raw C call)");
    std::fprintf(stderr, "      H5Dwrite err = %d\n", (int)werr);

    H5Sclose(file_space_id);

    probe("done (unrolled)");
    std::filesystem::remove(fname, ec);
    return werr < 0 ? 1 : 0;
}
