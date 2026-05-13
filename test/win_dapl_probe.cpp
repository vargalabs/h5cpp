// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// TEMP DIAGNOSTIC — wide-net Windows MSVC probes for h5::write segfault.
//
// Hypothesis under test:
//   Every failing wrapper write funnels through
//       h5::get_access_plist(ds)  ->  ds.dapl
//       H5Pexist(ds.dapl, H5CPP_DAPL_HIGH_THROUGHPUT)            (H5Dwrite.hpp:86,89)
//   The one passing wrapper (low_level_explicit_spaces) bypasses that path.
//   Both h5::open (H5Dopen.hpp:56) and h5::createds (H5capi.hpp:191) explicitly
//   assign ds_.dapl, but on MSVC the value may not survive the chain
//   h5::create<T>(fd,path,...) -> createds(...) -> caller.
//   These probes confirm or refute by inspecting ds.dapl directly, by trying a
//   forced fix (ds.dapl = H5P_DEFAULT) and by contrasting create-path vs open-path.
//
// Delete this file once root-caused.

#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>
#include <hdf5.h>

#include <h5cpp/all>

#ifndef H5CPP_WIN_DAPL_PROBE_MODE
#define H5CPP_WIN_DAPL_PROBE_MODE 0
#endif

namespace {
    void probe(const char* label) {
        std::fprintf(stderr, "[WINDAPL] %s\n", label);
        std::fflush(stderr);
    }

    void dump(const char* label, hid_t h) {
        std::fprintf(stderr, "[WINDAPL]   %-26s id=%lld  H5Iis_valid=%d\n",
            label, (long long)h, (int)H5Iis_valid(h));
        std::fflush(stderr);
    }

    std::string filename(const char* suffix) {
        return std::string("win_dapl_") + suffix + ".h5";
    }

    void remove_file(const std::string& path) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }

    // --- 1: dump ds.dapl after h5::create<double> -----------------------------
    int inspect_create_ds_dapl() {
        auto path = filename("inspect_create");
        remove_file(path);
        h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
        probe("inspect_create_ds_dapl: before h5::create<double>");
        h5::ds_t ds = h5::create<double>(fd, "vec", h5::current_dims{3});
        probe("inspect_create_ds_dapl: after  h5::create<double>");
        dump("ds.handle", static_cast<hid_t>(ds));
        dump("ds.dapl",   ds.dapl);
        // Only attempt H5Pexist if ds.dapl looks plausible.
        if (ds.dapl == H5P_DEFAULT) {
            probe("ds.dapl == H5P_DEFAULT (no-op for H5Pexist)");
        } else if (H5Iis_valid(ds.dapl) > 0) {
            int r = H5Pexist(ds.dapl, H5CPP_DAPL_HIGH_THROUGHPUT);
            std::fprintf(stderr, "[WINDAPL]   H5Pexist(ds.dapl,...) -> %d\n", r);
        } else {
            probe("ds.dapl INVALID — this is the segfault trigger");
        }
        remove_file(path);
        return 0;
    }

    // --- 2: dump ds.dapl after h5::open (control: should be valid) ------------
    int inspect_open_ds_dapl() {
        auto path = filename("inspect_open");
        remove_file(path);
        {
            hid_t fd  = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            hsize_t d[1] = {3};
            hid_t sp  = H5Screate_simple(1, d, nullptr);
            hid_t ds  = H5Dcreate2(fd, "vec", H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dclose(ds); H5Sclose(sp); H5Fclose(fd);
        }
        h5::fd_t fd = h5::open(path, H5F_ACC_RDWR);
        probe("inspect_open_ds_dapl: before h5::open(fd, \"vec\")");
        h5::ds_t ds = h5::open(fd, "vec");
        probe("inspect_open_ds_dapl: after  h5::open");
        dump("ds.handle", static_cast<hid_t>(ds));
        dump("ds.dapl",   ds.dapl);
        remove_file(path);
        return 0;
    }

    // --- 3: forced fix — overwrite ds.dapl=H5P_DEFAULT, then write ------------
    int write_with_normalized_dapl() {
        auto path = filename("normalized_dapl");
        remove_file(path);
        h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
        h5::ds_t ds = h5::create<double>(fd, "vec", h5::current_dims{3});
        std::fprintf(stderr, "[WINDAPL]   pre-fix ds.dapl=%lld valid=%d\n",
            (long long)ds.dapl, (int)H5Iis_valid(ds.dapl));
        ds.dapl = H5P_DEFAULT;
        std::vector<double> data{1.0, 2.0, 3.0};
        probe("write_with_normalized_dapl: before h5::write(ds, data)");
        h5::write(ds, data);
        probe("write_with_normalized_dapl: after  h5::write — SURVIVED");
        remove_file(path);
        return 0;
    }

    // --- 4: write after h5::open (write path with KNOWN-good ds.dapl) ---------
    int write_after_open() {
        auto path = filename("write_after_open");
        remove_file(path);
        {
            hid_t fd  = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            hsize_t d[1] = {3};
            hid_t sp  = H5Screate_simple(1, d, nullptr);
            hid_t ds  = H5Dcreate2(fd, "vec", H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dclose(ds); H5Sclose(sp); H5Fclose(fd);
        }
        h5::fd_t fd = h5::open(path, H5F_ACC_RDWR);
        h5::ds_t ds = h5::open(fd, "vec");
        std::fprintf(stderr, "[WINDAPL]   open ds.dapl=%lld valid=%d\n",
            (long long)ds.dapl, (int)H5Iis_valid(ds.dapl));
        std::vector<double> data{1.0, 2.0, 3.0};
        probe("write_after_open: before h5::write(ds, data)");
        h5::write(ds, data);
        probe("write_after_open: after  h5::write");
        remove_file(path);
        return 0;
    }

    // --- 5: read after h5::create (does ds.dapl break read too?) --------------
    int read_after_create() {
        auto path = filename("read_after_create");
        remove_file(path);
        h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
        h5::ds_t ds = h5::create<double>(fd, "vec", h5::current_dims{3});
        std::fprintf(stderr, "[WINDAPL]   ds.dapl=%lld valid=%d\n",
            (long long)ds.dapl, (int)H5Iis_valid(ds.dapl));
        // Populate via raw HDF5 so reads have data without touching h5::write.
        double zeros[3] = {0.0, 0.0, 0.0};
        H5Dwrite(static_cast<hid_t>(ds), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, zeros);
        std::vector<double> out(3);
        probe("read_after_create: before h5::read(ds, vector)");
        h5::read(ds, out);
        probe("read_after_create: after  h5::read");
        std::fprintf(stderr, "[WINDAPL]   read size=%zu\n", out.size());
        remove_file(path);
        return 0;
    }

    // --- 6: inspect global default_dapl / dapl singletons --------------------
    int default_dapl_singleton() {
        probe("default_dapl_singleton: inspect h5::default_dapl and h5::dapl");
        dump("h5::default_dapl", static_cast<hid_t>(h5::default_dapl));
        dump("h5::dapl",         static_cast<hid_t>(h5::dapl));
        int r = H5Pexist(static_cast<hid_t>(h5::default_dapl), H5CPP_DAPL_HIGH_THROUGHPUT);
        std::fprintf(stderr, "[WINDAPL]   H5Pexist(default_dapl,...) -> %d\n", r);
        return 0;
    }

    // --- 7: default-construct resolved_type_t<double> ------------------------
    int inspect_resolved_type() {
        probe("inspect_resolved_type: default-ctor resolved_type_t<double>");
        h5::meta::resolved_type_t<double> t;
        dump("resolved_type_t<double>", static_cast<hid_t>(t));
        return 0;
    }

    // --- 8: H5Dget_access_plist via h5::get_dapl (CAPI round-trip) -----------
    int get_dapl_via_capi() {
        auto path = filename("get_dapl_capi");
        remove_file(path);
        h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
        h5::ds_t ds = h5::create<double>(fd, "vec", h5::current_dims{3});
        probe("get_dapl_via_capi: h5::get_dapl(ds) (H5Dget_access_plist)");
        h5::dapl_t dapl = h5::get_dapl(ds);
        dump("h5::get_dapl(ds)", static_cast<hid_t>(dapl));
        remove_file(path);
        return 0;
    }

    // --- 9: ds.dapl survival across std::move --------------------------------
    int ds_dapl_through_move() {
        auto path = filename("dapl_through_move");
        remove_file(path);
        h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
        h5::ds_t ds1 = h5::create<double>(fd, "vec", h5::current_dims{3});
        std::fprintf(stderr, "[WINDAPL]   pre-move  ds1.dapl=%lld\n", (long long)ds1.dapl);
        h5::ds_t ds2 = std::move(ds1);
        std::fprintf(stderr, "[WINDAPL]   post-move ds2.dapl=%lld\n", (long long)ds2.dapl);
        std::fprintf(stderr, "[WINDAPL]   post-move ds1.dapl=%lld\n", (long long)ds1.dapl);
        remove_file(path);
        return 0;
    }

    // --- 10: low-level write via INNER overload using forced default plists --
    // Mirrors H5Dwrite.hpp:21–28 with explicit dxpl=H5P_DEFAULT — should pass.
    int inner_write_with_default_dxpl() {
        auto path = filename("inner_write_default_dxpl");
        remove_file(path);
        h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
        h5::ds_t ds = h5::create<double>(fd, "vec", h5::current_dims{3});
        h5::sp_t mem  = h5::create_simple(static_cast<hsize_t>(3));
        h5::sp_t file{H5Dget_space(static_cast<hid_t>(ds))};
        h5::dxpl_t dxpl{H5P_DEFAULT};
        h5::select_all(mem);
        h5::select_all(file);
        double data[3] = {1.0, 2.0, 3.0};
        const double* ptr = data;
        // Pin the exact innermost overload (H5Dwrite.hpp:21) via function-pointer cast
        // so MSVC can't ambiguate with the variadic forwarders.
        void (*inner)(const h5::ds_t&, const h5::sp_t&, const h5::sp_t&,
                      const h5::dxpl_t&, const double*) = h5::write<double>;
        probe("inner_write_with_default_dxpl: before inner h5::write<double>");
        inner(ds, mem, file, dxpl, ptr);
        probe("inner_write_with_default_dxpl: after");
        remove_file(path);
        return 0;
    }
}

int main() {
#if H5CPP_WIN_DAPL_PROBE_MODE == 1
    return inspect_create_ds_dapl();
#elif H5CPP_WIN_DAPL_PROBE_MODE == 2
    return inspect_open_ds_dapl();
#elif H5CPP_WIN_DAPL_PROBE_MODE == 3
    return write_with_normalized_dapl();
#elif H5CPP_WIN_DAPL_PROBE_MODE == 4
    return write_after_open();
#elif H5CPP_WIN_DAPL_PROBE_MODE == 5
    return read_after_create();
#elif H5CPP_WIN_DAPL_PROBE_MODE == 6
    return default_dapl_singleton();
#elif H5CPP_WIN_DAPL_PROBE_MODE == 7
    return inspect_resolved_type();
#elif H5CPP_WIN_DAPL_PROBE_MODE == 8
    return get_dapl_via_capi();
#elif H5CPP_WIN_DAPL_PROBE_MODE == 9
    return ds_dapl_through_move();
#elif H5CPP_WIN_DAPL_PROBE_MODE == 10
    return inner_write_with_default_dxpl();
#else
    probe("unknown H5CPP_WIN_DAPL_PROBE_MODE");
    return 2;
#endif
}
