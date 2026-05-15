/* Copyright (c) 2018-2026 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargalabs.com> */

/**
 * @file s3.cpp
 * @brief Open HDF5 files on S3 using the Read-Only S3 VFD (ROS3).
 *
 * Requires HDF5 built with ROS3 support (H5_HAVE_ROS3_VFD).
 * Build: cmake -DH5CPP_BUILD_EXAMPLES=ON .
 * Run:   AWS_REGION=us-east-1 ./s3
 */

#include <h5cpp/all>
#include <iostream>
#include <stdexcept>

#ifndef H5_HAVE_ROS3_VFD
int main() {
    std::cerr << "This example requires HDF5 built with ROS3 VFD support.\n"
              << "Rebuild HDF5 with --enable-ros3-vfd (and libcurl + OpenSSL).\n";
    return 1;
}
#else
#include <H5FDros3.h>

// ---------------------------------------------------------------------------
// Helper: print top-level dataset names in a file
// ---------------------------------------------------------------------------
static void list_datasets(h5::fd_t& fd, const std::string& label) {
    std::cout << label << "\n";
    // Iterate root group — using low-level CAPI for brevity
    hsize_t idx = 0;
    H5Literate(static_cast<hid_t>(fd), H5_INDEX_NAME, H5_ITER_NATIVE, &idx,
        [](hid_t, const char* name, const H5L_info_t*, void*) -> herr_t {
            std::cout << "  /" << name << "\n";
            return 0;
        }, nullptr);
}

int main() {
    try {
        // ── 1. Unauthenticated: public bucket ───────────────────────────────
        {
            auto fd = h5::open("s3://rhdf5-public/h5ex_t_array.h5",
                               H5F_ACC_RDONLY, h5::ros3{false, "eu-central-1", "", ""});
            list_datasets(fd, "Public bucket (no auth):");
        }

        // ── 2. Authenticated: long-term credentials ─────────────────────────
        // Read credentials from environment — never hard-code.
        const char* region = std::getenv("AWS_REGION");
        const char* key_id = std::getenv("AWS_ACCESS_KEY_ID");
        const char* secret = std::getenv("AWS_SECRET_ACCESS_KEY");

        if (region && key_id && secret) {
            auto fd = h5::open("s3://my-private-bucket/archive/run42.h5",
                H5F_ACC_RDONLY, h5::ros3{true, region, key_id, secret});
            list_datasets(fd, "Private bucket (long-term creds):");

            // Normal h5cpp I/O — no special S3 code after open
            auto temperatures = h5::read<std::vector<double>>(fd, "/sensor/temperature");
            std::cout << "  read " << temperatures.size() << " temperature values\n";
        } else {
            std::cout << "Set AWS_REGION, AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY "
                         "to test authenticated access.\n";
        }

#if H5_VERSION_GE(1,12,1)
        // ── 3. Temporary credentials: STS / IAM role ────────────────────────
        const char* token = std::getenv("AWS_SESSION_TOKEN");
        if (region && key_id && secret && token) {
            auto fd = h5::open("s3://my-secure-bucket/data.h5",
                H5F_ACC_RDONLY, h5::ros3{true, region, key_id, secret, token});
            list_datasets(fd, "Secure bucket (STS session token):");
        }
#endif

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
