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

int main() {
    try {
        // ── 1. Unauthenticated: public bucket ───────────────────────────────
        {
            auto fd      = h5::open("s3://rhdf5-public/h5ex_t_array.h5",
                                    H5F_ACC_RDONLY, h5::ros3{false, "eu-central-1", "", ""});
            auto entries = h5::ls(fd, "/");
            std::cout << "Public bucket (no auth):\n";
            for (auto& name : entries) std::cout << "  /" << name << "\n";
        }

        // ── 2. Authenticated: long-term credentials ─────────────────────────
        // Read credentials from environment — never hard-code.
        const char* region = std::getenv("AWS_REGION");
        const char* key_id = std::getenv("AWS_ACCESS_KEY_ID");
        const char* secret = std::getenv("AWS_SECRET_ACCESS_KEY");

        if (region && key_id && secret) {
            auto fd      = h5::open("s3://my-private-bucket/archive/run42.h5",
                                    H5F_ACC_RDONLY, h5::ros3{true, region, key_id, secret});
            auto entries = h5::ls(fd, "/");
            std::cout << "Private bucket (long-term creds):\n";
            for (auto& name : entries) std::cout << "  /" << name << "\n";

            // Normal h5cpp I/O — no S3-specific code needed after open
            auto temps = h5::read<std::vector<double>>(fd, "/sensor/temperature");
            std::cout << "  read " << temps.size() << " temperature values\n";
        } else {
            std::cout << "Set AWS_REGION, AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY "
                         "to test authenticated access.\n";
        }

#if H5_VERSION_GE(1,12,1)
        // ── 3. Temporary credentials: STS / IAM role ────────────────────────
        const char* token = std::getenv("AWS_SESSION_TOKEN");
        if (region && key_id && secret && token) {
            auto fd      = h5::open("s3://my-secure-bucket/data.h5",
                                    H5F_ACC_RDONLY, h5::ros3{true, region, key_id, secret, token});
            auto entries = h5::ls(fd, "/");
            std::cout << "Secure bucket (STS session token):\n";
            for (auto& name : entries) std::cout << "  /" << name << "\n";
        }
#endif

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
