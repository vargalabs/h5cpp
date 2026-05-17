#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/H5Pall.hpp>

TEST_CASE("instantiate commonly used property types") {
    h5::chunk c{10};
    h5::gzip g{6};
    h5::deflate d{6};
    h5::flag::fletcher32 f;
    h5::flag::shuffle s;
    h5::flag::nbit n;
    h5::fill_time ft{H5D_FILL_TIME_ALLOC};
    h5::alloc_time at{H5D_ALLOC_TIME_EARLY};
    h5::layout l{H5D_CHUNKED};
    h5::libver_bounds lv({H5F_LIBVER_LATEST, H5F_LIBVER_LATEST});
    h5::fclose_degree fd_{H5F_CLOSE_WEAK};
    h5::char_encoding ce{H5T_CSET_UTF8};
    h5::create_intermediate_group cig{1};
    CHECK(true);
}

TEST_CASE("instantiate dead fcpl property types") {
    h5::sizes sz;
    h5::sym_k sk;
    h5::istore_k ik;
    h5::shared_mesg_nindexes smn;
    h5::shared_mesg_index smi;
    h5::shared_mesg_phase_change smpc;
    h5::userblock ub;
    CHECK(true);
}

TEST_CASE("instantiate dead fapl property types") {
    h5::cache ca;
    h5::alignment al;
    h5::meta_block_size mbs;
    h5::sieve_buf_size sbs;
    h5::elink_file_cache_size efcs;
    CHECK(true);
}

TEST_CASE("instantiate dead gcpl property types") {
    h5::local_heap_size_hint lhsh;
    h5::link_creation_order lco;
    h5::est_link_info eli;
    h5::link_phase_change lpc;
    CHECK(true);
}

TEST_CASE("property types can be daisy-chained with operator|") {
    h5::dcpl_t props = h5::chunk{10} | h5::gzip{6} | h5::flag::fletcher32{};
    CHECK(static_cast<hid_t>(props) > 0);
}

// regression guard for issue #239 — combining a property builder with an
// existing property-list handle (v1.10-era idiom) must remain supported.
TEST_CASE("property builder can be chained against an existing dcpl handle") {
    h5::dcpl_t base = h5::gzip{6};
    REQUIRE(static_cast<hid_t>(base) > 0);

    h5::dcpl_t combined = h5::chunk{64 * 1024} | base;
    REQUIRE(static_cast<hid_t>(combined) > 0);
    CHECK(static_cast<hid_t>(combined) != static_cast<hid_t>(base)); // deep copy

    // chunk dims set by the LHS builder survive the merge
    hsize_t chunk_dims[H5S_MAX_RANK] = {0};
    int rank = H5Pget_chunk(static_cast<hid_t>(combined), H5S_MAX_RANK, chunk_dims);
    CHECK(rank == 1);
    CHECK(chunk_dims[0] == 64 * 1024);

    // deflate filter inherited from base — scan all filters since pipeline order
    // depends on H5Pcopy + lazy-apply semantics.
    int nfilters = H5Pget_nfilters(static_cast<hid_t>(combined));
    CHECK(nfilters == 1);
    bool found_deflate = false;
    for (int i = 0; i < nfilters; ++i) {
        unsigned flags = 0, filter_config = 0;
        size_t cd_nelmts = 1;
        unsigned cd_values[1] = {0};
        char name[16] = {0};
        H5Z_filter_t f = H5Pget_filter2(static_cast<hid_t>(combined), i, &flags,
            &cd_nelmts, cd_values, sizeof(name), name, &filter_config);
        if (f == H5Z_FILTER_DEFLATE && cd_values[0] == 6)
            found_deflate = true;
    }
    CHECK(found_deflate);
}
