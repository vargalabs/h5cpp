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
