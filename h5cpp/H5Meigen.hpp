/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once
#include <hdf5.h>
#include "H5Tmeta.hpp"
#include <tuple>
#include <type_traits>
#include <array>

#if defined(EIGEN_CORE_H) || defined(H5CPP_USE_EIGEN3)
/*
	Matrix<typename Scalar,
       int RowsAtCompileTime,
       int ColsAtCompileTime,
       int Options = 0,          // ColMajor | RowMajor
       int MaxRowsAtCompileTime = RowsAtCompileTime,
       int MaxColsAtCompileTime = ColsAtCompileTime>
*/

namespace h5::meta {

    template<class T,int R,int C, int O> struct is_contiguous<::Eigen::Matrix<T,R,C,O>> : std::true_type {};
    template<class T,int R,int C, int O> struct is_contiguous<::Eigen::Array<T,R,C,O>> : std::true_type {};
}

namespace h5::impl {
	// 1.) object -> H5T_xxx
	// Register Eigen types in has_explicit_decay to prevent ambiguity with the
	// structural decay fallback (Eigen matrices expose value_type = Scalar).
	template<class T,int R,int C, int O>
	struct detail::has_explicit_decay<::Eigen::Matrix<T,R,C,O>> : std::true_type {};
	template<class T,int R,int C, int O>
	struct detail::has_explicit_decay<::Eigen::Array<T,R,C,O>>  : std::true_type {};

	// 1.) object -> H5T_xxx
	template<class T,int R,int C, int O> struct decay<::Eigen::Matrix<T,R,C,O>>{ using type = T; };
	template<class T,int R,int C, int O> struct decay<::Eigen::Array<T,R,C,O>>{ using type = T; };
	    

	// get read access to datastaore
	template<class T,int R,int C,int O, int MR=R,int MC=C>
	const T* data(const ::Eigen::Matrix<T,R,C,O,MR,MC>& ref ){
			return ref.data();
	}
	// read write access
	template<class T,int R,int C,int O, int MR=R,int MC=C>
	T* data(::Eigen::Matrix<T,R,C,O,MR,MC>& ref ){
			return ref.data();
	}
	// get read access to datastaore
	template<class T,int R,int C,int O, int MR=R,int MC=C>
	const T* data(const ::Eigen::Array<T,R,C,O,MR,MC>& ref ){
			return ref.data();
	}
	// read write access
	template<class T,int R,int C,int O, int MR=R,int MC=C>
	T* data(::Eigen::Array<T,R,C,O,MR,MC>& ref ){
			return ref.data();
	}
	// determine rank and dimensions
	// MATRICES
	template<class T,int R,int C,int MR=R,int MC=C>
	inline std::array<size_t,2> size( const ::Eigen::Matrix<T,R,C,::Eigen::RowMajor,MR,MC>& ref ){
		return {(hsize_t)ref.rows(),(hsize_t)ref.cols()};
	}
	template<class T,int R,int C,int MR=R,int MC=C>
	inline std::array<size_t,2> size( const ::Eigen::Matrix<T,R,C,::Eigen::ColMajor,MR,MC>& ref ){
		return {(hsize_t)ref.cols(), (hsize_t)ref.rows()};
	}
	// ARRAYS
	template<class T,int R,int C,int MR=R,int MC=C>
	inline std::array<size_t,2> size( const ::Eigen::Array<T,R,C,::Eigen::RowMajor,MR,MC>& ref ){
		return {(hsize_t)ref.rows(),(hsize_t)ref.cols()};
	}

	template<class T,int R,int C,int MR=R,int MC=C>
	inline std::array<size_t,2> size( const ::Eigen::Array<T,R,C,::Eigen::ColMajor,MR,MC>& ref ){
		return {(hsize_t)ref.cols(), (hsize_t)ref.rows()};
	}

	// rank
	template<class T,int R,int C,int O,int MR,int MC>
	struct rank<::Eigen::Matrix<T,R,C,O,MR,MC>> : public std::integral_constant<size_t,2>{};
	template<class T,int R,int C,int O,int MR,int MC>
	struct rank<::Eigen::Array<T,R,C,O,MR,MC>>  : public std::integral_constant<size_t,2>{};
	// CTOR-s
	// MATRICES
	template<class T,int R,int C>
	struct get<::Eigen::Matrix<T,R,C,::Eigen::RowMajor>> {
		static inline ::Eigen::Matrix<T,R,C,::Eigen::RowMajor> ctor( std::array<size_t,2> dims ){
			return ::Eigen::Matrix<T,R,C,::Eigen::RowMajor>( dims[0], dims[1] );
	}};
	template<class T,int R,int C, int MR, int MC>
	struct get<::Eigen::Matrix<T,R,C,::Eigen::RowMajor,MR,MC>> {
		static inline ::Eigen::Matrix<T,R,C,::Eigen::RowMajor,MR,MC> ctor( std::array<size_t,2> dims ){
			return ::Eigen::Matrix<T,R,C,::Eigen::RowMajor,MR,MC>( dims[0], dims[1] );
	}};
	template<class T,int R,int C>
	struct get<::Eigen::Matrix<T,R,C,::Eigen::ColMajor>> {
		static inline ::Eigen::Matrix<T,R,C,::Eigen::ColMajor> ctor( std::array<size_t,2> dims ){
			return ::Eigen::Matrix<T,R,C,::Eigen::ColMajor>( dims[1], dims[0] );
	}};
	template<class T,int R,int C, int MR, int MC>
	struct get<::Eigen::Matrix<T,R,C,::Eigen::ColMajor,MR,MC>> {
		static inline ::Eigen::Matrix<T,R,C,::Eigen::ColMajor,MR,MC> ctor( std::array<size_t,2> dims ){
			return ::Eigen::Matrix<T,R,C,::Eigen::ColMajor,MR,MC>( dims[1], dims[0] );
	}};
	// ARRAYS
	template<class T,int R,int C>
	struct get<::Eigen::Array<T,R,C,::Eigen::RowMajor>> {
		static inline ::Eigen::Array<T,R,C,::Eigen::RowMajor> ctor( std::array<size_t,2> dims ){
			return ::Eigen::Array<T,R,C,::Eigen::RowMajor>( dims[0], dims[1] );
	}};
	template<class T,int R,int C, int MR, int MC>
	struct get<::Eigen::Array<T,R,C,::Eigen::RowMajor,MR,MC>> {
		static inline ::Eigen::Array<T,R,C,::Eigen::RowMajor,MR,MC> ctor( std::array<size_t,2> dims ){
			return ::Eigen::Array<T,R,C,::Eigen::RowMajor,MR,MC>( dims[0], dims[1] );
	}};
	template<class T,int R,int C>
	struct get<::Eigen::Array<T,R,C,::Eigen::ColMajor>> {
		static inline ::Eigen::Array<T,R,C,::Eigen::ColMajor> ctor( std::array<size_t,2> dims ){
			return ::Eigen::Array<T,R,C,::Eigen::ColMajor>( dims[1], dims[0] );
	}};
	template<class T,int R,int C, int MR, int MC>
	struct get<::Eigen::Array<T,R,C,::Eigen::ColMajor,MR,MC>> {
		static inline ::Eigen::Array<T,R,C,::Eigen::ColMajor,MR,MC> ctor( std::array<size_t,2> dims ){
			return ::Eigen::Array<T,R,C,::Eigen::ColMajor,MR,MC>( dims[1], dims[0] );
	}};
}

namespace h5::meta {
    template<class T,int R,int C,int O, int MR=R,int MC=C>
    T* data(const ::Eigen::Matrix<T,R,C,O,MR,MC>& ref ){
            return const_cast<T*>( ref.data() );
    }
    template<class T,int R,int C,int O, int MR=R,int MC=C>
    T* data(const ::Eigen::Array<T,R,C,O,MR,MC>& ref ){
            return const_cast<T*>( ref.data() );
    }
    template<class T,int R,int C,int MR=R,int MC=C>
    inline std::array<size_t,2> size( const ::Eigen::Matrix<T,R,C,::Eigen::RowMajor,MR,MC>& ref ){
        return {(hsize_t)ref.rows(),(hsize_t)ref.cols()};
    }
    template<class T,int R,int C,int MR=R,int MC=C>
    inline std::array<size_t,2> size( const ::Eigen::Matrix<T,R,C,::Eigen::ColMajor,MR,MC>& ref ){
        return {(hsize_t)ref.cols(), (hsize_t)ref.rows()};
    }
    template<class T,int R,int C,int MR=R,int MC=C>
    inline std::array<size_t,2> size( const ::Eigen::Array<T,R,C,::Eigen::RowMajor,MR,MC>& ref ){
        return {(hsize_t)ref.rows(),(hsize_t)ref.cols()};
    }
    template<class T,int R,int C,int MR=R,int MC=C>
    inline std::array<size_t,2> size( const ::Eigen::Array<T,R,C,::Eigen::ColMajor,MR,MC>& ref ){
        return {(hsize_t)ref.cols(), (hsize_t)ref.rows()};
    }
}

namespace h5::meta {
    template<class T,int R,int C, int O> struct decay<::Eigen::Matrix<T,R,C,O>> : h5::impl::decay<::Eigen::Matrix<T,R,C,O>> {};
    template<class T,int R,int C, int O> struct decay<::Eigen::Array<T,R,C,O>> : h5::impl::decay<::Eigen::Array<T,R,C,O>> {};
    template<class T,int R,int C,int O,int MR,int MC> struct rank<::Eigen::Matrix<T,R,C,O,MR,MC>> : h5::impl::rank<::Eigen::Matrix<T,R,C,O,MR,MC>> {};
    template<class T,int R,int C,int O,int MR,int MC> struct rank<::Eigen::Array<T,R,C,O,MR,MC>> : h5::impl::rank<::Eigen::Array<T,R,C,O,MR,MC>> {};
    template<class T,int R,int C> struct get<::Eigen::Matrix<T,R,C,::Eigen::RowMajor>> : h5::impl::get<::Eigen::Matrix<T,R,C,::Eigen::RowMajor>> {};
    template<class T,int R,int C, int MR, int MC> struct get<::Eigen::Matrix<T,R,C,::Eigen::RowMajor,MR,MC>> : h5::impl::get<::Eigen::Matrix<T,R,C,::Eigen::RowMajor,MR,MC>> {};
    template<class T,int R,int C> struct get<::Eigen::Matrix<T,R,C,::Eigen::ColMajor>> : h5::impl::get<::Eigen::Matrix<T,R,C,::Eigen::ColMajor>> {};
    template<class T,int R,int C, int MR, int MC> struct get<::Eigen::Matrix<T,R,C,::Eigen::ColMajor,MR,MC>> : h5::impl::get<::Eigen::Matrix<T,R,C,::Eigen::ColMajor,MR,MC>> {};
    template<class T,int R,int C> struct get<::Eigen::Array<T,R,C,::Eigen::RowMajor>> : h5::impl::get<::Eigen::Array<T,R,C,::Eigen::RowMajor>> {};
    template<class T,int R,int C, int MR, int MC> struct get<::Eigen::Array<T,R,C,::Eigen::RowMajor,MR,MC>> : h5::impl::get<::Eigen::Array<T,R,C,::Eigen::RowMajor,MR,MC>> {};
    template<class T,int R,int C> struct get<::Eigen::Array<T,R,C,::Eigen::ColMajor>> : h5::impl::get<::Eigen::Array<T,R,C,::Eigen::ColMajor>> {};
    template<class T,int R,int C, int MR, int MC> struct get<::Eigen::Array<T,R,C,::Eigen::ColMajor,MR,MC>> : h5::impl::get<::Eigen::Array<T,R,C,::Eigen::ColMajor,MR,MC>> {};
}

#endif
