/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */

#pragma once
#include <string>
#include <vector>
#include <array>
#include <initializer_list>
#include <type_traits>

namespace h5::impl {
/*STL: */
	// 2.) filter is_xxx_type
	// 4.) write access
	// 5.) obtain dimensions of extents
	// 6.) ctor with right dimensions

	// 1.) object -> H5T_xxx
	template <class T, class...> struct decay{ using type = T; };

	template <class T> struct decay<const T>{ using type = T; };
	template <class T> struct decay<const T*>{ using type = T*; };
	template <class T> struct decay<std::basic_string<T>>{ using type = const T*; };
	template <class T, signed N> struct decay<const T[N]>{ using type = T*; };
	template <class T, signed N> struct decay<T[N]>{ using type = T*; };

	template <class T> struct decay<std::initializer_list<const T*>>{ using type = const T*; };
	template <class T> struct decay<std::initializer_list<T*>>{ using type = T*; };
	template <class T> struct decay<std::initializer_list<T>>{ using type = T; };

	template <class T> struct decay<std::vector<const T*>>{ using type = const T*; };
	template <class T> struct decay<std::vector<T*>>{ using type = T*; };
	template <class T> struct decay<std::vector<T>>{ using type = T; };

	// helpers
	template <class T>
		using is_scalar = std::integral_constant<bool,
			std::is_integral<T>::value || (std::is_standard_layout_v<T> && std::is_trivial_v<T>) || std::is_same<T,std::string>::value>;
	template <class T, class B = typename impl::decay<T>::type>
		using is_rank01 = std::integral_constant<bool,
			std::is_same<T,std::initializer_list<B>>::value || 
			std::is_same<T,std::vector<B>>::value >;

	template<class T> struct rank : public std::integral_constant<size_t,0>{};
	template<> struct rank<std::string>: public std::integral_constant<size_t,1>{};
	template<class T> struct rank<T*>: public std::integral_constant<size_t,1>{};
	template<class T> struct rank<std::vector<T>>: public std::integral_constant<size_t,1>{};

	// 3.) read access
	template <class T> inline typename std::enable_if<std::is_integral<T>::value || (std::is_standard_layout_v<T> && std::is_trivial_v<T>),
		const T*>::type data( const T& ref ){ return &ref; }
	template<class T> inline typename std::enable_if< impl::is_scalar<T>::value,
		const T*>::type data( const std::initializer_list<T>& ref ){ return ref.begin(); }
	inline const char* const* data( const std::initializer_list<const char*>& ref ){ return ref.begin(); }
	inline const char* data( const std::string& ref ){ return ref.c_str(); }
	template <class T, class A> inline const T* data( const std::vector<T, A>& ref ){
		return ref.data();
	}
	template <class T, class A> inline T* data( std::vector<T, A>& ref ){
		return ref.data();
	}
	// 4.) write access
	template <class T> inline std::enable_if_t<std::is_integral_v<T>,
	T*> data( T& ref ){ return &ref; }
	//template <class T> inline std::enable_if_t<std::is_integral_v<T>,
	//	const T*> data( const T& ref ){ return &ref; }
		// 5.) obtain dimensions of extents
		template <class T> inline constexpr std::enable_if_t< impl::is_scalar<T>::value,
			std::array<size_t,0>> size( T ){ return{}; }
	template <class T, class A> inline std::array<size_t,1> size( const std::vector<T, A>& ref ){ return {ref.size()}; }
	template <class T> inline std::array<size_t,1> size( const std::initializer_list<T>& ref ){ return {ref.size()}; }
	template <class T> inline constexpr std::enable_if_t<!impl::is_scalar<T>::value,
		std::array<size_t,0>> size( const T& ){ return{}; }
		// 6.) ctor with right dimensions
		template <class T> struct get {
		static inline T ctor( std::array<size_t,impl::rank<T>::value> ){
				return T(); }};
	template<class T>
	struct get<std::vector<T>> {
		static inline std::vector<T> ctor( std::array<size_t,1> dims ){
			return std::vector<T>( dims[0] );
	}};
}
