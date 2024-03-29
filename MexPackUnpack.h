#ifndef __MexMatrix
#define __MexMatrix

#include "boost/pfr.hpp"
#include "mex.h"
#include <Eigen>
#include <MexPackUnpackTypes.h>
#include <map>
#include <string>
#include <tuple>
#include <typeindex>
#include <utility>
#include <variant>
#include <algorithm>
#include <type_traits>

#ifdef USE_MDSPAN
#include <experimental/mdspan>

namespace stdex = std::experimental;
//using stdex = std::experimental;

/*
template <class T>
using span_1d_dynamic_left = stdex::mdspan<T, stdex::extents<stdex::dynamic_extent>, stdex::layout_left>;

template <class T>
using span_2d_dynamic_left = stdex::mdspan<T, stdex::extents<stdex::dynamic_extent, stdex::dynamic_extent>, stdex::layout_left>;

template <class T>
using span_3d_dynamic_left = stdex::mdspan<T, stdex::extents<stdex::dynamic_extent, stdex::dynamic_extent,stdex::dynamic_extent>, stdex::layout_left>;

template <class T>
using span_4d_dynamic_left = stdex::mdspan<T, stdex::extents<stdex::dynamic_extent, stdex::dynamic_extent,stdex::dynamic_extent,stdex::dynamic_extent>, stdex::layout_left>;
*/

#endif

// COMPLEX_SPLIT determines whether complex data viewed as interleaved or split
// since matlab 2018a+  (when compiled with -R2018a) api for complex data requires mex functions not present in the octave or earlier matlab verisons mex.h. We have #ifdefs in parts of code for
// complex data to let the code work with earlier or later versions (this could probably be done more cleanly). 

#define COMPLEX_SPLIT
#ifdef MX_HAS_INTERLEAVED_COMPLEX // Should be set if using matlab 2018a or later
#if MX_HAS_INTERLEAVED_COMPLEX
#undef COMPLEX_SPLIT
#endif
#endif

// Helper Metafunction to extract Nth element of type pack

template <int i, typename... Args> struct NthElement;

template <int i, typename As, typename... Args> struct NthElement<i, As, Args...> { using type = typename NthElement<i - 1, Args...>::type; };

template <typename As, typename... Args> struct NthElement<0, As, Args...> { using type = As; };

using namespace MexPackUnpackTypes;

// adapted  from https://stackoverflow.com/questions/34111060/c-check-if-the-template-type-is-one-of-the-variadic-template-types
// checks if first type in type list matches any remaining ones
template <typename Kind, typename... Kinds> constexpr bool in_type_list() {
  /* The following expands to :
   * std::is_same_v<Kind, Kind0> || std::is_same_v<Kind, Kind1> || ... */
  if constexpr ((std::is_same_v<Kind, Kinds> || ...)) {
    return true;
  } else {
    return false;
  }
};

// tuple to struct utility from https://stackoverflow.com/questions/38561242/struct-to-from-stdtuple-conversion
namespace details {

template <typename result_type, typename... types, std::size_t... indices>
result_type make_struct(std::tuple<types...> t, std::index_sequence<indices...>) // &, &&, const && etc.
{
  return {std::get<indices>(t)...};
}

} // namespace details

template <typename result_type, typename... types>
result_type make_struct(std::tuple<types...> t) // &, &&, const && etc.
{
  return details::make_struct<result_type, types...>(t, std::index_sequence_for<types...>{}); // if there is repeated types, then the change for using std::index_sequence_for is trivial
}

//Detect tuple (or pair)
template <typename> struct is_tuple : std::false_type {};
template <typename... T> struct is_tuple<std::tuple<T...>> : std::true_type {};
template <typename... T> struct is_tuple<std::pair<T...>> : std::true_type {};

//Get type of Nth entry in a tuple
template <int i, typename...> struct tuple_type_N;

template <int i, typename...T>  struct tuple_type_N<i,std::tuple<T...> >{
  using type = typename NthElement<i,T...>::type;
};
template <int i, typename...T>  struct tuple_type_N<i,std::pair<T...> >{
  using type = typename NthElement<i,T...>::type;
};

template < typename >
struct is_complex : std::false_type {};

template < typename T >
struct is_complex<std::complex<T>> : std::true_type {
  using real_type=T;
};


//https://stackoverflow.com/questions/31762958/check-if-class-is-a-template-specialization
template <class T, template <class...> class Template>
struct is_specialization : std::false_type {};

template <template <class...> class Template, class... Args>
struct is_specialization<Template<Args...>, Template> : std::true_type {};

//https://stackoverflow.com/questions/34672441/stdis-base-of-for-template-classes
template < template <typename...> class base,typename derived>
struct is_base_of_template_impl
{
    template<typename... Ts>
    static constexpr std::true_type  test(const base<Ts...> *);
    static constexpr std::false_type test(...);
    using type = decltype(test(std::declval<derived*>()));
};

template < template <typename...> class base,typename derived>
using is_base_of_template = typename is_base_of_template_impl<base,derived>::type;

// Class to extract Matlab/Octave arrays passed to a mex function

/*
Usage:
  MexUnpacker<double,int,EDCSM,EFRM,RDP> my_unpack(nrhs,prhs);
  auto [a,b,c_split,e,g] = my_unpack.unpackMex();

  Note that in general one must be careful using auto with eigen objects due to its use of
  expression templates. However, in this case auto correctly deduces returned eigen objects types.
*/
template <typename... T> class MexUnpacker {
public:
  int nrhs;
  const mxArray **prhs;

  MexUnpacker(int nrhs_, const mxArray *prhs_[]) : nrhs{nrhs_}, prhs{prhs_} {}
  MexUnpacker(int nrhs_, mxArray **prhs_) : nrhs{nrhs_}, prhs{prhs_} {}

  // Note for all the get_ helper functions, the second argument is a pointer of the type we are returning
  // this pointer is initialized to null and ignored and solely present to allow us to pick which overload to use.
  // There is probably a more elegant way to do this.

  // Deal with all scalar cases that are real (no imaginary component)
  template <typename S> std::enable_if_t<std::is_scalar<S>::value, S> get_(int i, S *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (std::any_of(dims,dims+ndims,[](auto x) {return x!=1;}))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a scalar\n");

    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string("\n");
    }
    return static_cast<S>(mxGetScalar(prhs[i]));
  }

  // Deal with all pointer tuples that are real (no imaginary component)
  template <typename S> std::enable_if_t<std::is_scalar<S>::value, ptr_tuple<S>> get_(int i, ptr_tuple<S> *ignored) {
    //   get_(int i, RFP *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 2)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 2 dimensional\n");
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");
    }

    return {reinterpret_cast<S *>(mxGetData(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i])};
  }
  
  //Cases where we extract to real 3D array to tuple of ptr plus dims
  template <typename S> std::enable_if_t<std::is_scalar<S>::value, ptr_tuple_3dim<S>> get_(int i, ptr_tuple_3dim<S> *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 3)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 3 dimensional\n");

    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");
    }
    return {reinterpret_cast<S *>(mxGetData(prhs[i])), dims[0], dims[1], dims[2]};
  }

  //Handle strings (actually character arrays)
  std::string get_(int i, std::string *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxCHAR_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not an string\n");
    char *c_string = mxArrayToString(prhs[i]);
    return std::string(c_string);
  }

  //Cases where we extract to real Eigen Map
  template <class S,typename U = typename S::Scalar>
  std::enable_if_t<isEigenMap<S>::value&& (Eigen::NumTraits<U>::IsComplex==0),S > get_(int i, S *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 2)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 2 dimensional\n");
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<U>::mxClass)){
      std::string identifier{mxClassTraits<U>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");      
    }    
    S new_map(reinterpret_cast<U *>(mxGetData(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    return new_map;
  }

#ifdef USE_MDSPAN

#ifdef COMPLEX_SPLIT


 template <typename S,size_t... U,typename W> 
 std::enable_if_t<std::is_scalar<S>::value, stdex::mdspan<S, stdex::extents<U...>, W> >
  get_(int i,  stdex::mdspan<S, stdex::extents<U...>, W>* ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (ndims != stdex::extents<U...>::rank() && !(stdex::extents<U...>::rank()==1&&ndims==2 && dims[1]==1) )
      throw std::string("Argument ") + std::to_string(i) + std::string(" not " +std::to_string(stdex::extents<U...>::rank()) +" dimensional\n");
    
    std::array<mwSize,stdex::extents<U...>::rank()> dims_array;
    std::copy(dims,dims+ndims,dims_array.begin());

    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");
    } 
    return stdex::mdspan<S, stdex::extents<U...>,W>{reinterpret_cast<S *>(mxGetData(prhs[i])), dims_array};
  }

 template <typename S,size_t... U,typename W> 
 std::pair< stdex::mdspan<S, stdex::extents<U...>, W> ,stdex::mdspan<S, stdex::extents<U...>, W>  >
  get_(int i,  std::pair< stdex::mdspan<S, stdex::extents<U...>, W> ,stdex::mdspan<S, stdex::extents<U...>, W>  >
* ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (ndims != stdex::extents<U...>::rank() && !(stdex::extents<U...>::rank()==1&&ndims==2 && dims[1]==1))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not " +std::to_string(stdex::extents<U...>::rank()) +" dimensional\n");

    std::array<mwSize,stdex::extents<U...>::rank()> dims_array;
    std::copy(dims,dims+ndims,dims_array.begin());

    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");
    } 
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    return {stdex::mdspan<S, stdex::extents<U...>, W>{reinterpret_cast<S *>(mxGetPr(prhs[i])), dims_array},
      stdex::mdspan<S, stdex::extents<U...>, W>{reinterpret_cast<S *>(mxGetPi(prhs[i])), dims_array}  };  
  }

#else
     

 template <typename S,size_t... U,typename W> 
 stdex::mdspan<S, stdex::extents<U...>, W>
  get_(int i,  stdex::mdspan<S, stdex::extents<U...>, W>* ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (ndims != stdex::extents<U...>::rank() && !(stdex::extents<U...>::rank()==1&&ndims==2 && dims[1]==1) )
      throw std::string("Argument ") + std::to_string(i) + std::string(" not " +std::to_string(stdex::extents<U...>::rank()) +" dimensional\n");

    std::array<mwSize,stdex::extents<U...>::rank()> dims_array;
    std::copy(dims,dims+ndims,dims_array.begin());

    if constexpr (!is_complex<S>::value){

      if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
        std::string identifier{mxClassTraits<S>::name};
        throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");
      } 
      return stdex::mdspan<S, stdex::extents<U...>,W>{reinterpret_cast<S *>(mxGetData(prhs[i])), dims_array};
    }else{
      if (!mxIsComplex(prhs[i])) {
        throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
      }
      if (!(mxGetClassID(prhs[i]) == mxClassTraits<typename is_complex<S>::real_type>::mxClass)) {
        std::string identifier{mxClassTraits<typename is_complex<S>::real_type>::name};
        throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");
      } 
      return stdex::mdspan<S, stdex::extents<U...>, W>{reinterpret_cast<S *>(mxGetData(prhs[i])), dims_array};
    }
  }


#endif

/*
 template <typename S> 
 std::enable_if_t<std::is_scalar<S>::value,span_2d_dynamic_left<S> >     
  get_(int i, span_2d_dynamic_left<S> *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (ndims != 2)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 2 dimensional\n");
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");
    }
    return span_2d_dynamic_left<S>{reinterpret_cast<S *>(mxGetData(prhs[i])), dims[0], dims[1]};
  }
  
 template <typename S> 
 std::enable_if_t<std::is_scalar<S>::value,span_3d_dynamic_left<S> >     
  get_(int i, span_3d_dynamic_left<S> *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (ndims != 3)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 2 dimensional\n");
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");
    }
    return span_3d_dynamic_left<S>{reinterpret_cast<S *>(mxGetData(prhs[i])), dims[0], dims[1],dims[2]};
  }

 template <typename S> 
 std::enable_if_t<std::is_scalar<S>::value,span_4d_dynamic_left<S> >     
  get_(int i, span_4d_dynamic_left<S> *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (ndims != 4)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 2 dimensional\n");
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");
    }
    return span_4d_dynamic_left<S>{reinterpret_cast<S *>(mxGetData(prhs[i])), dims[0], dims[1],dims[2],dims[3]};
  }
  */


#endif
  
  // Behavior of complex arrays depends on whether arrays are interleaved or split
  // Octave and MATLAB before 2018a use split real and complex
  // Matlab when compiled with -R2018a uses interleaved real and complex
  // We return Eigen Maps which wrap the underlying pointers from octave
  // If interleaved the Eigen Maps is complex, if split we return a pair of real eigen maps


#ifdef COMPLEX_SPLIT

  //Cases where we extract split complex matrix to a pair of real Eigen Maps
  template <class S,
	    typename W= std::enable_if_t<is_tuple<S>::value && (isEigenMap<typename tuple_type_N<0,S>::type>::value),typename tuple_type_N<0,S>::type>,
	    typename U=typename W::Scalar  >
      S get_(int i, S *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 2)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 2 dimensional\n");
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<U>::mxClass)){
      std::string identifier{mxClassTraits<U>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");      
    }    
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    W new_map_real(reinterpret_cast<U *>(mxGetPr(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    W new_map_imag(reinterpret_cast<U *>(mxGetPi(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    return {new_map_real, new_map_imag};
  }

  //Complex scalar split
  template<typename S> 
  std::complex<S> get_(int i, std::complex<S> *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (std::any_of(dims,dims+ndims,[](auto x) {return x!=1;}))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a scalar\n");

    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)){
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");      
    }    
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    return std::complex<S>(reinterpret_cast<S *>(mxGetPr(prhs[i]))[0], reinterpret_cast<S *>(mxGetPi(prhs[i]))[0]);
  }

  //Complex 3dim array split to pair of pointers and dims  
  template <typename S> ptr_tuple_3dim_CS<S> get_(int i, ptr_tuple_3dim_CS<S> *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 3)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 3 dimensional\n");
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    return {std::pair{reinterpret_cast<S *>(mxGetPr(prhs[i])), reinterpret_cast<S *>(mxGetPi(prhs[i]))}, dims[0], dims[1], dims[2]};
  }


  template<class S> 
  ptr_tuple_CS<S> get_(int i, ptr_tuple_CS<S>*ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 2)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 2 dimensional\n");

    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)){
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");      
    }    
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    std::pair<S *, S *> pointer_pair = {reinterpret_cast<S *>(mxGetPr(prhs[i])), reinterpret_cast<S *>(mxGetPi(prhs[i]))};
    return {pointer_pair, mxGetM(prhs[i]), mxGetN(prhs[i])};
  }

#else

  //Cases where we extract interleaved complex matrix to complex Eigen Map
  template <class S,typename U = typename S::Scalar>
  std::enable_if_t<isEigenMap<S>::value && (Eigen::NumTraits<U>::IsComplex==1),S > get_(int i, S *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 2)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 2 dimensional\n");
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<typename Eigen::NumTraits<U>::Real>::mxClass)){
      std::string identifier{mxClassTraits<typename Eigen::NumTraits<U>::Real>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");      
    }    
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    U* mex_ptr;
    if constexpr (std::is_same<typename Eigen::NumTraits<U>::Real,double>::value)
		   mex_ptr = reinterpret_cast<U*>(mxGetComplexDoubles(prhs[i]));
    if constexpr (std::is_same<typename Eigen::NumTraits<U>::Real,float>::value)
		   mex_ptr = reinterpret_cast<U*>(mxGetComplexSingles(prhs[i]));		   		   		 		   
    S new_map(mex_ptr, mxGetM(prhs[i]), mxGetN(prhs[i]));
    return new_map;
  }


  //Cases where we extract complex scalar (interleaved)
  template<class S>
  std::complex<S> get_(int i, std::complex<S> *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (std::any_of(dims,dims+ndims,[](auto x) {return x!=1;}))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a scalar\n");

    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)){
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");      
    }    
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    if constexpr (std::is_same<S,float>::value)
		   return reinterpret_cast<std::complex<float> *>(mxGetComplexSingles(prhs[i]))[0];
    if constexpr (std::is_same<S,double>::value)
		   return reinterpret_cast<std::complex<double> *>(mxGetComplexDoubles(prhs[i]))[0];    
  }


  //Cases where we extract complex matrix (interleaved) to pointer tuple 
  template<class S> 
  ptr_tuple<std::complex<S>> get_(int i, ptr_tuple<std::complex<S>> *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 2)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 2 dimensional\n");    

    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)){
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");      
    }    
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    std::complex<S> * mex_pointer;
    if constexpr (std::is_same<S,double>::value)
		   mex_pointer = reinterpret_cast<std::complex<S> *>(mxGetComplexDoubles(prhs[i]));
    if constexpr (std::is_same<S,float>::value)
		   mex_pointer = reinterpret_cast<std::complex<S> *>(mxGetComplexSingles(prhs[i]));
    return {mex_pointer, mxGetM(prhs[i]), mxGetN(prhs[i])};
  }

  //Cases where we extract 3 dim complex array (interleaved) to pointer tuple 
  template <typename S> ptr_tuple_3dim<std::complex<S>> get_(int i, ptr_tuple_3dim<std::complex<S>> *ignored) {
    mwSize ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 3)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 3 dimensional\n");
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string(" array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    if constexpr (std::is_same<S, float>::value)
      return {reinterpret_cast<std::complex<S> *>(mxGetComplexSingles(prhs[i])), dims[0], dims[1], dims[2]};
    if constexpr (std::is_same<S, double>::value)
      return {reinterpret_cast<std::complex<S> *>(mxGetComplexDoubles(prhs[i])), dims[0], dims[1], dims[2]};
  }

  
#endif
  
  
  // Recursively use unpacker to handle structs
  template <class S, std::size_t... Is> std::tuple<boost::pfr::tuple_element_t<Is, S>...> recurseUnpack(std::index_sequence<Is...>, int i, S *ignored) {    
    int num_fields = mxGetNumberOfFields(prhs[i]);
    std::unique_ptr<mxArray *[]> sub_rhs(new mxArray *[num_fields]);
    for (int field_ind = 0; field_ind < num_fields; field_ind++) {
      sub_rhs[field_ind] = mxGetFieldByNumber(prhs[i], 0, field_ind);
    }
    MexUnpacker<boost::pfr::tuple_element_t<Is, S>...> sub_unpack(num_fields, const_cast<const mxArray **>(sub_rhs.get()));
    return sub_unpack.unpackMex();
  }

  template <class S>
  std::enable_if_t<!in_type_list<S, EDRM, EDCIM, EDCSM, EDR, RDP, CDIP, CDSP, EFRM, EFCIM, EFCSM, EFR, RFP, CFIP, CFSP>()
   && !std::is_scalar<S>::value && !is_tuple<S>::value,
//   && !is_specialization<S,stdex::mdspan >::value, 
   S> 
   get_(int i,
                                                                                                                                                                                       S *ignored) {
    if (!mxIsStruct(prhs[i]))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a struct\n");
    S *S_ignored = nullptr;
    return make_struct<S>(recurseUnpack(std::make_index_sequence<boost::pfr::tuple_size<S>::value>(), i, S_ignored));
  }

  // Recursively use unpacker to handle struct arrays
  template <class S, std::size_t... Is> std::vector<S> recurseUnpackVec(std::index_sequence<Is...>, int i, S *ignored) {
    std::size_t num_row = mxGetM(prhs[i]);
    std::size_t num_col = mxGetN(prhs[i]);

    if (num_row > 1 && num_col > 1)
      throw std::string("Argument ") + std::to_string(i) + std::string(" is 2D struct array, not yet supported\n");

    num_row = num_row * num_col;

    int num_fields = mxGetNumberOfFields(prhs[i]);
    std::unique_ptr<mxArray *[]> sub_rhs(new mxArray *[num_fields]);
    std::vector<S> output_struct_vec;

    for (std::size_t row_ind = 0; row_ind < num_row; row_ind++) {
      for (int field_ind = 0; field_ind < num_fields; field_ind++) {
        sub_rhs[field_ind] = mxGetFieldByNumber(prhs[i], row_ind, field_ind);
      }
      MexUnpacker<boost::pfr::tuple_element_t<Is, S>...> sub_unpack(num_fields, const_cast<const mxArray **>(sub_rhs.get()));
      output_struct_vec.push_back(make_struct<S>(sub_unpack.unpackMex()));
    }
    return output_struct_vec;
  }

  //  template <class S> std::enable_if_t<!in_type_list<S, EDRM, EDCIM, EDCSM, EDR, RDP, CDIP, CDSP, EFRM, EFCIM, EFCSM, EFR, RFP, CFIP, CFSP>(), std::vector<S>> get_(int i, std::vector<S> *ignored) {
  template <class S> std::vector<S> get_(int i, std::vector<S> *ignored) {
    if (!mxIsStruct(prhs[i]))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a struct\n");
    S *S_ignored = nullptr;
    return recurseUnpackVec(std::make_index_sequence<boost::pfr::tuple_size<S>::value>(), i, S_ignored);
  }


  //Cell arrays map to vector of variants
  std::vector<std::string> get_(int i, std::vector<std::string> *ignored) {
    if (!mxIsCell(prhs[i]))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a Cell\n");

    std::size_t num_row = mxGetM(prhs[i]);
    std::size_t num_col = mxGetN(prhs[i]);

    if (num_row > 1 && num_col > 1)
      throw std::string("Argument ") + std::to_string(i) + std::string(" is 2D cell array, not yet supported\n");

    num_row = num_row * num_col;

//    int num_fields = mxGetNumberOfFields(prhs[i]);
//    std::unique_ptr<mxArray *[]> sub_rhs(new mxArray *[num_fields]);
//    std::vector<S> output_struct_vec;
    std::vector<std::string> output_string;

    for (std::size_t row_ind = 0; row_ind < num_row; row_ind++) {
        mxArray* cur_cell = mxGetCell(prhs[i],row_ind);
        if (!(mxGetClassID(cur_cell) == mxCHAR_CLASS))
           throw std::string("Argument ") + std::to_string(i) + std::string(" not an string\n");
        char *c_string = mxArrayToString(cur_cell);
        output_string.push_back(std::string(c_string));

    }
    return output_string;
  }


  // Calls get_ with index of each type and
  // (Null) pointer to appropriate type to select the correect overload (probably more elegant way to do this)
  template <int i> typename NthElement<i, T...>::type get() {
    if (i >= nrhs)
      throw std::string("Can't extract argument ") + std::to_string(i) + std::string(", only ") + std::to_string(nrhs) + std::string(" arguments on right side.\n");
    typename NthElement<i, T...>::type *ignored = nullptr;
    return get_(i, ignored);
  }

  // Calls get for each index in type list parameter pack
  template <std::size_t... Is> std::tuple<T...> unpackIndirect(std::index_sequence<Is...>) { return {get<Is>()...}; }

  // User function to unpack all the inputs from Matlab/Octave
  // Need to get index of each type in the pack which is dones by calling unpackIndirect
  std::tuple<T...> unpackMex() { return unpackIndirect(std::make_index_sequence<sizeof...(T)>()); }
};

// Needed to evaluate arguments in a parameter pack for MexPacker
template <class... T> void ignore(T &&...) {}

// Class to return data from c++ mex function back to matlab / octave
/*
Usage:
  MexPacker<double,int,Eigen::MatrixXcd,Eigen::MatrixXf,RDP> my_pack(nlhs,plhs);
  my_pack.PackMex(3.8,2,d,f,g);
*/

template <typename... T> class MexPacker {
public:
  int nlhs;
  mxArray **plhs;
  std::map<std::type_index, std::vector<std::string>> field_name_map;
  std::map<std::type_index, std::vector<const char *>> field_name_map_cstr;

  MexPacker(int nlhs_, mxArray *plhs_[]) : nlhs{nlhs_}, plhs{plhs_} {}

  // Constructor with optional map for struct array field names
  MexPacker(int nlhs_, mxArray *plhs_[], std::map<std::type_index, std::vector<std::string>> user_map) : nlhs{nlhs_}, plhs{plhs_} {

    for (auto [key, val] : user_map) {
      field_name_map[key] = val;
      std::vector<const char *> strings;

      for (int j = 0; j < val.size(); ++j)
        strings.push_back(field_name_map[key][j].c_str());
      field_name_map_cstr[key] = strings;
    }
  }


  //return a string
  template <int i> int put(const std::string &arg) {
    plhs[i] = mxCreateString(arg.c_str());
    return 0;
  }

  // Handle all real scalars (no imaginary component)
  template <int i, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const S &arg) {

    mwSize dims[2] = {1, 1};

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *cur_pointer = reinterpret_cast<S *>(mxGetData(m));
    cur_pointer[0] = arg;
    plhs[i] = m;
    return 0;
  }

  // Handle all real pointer tuples (no imaginary component) for 2 dim arrays
  template <int i, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const ptr_tuple<S> &arg) {

    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *mex_pointer = reinterpret_cast<S *>(mxGetData(m));
    S *return_pointer = std::get<0>(arg);
    std::copy_n(return_pointer,dims[0]*dims[1],mex_pointer);
    plhs[i] = m;

    return 0;
  }

  // Handle all real pointer tuples (no imaginary component) for 3 dim arrays
  template <int i, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const ptr_tuple_3dim<S> &arg) {
    mwSize dims[3];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);
    dims[2] = std::get<3>(arg);

    mxArray *m = mxCreateNumericArray(3, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *mex_pointer = reinterpret_cast<S *>(mxGetData(m));
    S *return_pointer = std::get<0>(arg);
    std::copy_n(return_pointer,dims[0]*dims[1]*dims[2],mex_pointer);
    plhs[i] = m;
    return 0;
  }



#ifdef USE_MDSPAN

#ifdef COMPLEX_SPLIT

template <int n, typename S,size_t... U,typename W> std::enable_if_t<std::is_scalar<S>::value, int> put(const stdex::mdspan<S, stdex::extents<U...>, W> &arg) {
  mwSize dims[stdex::extents<U...>::rank()];

  for(size_t i = 0; i<stdex::extents<U...>::rank(); i++){
    dims[i] = arg.extent(i);
  }

  mxArray *m = mxCreateNumericArray(stdex::extents<U...>::rank(), dims, mxClassTraits<S>::mxClass, mxREAL);
  S *mex_pointer = reinterpret_cast<S *>(mxGetData(m));

//  stdex::mdspan<S, stdex::extents<U...>, W> m_mdspan(mex_pointer,arg.extents());
  stdex::mdspan<S, stdex::extents<U...>, stdex::layout_left> m_mdspan(mex_pointer,arg.extents());

  if constexpr (stdex::extents<U...>::rank()==1 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      m_mdspan(i) = arg(i);
    }
  }

  if constexpr (stdex::extents<U...>::rank()==2 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      for(size_t j = 0; j<arg.extent(1);j++){
        m_mdspan(i,j) = arg(i,j);
      }
    }
  }

  if constexpr (stdex::extents<U...>::rank()==3 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      for(size_t j = 0; j<arg.extent(1);j++){
        for(size_t k = 0; k<arg.extent(2);k++){
            m_mdspan(i,j,k) = arg(i,j,k);
        }
      }
    }
  }

  if constexpr (stdex::extents<U...>::rank()==4 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      for(size_t j = 0; j<arg.extent(1);j++){
        for(size_t k = 0; k<arg.extent(2);k++){
          for(size_t l = 0; l<arg.extent(3);l++){
            m_mdspan(i,j,k,l) = arg(i,j,k,l);
          }
        }
      }
    }
  }
  plhs[n] = m;
  return 0;
}

template <int n, typename S,size_t... U,typename W> 
std::enable_if_t<std::is_same<S, double>::value ||std::is_same<S, float>::value,int>
put(const stdex::mdspan<std::complex<S>, stdex::extents<U...>, W> &arg) {
  mwSize dims[stdex::extents<U...>::rank()];

  for(size_t i = 0; i< stdex::extents<U...>::rank(); i++){
    dims[i] = arg.extent(i);
  }

//  mxArray *m = mxCreateNumericArray(stdex::extents<U...>::rank(), dims, mxClassTraits<S>::mxClass, mxCOMPLEX);
//  S *mex_pointer = reinterpret_cast<S *>(mxGetData(m));


  mxArray *m = mxCreateNumericArray(stdex::extents<U...>::rank(), dims, mxClassTraits<S>::mxClass, mxCOMPLEX);
  S *mex_pointer_real = reinterpret_cast<S *>(mxGetPr(m));
  S *mex_pointer_imag = reinterpret_cast<S *>(mxGetPi(m));


  stdex::mdspan<S, stdex::extents<U...>, stdex::layout_left> m_mdspan_r(mex_pointer_real,arg.extents());
  stdex::mdspan<S, stdex::extents<U...>, stdex::layout_left> m_mdspan_i(mex_pointer_imag,arg.extents());

  if constexpr (stdex::extents<U...>::rank()==1 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      m_mdspan_r(i) = arg(i).real();
      m_mdspan_i(i) = arg(i).imag();
    }
  }

  if constexpr (stdex::extents<U...>::rank()==2 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      for(size_t j = 0; j<arg.extent(1);j++){
        m_mdspan_r(i,j) = arg(i,j).real();
        m_mdspan_i(i,j) = arg(i,j).imag();
      }
    }
  }

  if constexpr (stdex::extents<U...>::rank()==3 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      for(size_t j = 0; j<arg.extent(1);j++){
        for(size_t k = 0; k<arg.extent(2);k++){
          m_mdspan_r(i,j,k) = arg(i,j,k).real();
          m_mdspan_i(i,j,k) = arg(i,j,k).imag();
        }
      }
    }
  }

  if constexpr (stdex::extents<U...>::rank()==4 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      for(size_t j = 0; j<arg.extent(1);j++){
        for(size_t k = 0; k<arg.extent(2);k++){
          for(size_t l = 0; l<arg.extent(3);l++){
            m_mdspan_r(i,j,k,l) = arg(i,j,k,l).real();
            m_mdspan_i(i,j,k,l) = arg(i,j,k,l).imag();
          }
        }
      }
    }
  }
  plhs[n] = m;
  return 0;
}


#else
     
template<int n  , typename S,size_t... U,typename W> 
//template<int n  , template <typename S,size_t... U,typename W>   stdex::mdspan<S, stdex::extents<U...>, W> >
int put(const stdex::mdspan<S, stdex::extents<U...>, W> &arg) {
  mwSize dims[stdex::extents<U...>::rank()];

  for(size_t i = 0; i<stdex::extents<U...>::rank(); i++){
    dims[i] = arg.extent(i);
  }
  

  mxArray *m;
  if constexpr (!is_complex<S>::value){
    m = mxCreateNumericArray(stdex::extents<U...>::rank(), dims, mxClassTraits<S>::mxClass, mxREAL);
  }else{
    static_assert(std::is_same<typename is_complex<S>::real_type,float>::value||std::is_same<typename is_complex<S>::real_type,double>::value,
    "Complex arrays must be complex floats or complex doubles");
    m = mxCreateNumericArray(stdex::extents<U...>::rank(), dims, mxClassTraits<typename is_complex<S>::real_type>::mxClass, mxCOMPLEX);
  }
  S *mex_pointer = reinterpret_cast<S *>(mxGetData(m));

  stdex::mdspan<S, stdex::extents<U...>, stdex::layout_left> m_mdspan(mex_pointer,arg.extents());

  if constexpr (stdex::extents<U...>::rank()==1 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      m_mdspan(i) = arg(i);
    }
  }

  if constexpr (stdex::extents<U...>::rank()==2 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      for(size_t j = 0; j<arg.extent(1);j++){
        m_mdspan(i,j) = arg(i,j);
      }
    }
  }

  if constexpr (stdex::extents<U...>::rank()==3 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      for(size_t j = 0; j<arg.extent(1);j++){
        for(size_t k = 0; k<arg.extent(2);k++){
            m_mdspan(i,j,k) = arg(i,j,k);
        }
      }
    }
  }

  if constexpr (stdex::extents<U...>::rank()==4 ){
    for(size_t i = 0; i<arg.extent(0); i++){
      for(size_t j = 0; j<arg.extent(1);j++){
        for(size_t k = 0; k<arg.extent(2);k++){
          for(size_t l = 0; l<arg.extent(3);l++){
            m_mdspan(i,j,k,l) = arg(i,j,k,l);
          }
        }
      }
    }
  }
  plhs[n] = m;
  return 0;
}

#endif
#endif

#ifdef COMPLEX_SPLIT

  //handle 3 dim arrays with split real/imaginary
  template <int i, typename S> int put(const ptr_tuple_3dim_CS<S> &arg) {

    mwSize dims[3];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);
    dims[2] = std::get<3>(arg);

    mxArray *m = mxCreateNumericArray(3, dims, mxClassTraits<S>::mxClass, mxCOMPLEX);
    S *mex_pointer_real = reinterpret_cast<S *>(mxGetPr(m));
    S *mex_pointer_imag = reinterpret_cast<S *>(mxGetPi(m));

    S *return_pointer_real = std::get<0>(arg).first;
    S *return_pointer_imag = std::get<0>(arg).second;

    std::copy_n(return_pointer_real,dims[0]*dims[1]*dims[2],mex_pointer_real);
    std::copy_n(return_pointer_imag,dims[0]*dims[1]*dims[2],mex_pointer_imag);    
    plhs[i] = m;
    return 0;
  }

#else
  //handle 3 dim arrays with interleaved real/imaginary
  template <int i, typename S> std::enable_if_t<std::is_same<S, float>::value || std::is_same<S, double>::value, int> put(const ptr_tuple_3dim<std::complex<S>> &arg) {

    mwSize dims[3];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);
    dims[2] = std::get<3>(arg);

    mxArray *m = mxCreateNumericArray(3, dims, mxClassTraits<S>::mxClass, mxCOMPLEX);

    std::complex<S> *mex_pointer = nullptr;
    if constexpr (std::is_same<S, double>::value)
      mex_pointer = reinterpret_cast<std::complex<S> *>(mxGetComplexDoubles(m));
    if constexpr (std::is_same<S, float>::value)
      mex_pointer = reinterpret_cast<std::complex<S> *>(mxGetComplexSingles(m));

    std::complex<S> *return_pointer = std::get<0>(arg);

    std::copy_n(return_pointer,dims[0]*dims[1]*dims[2],mex_pointer);
    
    plhs[i] = m;

    return 0;
  }

#endif


  //Real ptr tuple (unique pointer)
  template <int i, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const unique_ptr_tuple<S> &arg) {

    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *mex_pointer = reinterpret_cast<S *>(mxGetData(m));
    //    S *return_pointer = std::get<0>(arg);
    std::copy_n(std::get<0>(arg).get(),dims[0]*dims[1],mex_pointer);
    plhs[i] = m;

    return 0;
  }

  //Real ptr tuple (shared pointer)  
  template <int i, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const shared_ptr_tuple<S> &arg) {

    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *mex_pointer = reinterpret_cast<S *>(mxGetData(m));
    std::copy_n(std::get<0>(arg).get(),dims[0]*dims[1],mex_pointer);
    plhs[i] = m;

    return 0;
  }


  //Real array
  template <int i, std::size_t N, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const S (&arg)[N]) {

    mwSize dims[2];
    dims[0] = N;
    dims[1] = 1;

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *mex_pointer = reinterpret_cast<S *>(mxGetData(m));
    std::copy_n(arg,dims[0]*dims[1],mex_pointer);
    plhs[i] = m;
    return 0;
  }


  //Handle real Eigen matrix cases
  template <int i, typename Derived, typename U = typename Derived::Scalar, typename = std::enable_if_t<std::is_same<U, float>::value||std::is_same<U, double>::value>>
  int put(const Eigen::MatrixBase<Derived> &arg){
    
    mwSize dims[2];
    dims[0] = arg.rows();
    dims[1] = arg.cols();

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<U>::mxClass, mxREAL);
    U *mex_pointer = reinterpret_cast<U *>(mxGetData(m));
    plhs[i] = m;

    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        mex_pointer[jj * dims[0] + ii] = arg(ii, jj);
      }
    }
    return 0;
  }


  //Handle complex Eigen matrix cases  
  template <int i, typename Derived, typename U = typename Derived::Scalar,typename W = typename Eigen::NumTraits<U>::Real,
	    typename = std::enable_if_t<std::is_same<U, std::complex<double>>::value ||std::is_same<U, std::complex<float>>::value>>
  int put(const Eigen::MatrixBase<Derived> &arg) {

    mwSize dims[2];
    dims[0] = arg.rows();
    dims[1] = arg.cols();

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<W>::mxClass, mxCOMPLEX);

#ifdef COMPLEX_SPLIT

    
    W *mex_pointer_real = reinterpret_cast<W *>(mxGetPr(m));
    W *mex_pointer_imag = reinterpret_cast<W *>(mxGetPi(m));

    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {

        mex_pointer_real[jj * dims[0] + ii] = arg(ii, jj).real();
        mex_pointer_imag[jj * dims[0] + ii] = arg(ii, jj).imag();
      }
    }

#else

    std::complex<W> *mex_pointer;
    if constexpr (std::is_same<W,double>::value)
		   mex_pointer = reinterpret_cast<std::complex<W> *>(mxGetComplexDoubles(m));
    
    if constexpr (std::is_same<W,float>::value)
		   mex_pointer = reinterpret_cast<std::complex<W> *>(mxGetComplexSingles(m));

    
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        mex_pointer[jj * dims[0] + ii] = arg(ii, jj);
      }
    }

#endif

    plhs[i] = m;

    return 0;
  }



  //Handle complex pointer tuple cases
#ifdef COMPLEX_SPLIT

  template <int i,typename S> int put(const ptr_tuple_CS<S>&arg) {

    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    S *return_real = std::get<0>(arg).first;
    S *return_imag = std::get<0>(arg).second;

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxCOMPLEX);
    S *mex_pointer_real = reinterpret_cast<S *>(mxGetPr(m));
    S *mex_pointer_imag = reinterpret_cast<S *>(mxGetPi(m));
    std::copy_n(return_real,dims[0]*dims[1],mex_pointer_real);
    std::copy_n(return_imag,dims[0]*dims[1],mex_pointer_imag);    
    plhs[i] = m;
    return 0;
  }

#else

  template <int i,class S> 
  std::enable_if_t<std::is_same<S,std::complex<double>>::value||std::is_same<S,std::complex<float>>::value, int> put(const  ptr_tuple<S>&arg) {
    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);
    typedef typename Eigen::NumTraits<S>::Real U;

    std::complex<U> *return_complex = std::get<0>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<U>::mxClass, mxCOMPLEX);
    std::complex<U> *mex_pointer;
    if constexpr(std::is_same<U,double>::value) 
      mex_pointer= reinterpret_cast<std::complex<double> *>(mxGetComplexDoubles(m));
    if constexpr(std::is_same<U,float>::value) 
      mex_pointer= reinterpret_cast<std::complex<float> *>(mxGetComplexSingles(m));
    std::copy_n(return_complex,dims[0]*dims[1],mex_pointer);
    plhs[i] = m;
    return 0;
  }

#endif



  //Handle complex scalar cases
  template <int i, typename S> 
  std::enable_if_t<std::is_same<S, std::complex<double>>::value||std::is_same<S,std::complex<float>>::value, int> put(const S &arg) {
    mwSize dims[2];
    dims[0] = 1;
    dims[1] = 1;
    typedef typename Eigen::NumTraits<S>::Real U;
    
    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<U>::mxClass, mxCOMPLEX);

    
#ifdef COMPLEX_SPLIT

    U *mex_pointer_real = reinterpret_cast<U*>(mxGetPr(m));
    U *mex_pointer_imag = reinterpret_cast<U*>(mxGetPi(m));

    mex_pointer_real[0] = arg.real();
    mex_pointer_imag[0] = arg.imag();

#else

    std::complex<U> *mex_pointer;
    if constexpr (std::is_same<U,double>::value)
		   mex_pointer = reinterpret_cast<std::complex<U> *>(mxGetComplexDoubles(m));
    if constexpr (std::is_same<U,float>::value)
		   mex_pointer = reinterpret_cast<std::complex<U> *>(mxGetComplexSingles(m));    
    mex_pointer[0] = arg;

#endif

    plhs[i] = m;
    return 0;
  }



  // Get fieldnames for a give user struct type
  template <class S> const char **get_names() {
    if (field_name_map.find(typeid(S)) == field_name_map.end()) {
      int num_fields = boost::pfr::tuple_size<S>::value;
      std::vector<std::string> field_names;
      std::vector<const char *> strings;
      for (int j = 0; j < num_fields; ++j)
        field_names.push_back("field_" + std::to_string(j));

      field_name_map[typeid(S)] = field_names;

      for (int j = 0; j < num_fields; ++j)
        strings.push_back(field_name_map[typeid(S)][j].c_str());

      field_name_map_cstr[typeid(S)] = strings;
    }
    return field_name_map_cstr[typeid(S)].data();
  }

  // Recursively call Pack to handle structs
  template <class S, std::size_t... Is> void recursePack(std::index_sequence<Is...>, int i, S &arg) {

    int num_fields = boost::pfr::tuple_size<S>::value;

    plhs[i] = mxCreateStructMatrix(1, 1, num_fields, get_names<S>());

    std::unique_ptr<mxArray *[]> sub_lhs(new mxArray *[num_fields]);

    MexPacker<boost::pfr::tuple_element_t<Is, S>...> sub_pack(num_fields, sub_lhs.get(), field_name_map);
    sub_pack.PackMex(boost::pfr::get<Is, S>(arg)...);
    for (int field_ind = 0; field_ind < num_fields; field_ind++) {
      mxSetFieldByNumber(plhs[i], 0, field_ind, sub_lhs[field_ind]);
    }
  }

  //Assumes any unexpected type is a user struct to be returned as a struct
  template <int i, class S = void>
  std::enable_if_t<!in_type_list<S, EDRM, EDCIM, EDCSM, EDR, RDP, CDIP, CDSP, EFRM, EFCIM, EFCSM, EFR, RFP, CFIP, CFSP, EDC, EFC, std::complex<float>, std::complex<double>>() &&
                       !std::is_scalar<S>::value && !is_tuple<S>::value && 
                       !is_base_of_template<Eigen::MatrixBase,S>::value,
//                       &&!is_specialization<S,stdex::mdspan >::value,
                   int>
  put(const S &arg) {
    recursePack(std::make_index_sequence<boost::pfr::tuple_size<S>::value>(), i, arg);
    return 0;
  }

  // Recursively call Pack to handle struct arrays
  template <class S, std::size_t... Is> void recursePackVec(std::index_sequence<Is...>, int i, const std::vector<S> &arg) {

    int num_fields = boost::pfr::tuple_size<S>::value;
    plhs[i] = mxCreateStructMatrix(arg.size(), 1, num_fields, get_names<S>());

    std::unique_ptr<mxArray *[]> sub_lhs(new mxArray *[num_fields]);
    //     S tmp;
    for (std::size_t row_ind = 0; row_ind < arg.size(); row_ind++) {

      MexPacker<boost::pfr::tuple_element_t<Is, S>...> sub_pack(num_fields, sub_lhs.get(), field_name_map);
      sub_pack.PackMex(boost::pfr::get<Is, S>(arg[row_ind])...);
      for (int field_ind = 0; field_ind < num_fields; field_ind++) {
        mxSetFieldByNumber(plhs[i], row_ind, field_ind, sub_lhs[field_ind]);
      }
    }
  }

  // Handle struct arrays
  //  template <int i, class S = void> std::enable_if_t<!in_type_list<S, EDRM, EDCIM, EDCSM, EDR, RDP, CDIP, CDSP, EFRM, EFCIM, EFCSM, EFR, RFP, CFIP, CFSP,std::complex<float>,std::complex<double>
  //  >(), int>
  template <int i, typename S> std::enable_if_t<!std::is_scalar<S>::value, int>
  //template <int i, class S> int 
  put(const std::vector<S> &arg) {
    recursePackVec(std::make_index_sequence<boost::pfr::tuple_size<S>::value>(), i, arg);
    return 0;
  }



   //cell arrays
  template <int i> 
  void put(const std::vector<std::string> &arg) {
    std::size_t num_elements = arg.size();
    mxArray *cell_ptr = mxCreateCellMatrix(num_elements, 1);
    plhs[i] = cell_ptr;

    for (std::size_t j = 0; j < num_elements; ++j) {
      mxArray *cur_cell=mxCreateString(arg[j].c_str());
      mxSetCell(cell_ptr, j, cur_cell);
    }
//    return;
  }




  template <int i, typename S> std::enable_if_t<std::is_scalar<S>::value, int>
  put(const std::vector<S> &arg) {
    mwSize dims[2];
    dims[0] = arg.size();
    dims[1] = 1;

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *mex_pointer = reinterpret_cast<S *>(mxGetData(m));
    std::copy_n(arg.begin(),dims[0],mex_pointer);
    plhs[i] = m;

    return 0;
  }

  template <int i, typename S>
  std::enable_if_t<std::is_same<S, float>::value||std::is_same<S, double>::value,int>
  put(const std::vector<std::complex<S>> &arg) {
    mwSize dims[2];
    dims[0] = arg.size();
    dims[1] = 1;
    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxCOMPLEX);


#ifdef COMPLEX_SPLIT

    S *mex_pointer_real = reinterpret_cast<S *>(mxGetPr(m));
    S *mex_pointer_imag = reinterpret_cast<S *>(mxGetPi(m));
    std::transform(arg.begin(),arg.end(),mex_pointer_real,[](auto a){return a.real();});
    std::transform(arg.begin(),arg.end(),mex_pointer_imag,[](auto a){return a.imag();});  

#else

    std::complex<S> *mex_pointer;
    if constexpr(std::is_same<S,double>::value) 
      mex_pointer= reinterpret_cast<std::complex<double> *>(mxGetComplexDoubles(m));
    if constexpr(std::is_same<S,float>::value) 
      mex_pointer= reinterpret_cast<std::complex<float> *>(mxGetComplexSingles(m));
    std::copy_n(arg.begin(),dims[0]*dims[1],mex_pointer);  

#endif


    plhs[i] = m;
    return 0;
  }


  // Variants are turned into cell arrays
  template <int i, class... S> int put(const std::vector<std::variant<S...>> &arg) {
    std::size_t num_elements = arg.size();
    mxArray *cell_ptr = mxCreateCellMatrix(num_elements, 1);
    plhs[i] = cell_ptr;

    for (std::size_t j = 0; j < num_elements; ++j) {

      mxArray *cur_cell;
      auto pack_visitor = [&cur_cell, this](const auto &a) {
        MexPacker<decltype(a)> sub_pack(1, &cur_cell, field_name_map);
        sub_pack.PackMex(a);
      };
      std::visit(pack_visitor, arg[j]);
      mxSetCell(cell_ptr, j, cur_cell);
    }
    return 0;
  }

  template <int i, typename S> int check_and_put(const S &arg) {
    if (i >= nlhs) {
      std::string s = std::string("Can't return argument ") + std::to_string(i) + std::string(", only ") + std::to_string(nlhs) + std::string(" arguments on right side.\n");
      mexPrintf(s.data());
    } else {
      put<i>(arg);
    }
    return 0;
  }

  // ignore does nothing but lets us expand the parameter pack to evaluate all the put functions for each type
  template <std::size_t... Is> void packIndirect(const T &... args, std::index_sequence<Is...>) { ignore(check_and_put<Is, T>(args)...); }

  // Need to get index of each type on the pack
  void PackMex(const T &... args) { return packIndirect(args..., std::make_index_sequence<sizeof...(T)>()); }
};

#endif
