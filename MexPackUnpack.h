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

// COMPLEX_SPLIT determines whether complex data viewed as interleaved or split
// since matlab 2018a+  (when compiled with -R2018a) api for complex data requires mex functions not present in the octave or earlier matlab verisons mex.h. We have #ifdefs in parts of code for
// complex data to let the code work with earlier or later versions (this could probably be done better). This version of the code has not actually been tested with R2018a (due to not having a recent
// enough version of matlab to test on, so using with R2018a may require some tweaks)

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

template <typename> struct is_tuple : std::false_type {};
template <typename... T> struct is_tuple<std::tuple<T...>> : std::true_type {};

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
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string("\n");
    }
    return static_cast<S>(mxGetScalar(prhs[i]));
  }

  // Deal with all pointer tuples that are real (no imaginary component)
  template <typename S> std::enable_if_t<std::is_scalar<S>::value, ptr_tuple<S>> get_(int i, ptr_tuple<S> *ignored) {
    //   get_(int i, RFP *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string("array\n");
    }
    //      throw std::string("Argument ") + std::to_string(i) + std::string(" not a single array\n");
    return {reinterpret_cast<S *>(mxGetPr(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i])};
  }

  template <typename S> std::enable_if_t<std::is_scalar<S>::value, ptr_tuple_3dim<S>> get_(int i, ptr_tuple_3dim<S> *ignored) {
    int ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 3)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 3 dimensional\n");

    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string("array\n");
    }
    return {reinterpret_cast<S *>(mxGetPr(prhs[i])), dims[0], dims[1], dims[2]};
  }

#ifdef COMPLEX_SPLIT

  template <typename S> ptr_tuple_3dim_CI<S> get_(int i, ptr_tuple_3dim_CI<S> *ignored) {
    int ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 3)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 3 dimensional\n");
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string("array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    return {std::pair{reinterpret_cast<S *>(mxGetPr(prhs[i])), reinterpret_cast<S *>(mxGetPi(prhs[i]))}, dims[0], dims[1], dims[2]};
  }

#else

  template <typename S> ptr_tuple_3dim<std::complex<S>> get_(int i, ptr_tuple_3dim<std::complex<S>> *ignored) {
    int ndims = mxGetNumberOfDimensions(prhs[i]);
    if (ndims != 3)
      throw std::string("Argument ") + std::to_string(i) + std::string(" not 3 dimensional\n");
    const mwSize *dims = mxGetDimensions(prhs[i]);
    if (!(mxGetClassID(prhs[i]) == mxClassTraits<S>::mxClass)) {
      std::string identifier{mxClassTraits<S>::name};
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a ") + identifier + std::string("array\n");
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

  std::string get_(int i, std::string *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxCHAR_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not an string\n");
    char *c_string = mxArrayToString(prhs[i]);
    return std::string(c_string);
  }

  // Behavior of complex arrays depends on whether arrays are interleaved or split
  // Octave and MATLAB before 2018a use split real and complex
  // Matlab when compiled with -R2018a uses interleaved real and complex
  // We return Eigen Maps which wrap the underlying pointers from octave
  // If interleaved the Eigen Maps is complex, if split we return a pair of real eigen maps
  // Note: This has not yet been tested with Matlab 2018a (basic idea should work but may need to be tweaked)

#ifdef COMPLEX_SPLIT

  std::complex<double> get_(int i, std::complex<double> *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxDOUBLE_CLASS)) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a double array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    return std::complex<double>(static_cast<double *>(mxGetPr(prhs[i]))[0], static_cast<double *>(mxGetPi(prhs[i]))[0]);
  }

  EDCSM get_(int i, EDCSM *ignored) {

    if (!(mxGetClassID(prhs[i]) == mxDOUBLE_CLASS)) {

      throw std::string("Argument ") + std::to_string(i) + std::string(" not a double array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    Eigen::Map<Eigen::MatrixXd> new_map_real(static_cast<double *>(mxGetPr(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    Eigen::Map<Eigen::MatrixXd> new_map_imag(static_cast<double *>(mxGetPi(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    return {new_map_real, new_map_imag};
  }

  CDSP get_(int i, CDSP *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxDOUBLE_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a double array\n");
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    std::pair<double *, double *> pointer_pair = {static_cast<double *>(mxGetPr(prhs[i])), static_cast<double *>(mxGetPi(prhs[i]))};
    return {pointer_pair, mxGetM(prhs[i]), mxGetN(prhs[i])};
  }

#else

  std::complex<double> get_(int i, std::complex<double> *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxDOUBLE_CLASS)) {

      throw std::string("Argument ") + std::to_string(i) + std::string(" not a double array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    return reinterpret_cast<std::complex<double> *>(mxGetComplexDoubles(prhs[i]))[0];
  }

  EDCIM get_(int i, EDCIM *ignored) {

    if (!(mxGetClassID(prhs[i]) == mxDOUBLE_CLASS)) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a double array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    Eigen::Map<Eigen::MatrixXcd> new_map(reinterpret_cast<std::complex<double> *>(mxGetComplexDoubles(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    return new_map;
  }

  CDIP get_(int i, CDIP *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxDOUBLE_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a double array\n");
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    std::complex<double> *complex_pointer = reinterpret_cast<std::complex<double> *>(mxGetComplexDoubles(prhs[i]));
    return {complex_pointer, mxGetM(prhs[i]), mxGetN(prhs[i])};
  }

#endif

#ifdef COMPLEX_SPLIT

  std::complex<float> get_(int i, std::complex<float> *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxSINGLE_CLASS)) {

      throw std::string("Argument ") + std::to_string(i) + std::string(" not a single array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    return std::complex<float>(reinterpret_cast<float *>(mxGetPr(prhs[i]))[0], reinterpret_cast<float *>(mxGetPi(prhs[i]))[0]);
  }

  EFCSM get_(int i, EFCSM *ignored) {

    if (!(mxGetClassID(prhs[i]) == mxSINGLE_CLASS)) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a single array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    Eigen::Map<Eigen::MatrixXf> new_map_real(reinterpret_cast<float *>(mxGetPr(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    Eigen::Map<Eigen::MatrixXf> new_map_imag(reinterpret_cast<float *>(mxGetPi(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    return {new_map_real, new_map_imag};
  }

  CFSP get_(int i, CFSP *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxSINGLE_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a single array\n");
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    std::pair<float *, float *> pointer_pair = {reinterpret_cast<float *>(mxGetPr(prhs[i])), reinterpret_cast<float *>(mxGetPi(prhs[i]))};
    return {pointer_pair, mxGetM(prhs[i]), mxGetN(prhs[i])};
  }

#else

  std::complex<float> get_(int i, std::complex<float> *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxSINGLE_CLASS)) {

      throw std::string("Argument ") + std::to_string(i) + std::string(" not a single array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    return reinterpret_cast<std::complex<float> *>(mxGetComplexSingles(prhs[i]))[0];
  }

  EFCIM get_(int i, EFCIM *ignored) {

    if (!(mxGetClassID(prhs[i]) == mxSINGLE_CLASS)) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a single array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    Eigen::Map<Eigen::MatrixXcf> new_map(reinterpret_cast<std::complex<float> *>(mxGetComplexSingles(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    return new_map;
  }

  CFIP get_(int i, CFIP *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxSINGLE_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a single array\n");
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    std::complex<float> *complex_pointer = reinterpret_cast<std::complex<float> *>(mxGetComplexSingles(prhs[i]));
    return {complex_pointer, mxGetM(prhs[i]), mxGetN(prhs[i])};
  }

#endif

  Eigen::Map<Eigen::MatrixXf> get_(int i, Eigen::Map<Eigen::MatrixXf> *ignored) {

    if (!(mxGetClassID(prhs[i]) == mxSINGLE_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a single array\n");
    Eigen::Map<Eigen::MatrixXf> new_map(reinterpret_cast<float *>(mxGetPr(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    return new_map;
  }

  Eigen::Map<Eigen::MatrixXd> get_(int i, Eigen::Map<Eigen::MatrixXd> *ignored) {

    if (!(mxGetClassID(prhs[i]) == mxDOUBLE_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a double array\n");
    Eigen::Map<Eigen::MatrixXd> new_map(reinterpret_cast<double *>(mxGetPr(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    return new_map;
  }

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
  std::enable_if_t<!in_type_list<S, EDRM, EDCIM, EDCSM, EDR, RDP, CDIP, CDSP, EFRM, EFCIM, EFCSM, EFR, RFP, CFIP, CFSP>() && !std::is_scalar<S>::value && !is_tuple<S>::value, S> get_(int i,
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

  template <int i> int put(const std::string &arg) {
    plhs[i] = mxCreateString(arg.c_str());
    return 0;
  }

  // Handle all real scalars (no imaginary component)

  template <int i, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const S &arg) {
    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2] = {1, 1};

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *cur_pointer = reinterpret_cast<S *>(mxGetPr(m));
    cur_pointer[0] = arg;
    plhs[i] = m;
    return 0;
  }

  // Handle all real pointer tuples (no imaginary component)

  template <int i, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const ptr_tuple<S> &arg) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *cur_pointer = reinterpret_cast<S *>(mxGetPr(m));
    S *return_pointer = std::get<0>(arg);

    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = return_pointer[jj * dims[0] + ii];
      }
    }
    return 0;
  }

  template <int i, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const ptr_tuple_3dim<S> &arg) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[3];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);
    dims[2] = std::get<3>(arg);

    mxArray *m = mxCreateNumericArray(3, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *cur_pointer = reinterpret_cast<S *>(mxGetPr(m));
    S *return_pointer = std::get<0>(arg);

    plhs[i] = m;
    for (std::size_t kk = 0; kk < dims[0] * dims[1] * dims[2]; kk++) {
      cur_pointer[kk] = return_pointer[kk];
    }

    return 0;
  }

#ifdef COMPLEX_SPLIT

  template <int i, typename S> int put(const ptr_tuple_3dim_CI<S> &arg) {

    mwSize dims[3];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);
    dims[2] = std::get<3>(arg);

    mxArray *m = mxCreateNumericArray(3, dims, mxClassTraits<S>::mxClass, mxCOMPLEX);
    S *cur_pointer_real = reinterpret_cast<S *>(mxGetPr(m));
    S *cur_pointer_imag = reinterpret_cast<S *>(mxGetPi(m));

    S *return_pointer_real = std::get<0>(arg).first;
    S *return_pointer_imag = std::get<0>(arg).second;

    plhs[i] = m;
    for (std::size_t kk = 0; kk < dims[0] * dims[1] * dims[2]; kk++) {
      cur_pointer_real[kk] = return_pointer_real[kk];
      cur_pointer_imag[kk] = return_pointer_imag[kk];
    }

    return 0;
  }

#else

  template <int i, typename S> std::enable_if_t<std::is_same<S, float>::value || std::is_same<S, double>::value, int> put(const ptr_tuple_3dim<std::complex<S>> &arg) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[3];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);
    dims[2] = std::get<3>(arg);

    mxArray *m = mxCreateNumericArray(3, dims, mxClassTraits<S>::mxClass, mxCOMPLEX);

    std::complex<S> *cur_pointer = nullptr;
    if constexpr (std::is_same<S, double>::value)
      std::complex<S> *cur_pointer = reinterpret_cast<S *>(mxGetComplexDoubles(m));
    if constexpr (std::is_same<S, float>::value)
      std::complex<S> *cur_pointer = reinterpret_cast<S *>(mxGetComplexSingles(m));

    std::complex<S> *return_pointer = std::get<0>(arg);

    plhs[i] = m;
    for (std::size_t kk = 0; kk < dims[0] * dims[1] * dims[2]; kk++) {
      cur_pointer[kk] = return_pointer[kk];
    }

    return 0;
  }

#endif

  template <int i, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const unique_ptr_tuple<S> &arg) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *cur_pointer = reinterpret_cast<S *>(mxGetPr(m));
    //    S *return_pointer = std::get<0>(arg);

    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = std::get<0>(arg)[jj * dims[0] + ii];
      }
    }
    return 0;
  }

  template <int i, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const shared_ptr_tuple<S> &arg) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *cur_pointer = reinterpret_cast<S *>(mxGetPr(m));

    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = std::get<0>(arg)[jj * dims[0] + ii];
      }
    }
    return 0;
  }

  template <int i, std::size_t N, typename S> std::enable_if_t<std::is_scalar<S>::value, int> put(const S (&arg)[N]) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = N;
    dims[1] = 1;

    mxArray *m = mxCreateNumericArray(2, dims, mxClassTraits<S>::mxClass, mxREAL);
    S *cur_pointer = reinterpret_cast<S *>(mxGetPr(m));

    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = arg[jj * dims[0] + ii];
      }
    }
    return 0;
  }

  template <int i, typename Derived, typename = std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value>> int put(const Eigen::MatrixBase<Derived> &arg, double *ignored = nullptr) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = arg.rows();
    dims[1] = arg.cols();

    mxArray *m = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double *cur_pointer = static_cast<double *>(mxGetPr(m));
    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = arg(ii, jj);
      }
    }
    return 0;
  }

  template <int i, typename Derived, typename = std::enable_if_t<std::is_same<typename Derived::Scalar, float>::value>> int put(const Eigen::MatrixBase<Derived> &arg, float *x = nullptr) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = arg.rows();
    dims[1] = arg.cols();

    mxArray *m = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    float *cur_pointer = reinterpret_cast<float *>(mxGetPr(m));
    plhs[i] = m;

    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = arg(ii, jj);
      }
    }
    return 0;
  }

#ifdef COMPLEX_SPLIT

  template <int i> int put(const CDSP &arg) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    double *return_real = std::get<0>(arg).first;
    double *return_imag = std::get<0>(arg).second;

    mxArray *m = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
    double *cur_pointer_real = static_cast<double *>(mxGetPr(m));
    double *cur_pointer_imag = static_cast<double *>(mxGetPi(m));

    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer_real[jj * dims[0] + ii] = return_real[jj * dims[0] + ii];
        cur_pointer_imag[jj * dims[0] + ii] = return_imag[jj * dims[0] + ii];
      }
    }
    return 0;
  }

#else

  template <int i> int put(const CDIP &arg) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    std::complex<double> *return_complex = std::get<0>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
    std::complex<double> *cur_pointer = reinterpret_cast<std::complex<double> *>(mxGetComplexDoubles(m));

    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = return_complex[jj * dims[0] + ii];
      }
    }
    return 0;
  }

#endif

  // Input is complex eigen array stored as interleaved real / complex
  // How it is returned to MATLAB depends on whether COMPLEX_SPLIT is set
  // Note: This has not yet been tested with Matlab 2018a (basic idea should work but may need to be tweaked)
  template <int i, typename Derived, typename = std::enable_if_t<std::is_same<typename Derived::Scalar, std::complex<double>>::value>>
  int put(const Eigen::MatrixBase<Derived> &arg, std::complex<double> *ignored = nullptr) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = arg.rows();
    dims[1] = arg.cols();

    mxArray *m = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);

#ifdef COMPLEX_SPLIT

    double *cur_pointer_real = static_cast<double *>(mxGetPr(m));
    double *cur_pointer_imag = static_cast<double *>(mxGetPi(m));

    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {

        cur_pointer_real[jj * dims[0] + ii] = arg(ii, jj).real();
        cur_pointer_imag[jj * dims[0] + ii] = arg(ii, jj).imag();
      }
    }

#else

    std::complex<double> *cur_pointer = reinterpret_cast<std::complex<double> *>(mxGetComplexDoubles(m));

    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = arg(ii, jj);
      }
    }

#endif

    plhs[i] = m;

    return 0;
  }

  template <int i, typename Derived, typename = std::enable_if_t<std::is_same<typename Derived::Scalar, std::complex<float>>::value>>
  int put(const Eigen::MatrixBase<Derived> &arg, std::complex<float> *ignored = nullptr) {

    mwSize dims[2];
    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    dims[0] = arg.rows();
    dims[1] = arg.cols();

    mxArray *m = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxCOMPLEX);

#ifdef COMPLEX_SPLIT

    float *cur_pointer_real = reinterpret_cast<float *>(mxGetPr(m));
    float *cur_pointer_imag = reinterpret_cast<float *>(mxGetPi(m));

    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer_real[jj * dims[0] + ii] = arg(ii, jj).real();
        cur_pointer_imag[jj * dims[0] + ii] = arg(ii, jj).imag();
      }
    }

#else

    std::complex<float> *cur_pointer = reinterpret_cast<std::complex<float> *>(mxGetComplexSingles(m));

    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = arg(ii, jj);
      }
    }

#endif

    plhs[i] = m;

    return 0;
  }

#ifdef COMPLEX_SPLIT

  template <int i> int put(const CFSP &arg) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    float *return_real = std::get<0>(arg).first;
    float *return_imag = std::get<0>(arg).second;

    mxArray *m = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxCOMPLEX);
    float *cur_pointer_real = reinterpret_cast<float *>(mxGetPr(m));
    float *cur_pointer_imag = reinterpret_cast<float *>(mxGetPi(m));

    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer_real[jj * dims[0] + ii] = return_real[jj * dims[0] + ii];
        cur_pointer_imag[jj * dims[0] + ii] = return_imag[jj * dims[0] + ii];
      }
    }
    return 0;
  }

#else

  template <int i> int put(const CFIP &arg) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    std::complex<float> *return_complex = std::get<0>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxCOMPLEX);
    std::complex<float> *cur_pointer = reinterpret_cast<std::complex<float> *>(mxGetComplexSingles(m));

    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = return_complex[jj * dims[0] + ii];
      }
    }
    return 0;
  }

#endif

  template <int i, typename S> std::enable_if_t<std::is_same<S, std::complex<double>>::value, int> put(const S &arg) {
    //    int put(const std::complex<double> &arg) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = 1;
    dims[1] = 1;
    mxArray *m = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);

#ifdef COMPLEX_SPLIT

    double *cur_pointer_real = static_cast<double *>(mxGetPr(m));
    double *cur_pointer_imag = static_cast<double *>(mxGetPi(m));

    cur_pointer_real[0] = arg.real();
    cur_pointer_imag[0] = arg.imag();

#else

    std::complex<double> *cur_pointer = reinterpret_cast<std::complex<double> *>(mxGetComplexDoubles(m));
    cur_pointer[0] = arg;

#endif

    plhs[i] = m;
    return 0;
  }

  template <int i, typename S> std::enable_if_t<std::is_same<S, std::complex<float>>::value, int> put(const S &arg) {
    //  template <int i> int put(const std::complex<float> &arg) {

    //    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    mwSize dims[2];
    dims[0] = 1;
    dims[1] = 1;
    mxArray *m = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxCOMPLEX);

#ifdef COMPLEX_SPLIT

    float *cur_pointer_real = reinterpret_cast<float *>(mxGetPr(m));
    float *cur_pointer_imag = reinterpret_cast<float *>(mxGetPi(m));

    cur_pointer_real[0] = arg.real();
    cur_pointer_imag[0] = arg.imag();

#else

    std::complex<float> *cur_pointer = reinterpret_cast<std::complex<float> *>(mxGetComplexSingles(m));
    cur_pointer[0] = arg;

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

  template <int i, class S = void>
  std::enable_if_t<!in_type_list<S, EDRM, EDCIM, EDCSM, EDR, RDP, CDIP, CDSP, EFRM, EFCIM, EFCSM, EFR, RFP, CFIP, CFSP, EDC, EFC, std::complex<float>, std::complex<double>>() &&
                       !std::is_scalar<S>::value && !is_tuple<S>::value,
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
  template <int i, class S> int put(const std::vector<S> &arg) {
    recursePackVec(std::make_index_sequence<boost::pfr::tuple_size<S>::value>(), i, arg);
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
