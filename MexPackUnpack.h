#ifndef __MexMatrix
#define __MexMatrix

#include "boost/pfr.hpp"
#include "mex.h"
#include <Eigen>
#include <string>
#include <tuple>
#include <utility>

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

namespace MexPackUnpackTypes {

using EDRM = Eigen::Map<Eigen::MatrixXd>;                                          // Eigen Double Real Map
using EDCIM = Eigen::Map<Eigen::MatrixXcd>;                                        // Eigen double complex interleaved map
using EDCSM = std::pair<Eigen::Map<Eigen::MatrixXd>, Eigen::Map<Eigen::MatrixXd>>; // Eigen double complex split map
using EDR = Eigen::MatrixXd;                                                       // Eigen Double Real Matrix
using RDP = std::tuple<double *, std::size_t, std::size_t>;                        // Pointer to double of underlying data (real) with dimensions (row,col)
using CDIP = std::tuple<std::complex<double> *, std::size_t, std::size_t>;         // Pointer to complex<double> of underlying data (interleaved complex) with dimensions (row,col)
using CDSP = std::tuple<std::pair<double *, double *>, std::size_t, std::size_t>;  // Pair of pointers to doubles of underlying data (split complex) with dimensions (row,col)

using EFRM = Eigen::Map<Eigen::MatrixXf>;                                          // Eigen Single prec Real Map
using EFCIM = Eigen::Map<Eigen::MatrixXcf>;                                        // Eigen Single prec complex interleaved map
using EFCSM = std::pair<Eigen::Map<Eigen::MatrixXf>, Eigen::Map<Eigen::MatrixXf>>; // Eigen Single prec complex split map
using EFR = Eigen::MatrixXf;                                                       // Eigen Single Prec Real Matrix
using RFP = std::tuple<float *, std::size_t, std::size_t>;                         // Pointer to float of underlying data (real) with dimensions (row,col)
using CFIP = std::tuple<std::complex<float> *, std::size_t, std::size_t>;          // Pointer to complex<double> of underlying data (interleaved complex) with dimensions (row,col)
using CFSP = std::tuple<std::pair<float *, float *>, std::size_t, std::size_t>;    // Pair of pointers to floats of underlying data (interleaved complex) with dimensions (row,col)

// Wrapper for dealing with MATLAB/Octave Structs
template <class T> struct MXStruct {
  const T &obj;
  MXStruct(T &in) : obj{in} {}
};

// Wrapper for dealing with MATLAB/Octave Struct array
template <class T> struct MXVecStruct {
  const std::vector<T> &obj;
  MXVecStruct(std::vector<T> &in) : obj{in} {}
};

} // namespace MexPackUnpackTypes

using namespace MexPackUnpackTypes;

template <class T, class S = void> struct TypeMap { using type = T; };

template <class S> struct TypeMap<MXStruct<S>> { using type = S; };

template <class S> struct TypeMap<MXVecStruct<S>> { using type = std::vector<S>; };

// From https://stackoverflow.com/questions/38561242/struct-to-from-stdtuple-conversion
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

template <typename... T> class MexUnpacker;

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

  double get_(int i, double *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxDOUBLE_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not an double\n");
    return static_cast<double>(mxGetScalar(prhs[i]));
  }

  int get_(int i, int *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxINT32_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not an int32\n");
    return static_cast<int>(mxGetScalar(prhs[i]));
  }

  float get_(int i, float *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxSINGLE_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a single\n");
    return static_cast<float>(mxGetScalar(prhs[i]));
  }

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

    return static_cast<std::complex<double> *>(mxGetComplexDoubles(prhs[i]))[0];
  }

  EDCIM get_(int i, EDCIM *ignored) {

    if (!(mxGetClassID(prhs[i]) == mxDOUBLE_CLASS)) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a double array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    Eigen::Map<Eigen::MatrixXcd> new_map(static_cast<std::complex<double> *>(mxGetComplexDoubles(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    return new_map;
  }

  CDIP get_(int i, CDIP *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxDOUBLE_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a double array\n");
    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }
    std::complex<double> *complex_pointer = static_cast<std::complex<double> *>(mxGetComplexDoubles);
    return {complex_pointer, mxGetM(prhs[i]), mxGetN(prhs[i])};
  }

#endif

  RDP get_(int i, RDP *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxDOUBLE_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a double array\n");

    return {static_cast<double *>(mxGetPr(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i])};
  }

  // Same as above but for single (32 bit float) cases

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

    return static_cast<std::complex<float> *>(mxGetComplexSingles(prhs[i]))[0];
  }

  EFCIM get_(int i, EFCIM *ignored) {

    if (!(mxGetClassID(prhs[i]) == mxSINGLE_CLASS)) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a single array\n");
    }

    if (!mxIsComplex(prhs[i])) {
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a complex array\n");
    }

    Eigen::Map<Eigen::MatrixXcf> new_map(static_cast<std::complex<float> *>(mxGetComplexSingles(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i]));
    return new_map;
  }

#endif

  RFP get_(int i, RFP *ignored) {
    if (!(mxGetClassID(prhs[i]) == mxSINGLE_CLASS))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a single array\n");
    return {reinterpret_cast<float *>(mxGetPr(prhs[i])), mxGetM(prhs[i]), mxGetN(prhs[i])};
  }

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

   //Recursively use unpacker to handle structs
  template <class S, std::size_t... Is> std::tuple<boost::pfr::tuple_element_t<Is, S>...> recurseUnpack(std::index_sequence<Is...>, int i, S *ignored) {
    int num_fields = mxGetNumberOfFields(prhs[i]);
    std::unique_ptr<mxArray *[]> sub_rhs(new mxArray *[num_fields]);
    for (int field_ind = 0; field_ind < num_fields; field_ind++) {
      sub_rhs[field_ind] = mxGetFieldByNumber(prhs[i], 0, field_ind);
    }
    MexUnpacker<boost::pfr::tuple_element_t<Is, S>...> sub_unpack(num_fields, const_cast<const mxArray **>(sub_rhs.get()));
    return sub_unpack.unpackMex();
  }

  template <class S> S get_(int i, MXStruct<S> *ignored) {
    if (!mxIsStruct(prhs[i]))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a struct\n");
    S *S_ignored = nullptr;
    return make_struct<S>(recurseUnpack(std::make_index_sequence<boost::pfr::tuple_size<S>::value>(), i, S_ignored));
  }



   //Recursively use unpacker to handle struct arrays
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

  template <class S> std::vector<S> get_(int i, MXVecStruct<S> *ignored) {
    if (!mxIsStruct(prhs[i]))
      throw std::string("Argument ") + std::to_string(i) + std::string(" not a struct\n");
    S *S_ignored = nullptr;
    return recurseUnpackVec(std::make_index_sequence<boost::pfr::tuple_size<S>::value>(), i, S_ignored);
  }

  // Calls get_ with index of each type and
  // (Null) pointer to appropriate type to select the correect overload (probably more elegant way to do this)
  template <int i> typename NthElement<i, typename TypeMap<T>::type...>::type get() {
    if (i >= nrhs)
      throw std::string("Can't extract argument ") + std::to_string(i) + std::string(", only ") + std::to_string(nrhs) + std::string(" arguments on right side.\n");
    typename NthElement<i, T...>::type *ignored = nullptr;
    return get_(i, ignored);
  }

  // Calls get for each index in type list parameter pack
  template <std::size_t... Is> std::tuple<typename TypeMap<T>::type...> unpackIndirect(std::index_sequence<Is...>) { return {get<Is>()...}; }

  // User function to unpack all the inputs from Matlab/Octave
  // Need to get index of each type in the pack which is dones by calling unpackIndirect
  std::tuple<typename TypeMap<T>::type...> unpackMex() { return unpackIndirect(std::make_index_sequence<sizeof...(T)>()); }
};

// Needed to evaluate arguments in a parameter pack
template <class... T> void ignore(T &&...) {}

// Class to return data from c++ mex function back to matlab / octave
/*
Usage:
  MexPacker<double,int,Eigen::MatrixXcd,Eigen::MatrixXf,RDP> my_pack(nlhs,plhs);
  my_pack.PackMex(3.8,2,d,f,g);
*/

template <typename... T> class MexPacker;

template <typename... T> class MexPacker {
public:
  int nlhs;
  mxArray **plhs;

  MexPacker(int nlhs_, mxArray *plhs_[]) : nlhs{nlhs_}, plhs{plhs_} {}

#ifdef COMPLEX_SPLIT

  template <int i> int put(const CDSP &arg) {

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
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

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    std::complex<double> *return_complex = std::get<0>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
    std::complex<double> *cur_pointer = reinterpret_cast<std::complex<double> *>(mxGetComplexDoubles(m));

    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer_real[jj * dims[0] + ii] = return_complex[jj * dims[0] + ii];
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

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
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

  //Recursively call Pack to handle structs
  template <class S, std::size_t... Is> void recursePack(std::index_sequence<Is...>, int i, S &arg) {
    int num_fields = boost::pfr::tuple_size<S>::value;
    std::vector<std::string> field_names;
    std::vector<const char *> strings;
    for (int j = 0; j < num_fields; ++j)
      field_names.push_back("field_" + std::to_string(j));

    for (int j = 0; j < num_fields; ++j)
      strings.push_back(field_names[j].c_str());

    plhs[i] = mxCreateStructMatrix(1, 1, num_fields, strings.data());

    std::unique_ptr<mxArray *[]> sub_lhs(new mxArray *[num_fields]);

    MexPacker<boost::pfr::tuple_element_t<Is, S>...> sub_pack(num_fields, sub_lhs.get());
    sub_pack.PackMex(boost::pfr::get<Is, S>(arg)...);
    for (int field_ind = 0; field_ind < num_fields; field_ind++) {
      mxSetFieldByNumber(plhs[i], 0, field_ind, sub_lhs[field_ind]);
    }
  }

  template <int i, class S = void> int put(const MXStruct<S> &arg) {
    recursePack(std::make_index_sequence<boost::pfr::tuple_size<S>::value>(), i, arg.obj);
    return 0;
  }


  //Recursively call Pack to handle struct arrays
  template <class S, std::size_t... Is> void recursePackVec(std::index_sequence<Is...>, int i, const std::vector<S> &arg) {
    int num_fields = boost::pfr::tuple_size<S>::value;
    std::vector<std::string> field_names;
    std::vector<const char *> strings;

    for (int j = 0; j < num_fields; ++j)
      field_names.push_back("field_" + std::to_string(j));

    for (int j = 0; j < num_fields; ++j)
      strings.push_back(field_names[j].c_str());

    plhs[i] = mxCreateStructMatrix(arg.size(), 1, num_fields, strings.data());

    std::unique_ptr<mxArray *[]> sub_lhs(new mxArray *[num_fields]);
    //     S tmp;
    for (std::size_t row_ind = 0; row_ind < arg.size(); row_ind++) {

      MexPacker<boost::pfr::tuple_element_t<Is, S>...> sub_pack(num_fields, sub_lhs.get());
      sub_pack.PackMex(boost::pfr::get<Is, S>(arg[row_ind])...);
      for (int field_ind = 0; field_ind < num_fields; field_ind++) {
        mxSetFieldByNumber(plhs[i], row_ind, field_ind, sub_lhs[field_ind]);
      }
    }
  }

  template <int i, class S = void> int put(const MXVecStruct<S> &arg) {
    recursePackVec(std::make_index_sequence<boost::pfr::tuple_size<S>::value>(), i, arg.obj);

    return 0;
  }

  template <int i, typename Derived, typename = std::enable_if_t<std::is_same<typename Derived::Scalar, double>::value>> int put(const Eigen::MatrixBase<Derived> &arg, double *ignored = nullptr) {

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
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

  template <int i> int put(const RDP &arg) {

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double *cur_pointer = static_cast<double *>(mxGetPr(m));
    double *return_pointer = std::get<0>(arg);

    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = return_pointer[jj * dims[0] + ii];
      }
    }
    return 0;
  }

  // Versions for float arrays

  template <int i, typename Derived, typename = std::enable_if_t<std::is_same<typename Derived::Scalar, std::complex<float>>::value>>
  int put(const Eigen::MatrixBase<Derived> &arg, std::complex<float> *ignored = nullptr) {

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
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

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
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

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
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

  template <int i> int put(const RFP &arg) {

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    dims[0] = std::get<1>(arg);
    dims[1] = std::get<2>(arg);

    mxArray *m = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    float *cur_pointer = reinterpret_cast<float *>(mxGetPr(m));
    float *return_pointer = std::get<0>(arg);

    plhs[i] = m;
    for (std::size_t jj = 0; jj < dims[1]; jj++) {
      for (std::size_t ii = 0; ii < dims[0]; ii++) {
        cur_pointer[jj * dims[0] + ii] = return_pointer[jj * dims[0] + ii];
      }
    }
    return 0;
  }

  template <int i, typename Derived, typename = std::enable_if_t<std::is_same<typename Derived::Scalar, float>::value>> int put(const Eigen::MatrixBase<Derived> &arg, float *x = nullptr) {

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
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

  template <int i> int put(const float &arg) {

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    dims[0] = 1;
    dims[1] = 1;

    mxArray *m = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    float *cur_pointer = reinterpret_cast<float *>(mxGetPr(m));
    cur_pointer[0] = arg;
    plhs[i] = m;
    return 0;
  }

  template <int i> int put(const double &arg) {

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    dims[0] = 1;
    dims[1] = 1;

    mxArray *m = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double *cur_pointer = static_cast<double *>(mxGetPr(m));
    cur_pointer[0] = arg;
    plhs[i] = m;
    return 0;
  }

  template <int i> int put(const std::complex<double> &arg) {

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
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

  template <int i> int put(const std::complex<float> &arg) {

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    dims[0] = 1;
    dims[1] = 1;
    mxArray *m = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxCOMPLEX);

#ifdef COMPLEX_SPLIT

    float *cur_pointer_real = reinterpret_cast<float *>(mxGetPr(m));
    float *cur_pointer_imag = reinterpret_cast<float *>(mxGetPi(m));

    cur_pointer_real[0] = arg.real();
    cur_pointer_imag[0] = arg.imag();

#else

    std::complex<float> *cur_pointer = static_cast<std::complex<float> *>(mxGetComplexSingles(m));
    cur_pointer[0] = arg;

#endif

    plhs[i] = m;
    return 0;
  }

  template <int i> int put(const std::string &arg) {
    plhs[i] = mxCreateString(arg.c_str());
    return 0;
  }

  template <int i> int put(const int &arg) {

    mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
    dims[0] = 1;
    dims[1] = 1;

    mxArray *m = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    int *cur_pointer = reinterpret_cast<int *>(mxGetPr(m));
    cur_pointer[0] = arg;
    plhs[i] = m;
    return 0;
  }

  template <int i, typename S> int check_and_put(const S &arg) {
    if (i >= nlhs) {
      std::string s = std::string("Can't return argument ") + std::to_string(i) + std::string(", only ") + std::to_string(nlhs) + std::string(" arguments on right side.\n");
      mexPrintf(s.data());
    } else {
      put<i>(arg);
    }
  }

  // ignore does nothing but lets us expand the parameter pack to evaluate all the put functions for each type
  template <std::size_t... Is> void packIndirect(const T &... args, std::index_sequence<Is...>) {
    ignore(check_and_put<Is, T>(args)...);
  }

  // Need to get index of each type on the pack
  void PackMex(const T &... args) { return packIndirect(args..., std::make_index_sequence<sizeof...(T)>()); }
};

#endif
