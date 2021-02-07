#ifndef __MexPackUnpackTypes
#define __MexPackUnpackTypes

#include "mex.h"
#include <Eigen>
#include <string>
#include <tuple>
#include <utility>
#include <memory>

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
using EDC = Eigen::MatrixXcd;
using EFC = Eigen::MatrixXcf;

template<typename T>
using ptr_tuple = std::tuple<T*,std::size_t,std::size_t>;

template<typename T>
using unique_ptr_tuple = std::tuple<std::unique_ptr<T[]>,std::size_t,std::size_t>;

template<typename T>
using shared_ptr_tuple = std::tuple<std::shared_ptr<T[]>,std::size_t,std::size_t>;


template<typename T>
using ptr_tuple_3dim = std::tuple<T*,std::size_t,std::size_t,std::size_t>;

//For 3 dimensional arrays with split complex representation
template<typename T>
using ptr_tuple_3dim_CS = std::tuple<std::pair<T*,T*>,std::size_t,std::size_t,std::size_t>;

//Future?
//template<typename T>
//using unique_ptr_tuple_3dim = std::tuple<std::unique_ptr<T[]>,std::size_t,std::size_t,std::size_t>;

//template<typename T>
//using shared_ptr_tuple_3dim = std::tuple<std::shared_ptr<T[]>,std::size_t,std::size_t,std::size_t>;



} // namespace MexPackUnpackTypes


/*
template <class T> struct EigenType;

template <> struct EigenType<double>{
  using EigenMap = Eigen::Map<Eigen::MatrixXd>;
};
  

template <> struct EigenType<float>{
  using EigenMap = Eigen::Map<Eigen::MatrixXf>;
};
*/


template <class T> struct isEigenMap{
  constexpr static bool value=false;
};

template <> struct isEigenMap<Eigen::Map<Eigen::MatrixXd> > {
  constexpr static bool value=true; 
};

template <> struct isEigenMap<Eigen::Map<Eigen::MatrixXf> > {
  constexpr static bool value=true; 
};

template <> struct isEigenMap<Eigen::Map<Eigen::MatrixXcd> > {
  constexpr static bool value=true; 
};

template <> struct isEigenMap<Eigen::Map<Eigen::MatrixXcf> > {
  constexpr static bool value=true; 
};

  
//  using EigenDyn = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;

template <class T> struct mxClassTraits;

template <> struct mxClassTraits<double>
{ 
  constexpr static mxClassID mxClass = mxDOUBLE_CLASS;
  constexpr static const char* name = "double";
};

template <> struct mxClassTraits<int>
{ 
  constexpr static mxClassID mxClass = mxINT32_CLASS;
  constexpr static const char* name = "int";
};


template <> struct mxClassTraits<uint32_t>
{ 
  constexpr static mxClassID mxClass = mxUINT32_CLASS;
  constexpr static const char* name = "uint32";
};


template <> struct mxClassTraits<int16_t>
{ 
  constexpr static mxClassID mxClass = mxINT16_CLASS;
  constexpr static const char* name = "int16";
};


template <> struct mxClassTraits<uint16_t>
{ 
  constexpr static mxClassID mxClass = mxUINT16_CLASS;
  constexpr static const char* name = "uint16";
};

template <> struct mxClassTraits<int8_t>
{ 
  constexpr static mxClassID mxClass = mxINT8_CLASS;
  constexpr static const char* name = "int8";
};


template <> struct mxClassTraits<uint8_t>
{ 
  constexpr static mxClassID mxClass = mxUINT8_CLASS;
  constexpr static const char* name = "int8";
};


template <> struct mxClassTraits<float>
{ 
  constexpr static mxClassID mxClass = mxSINGLE_CLASS;
  constexpr static const char* name = "float";
};


#endif
