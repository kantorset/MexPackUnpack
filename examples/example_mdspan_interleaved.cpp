#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <iostream>

using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Create an unpacker with template parameters corresponding to what we expect from matlab
  // argument 1: Scalar double
  // argument 2: EDRM (Eigen map which is basically an Eigen matrix wrapping an existing pointer)

  MexUnpacker<
    span_2d_dynamic_left<double>,
    span_2d_dynamic_left<std::complex<double>
    >
  > my_unpack(nrhs, prhs);//We need the pointers to the inputs from matlab/octave

  try {

    auto [a,b] = my_unpack.unpackMex(); //  a is double, b is Eigen::Map<Eigen::MatrixXd> (which EDRM aliases)
                                         //  Note that in general one must be careful using auto with eigen 
                                         //  objects due to its use of expression templates. However, in 
                                         //  this case auto correctly deduces returned eigen objects types.

  auto c = stdex::submdspan(b,std::pair{2ul,4ul},std::pair{3ul,5ul});
// unpackMex can throw if the matlab objects passed at run time
// are not compatible with what is expected
    MexPacker<  decltype(c) > my_pack(nlhs, plhs); //Create a packing object
//    MexPacker<  span_2d_dynamic_left<std::complex<double> > > my_pack(nlhs, plhs); //Create a packing object
//    a(0,0)=10;
    my_pack.PackMex(c); //Return this object to matlab

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
