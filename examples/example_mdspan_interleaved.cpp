#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <iostream>

using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Create an unpacker with template parameters corresponding to what we expect from matlab
  // argument 1: 4 dimensional mdspan

  MexUnpacker<span_4d_dynamic_left<std::complex<double>>  > my_unpack(nrhs, prhs); //We need the pointers to the inputs from matlab/octave

  try {

    auto [a] = my_unpack.unpackMex(); 
    if(a.extent(1)<4||a.extent(2)<3||a.extent(3)<4){
      throw std::string("Input is too small for slice indices");
    }

    std::complex<double> b = a(0,1,2,3);

    auto c = stdex::submdspan(a,1,std::pair{2ul,4ul},stdex::full_extent,3); //Create a slice (submdspan)

    MexPacker<std::complex<double> , decltype(c) > my_pack(nlhs, plhs); //Create a packing object, we use decltype to get the strides from submdspan correct
    my_pack.PackMex(b,c); //Return this object to matlab

  } catch (std::string s) {
    mexPrintf(s.data());
  }
  
}
