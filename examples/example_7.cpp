#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>

using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


  MexUnpacker<uint32_t,  ptr_tuple<uint32_t>,ptr_tuple<uint8_t>,ptr_tuple<int16_t>> my_unpack(nrhs, prhs);//We need the pointers to the inputs from matlab/octave

  try {

    auto [a,b,c,d] = my_unpack.unpackMex(); 

    a = 16;

    MexPacker<uint32_t,  ptr_tuple<uint32_t>,ptr_tuple<uint8_t>,ptr_tuple<int16_t> > my_pack(nlhs, plhs); //Create a packing object
    my_pack.PackMex(a,b,c,d); //Return these objects to matlab

// unpackMex can throw if the matlab objects passed at run time
// are not compatible with what is expected
  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
