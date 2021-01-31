#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>

using namespace MexPackUnpackTypes;

struct userStruct{
  int a;
  double b;
  unique_ptr_tuple<uint32_t> c;
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  try {

//    MexUnpacker<ptr_tuple_3dim<double>  > my_unpack(nrhs, prhs);//We need the pointers to the inputs from matlab/octave
//    auto [a] = my_unpack.unpackMex();

    MexUnpacker<ptr_tuple_3dim<double>,ptr_tuple_3dim_CI<double> > my_unpack(nrhs, prhs);//We need the pointers to the inputs from matlab/octave
    auto [a,b] = my_unpack.unpackMex();


    MexPacker<ptr_tuple_3dim<double> ,ptr_tuple_3dim_CI<double> > my_pack(nlhs, plhs); //Create a packing object
    my_pack.PackMex(a,b); 

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
