#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <typeindex>

using namespace MexPackUnpackTypes;
struct userStruct{
  int a;
  double b;
  EDRM c;
  std::string d;
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


  MexUnpacker<userStruct> my_unpack(nrhs, prhs);

  try {

    auto [d] = my_unpack.unpackMex();
    d.a=4;

    MexPacker<userStruct >my_pack(nlhs, plhs); 
    my_pack.PackMex(d); 

  } catch (std::string s) {
    mexPrintf(s.data());
  }
}