#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <iostream>

using namespace MexPackUnpackTypes;
struct userStruct{
  int a;
  double b;
  EDRM c;
  std::string d;
};


struct userStruct2{
  int a;
  double b;
  std::string c; 
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


  //To indicate that a struct type will be used to receive a matlab struct it needs to be "wrapped"
  //in MXStruct
  MexUnpacker<MXStruct<userStruct>> my_unpack(nrhs, prhs);

  try {

    auto [d] = my_unpack.unpackMex();
    d.a=4;

    MexPacker<MXStruct<userStruct> >my_pack(nlhs, plhs); 
    my_pack.PackMex(MXStruct(d)); //MXStruct just stores a reference to the actual struct d to be returned

  } catch (std::string s) {
    mexPrintf(s.data());
  }
}