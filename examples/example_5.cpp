#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>

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
  MXVecStruct<userStruct> d; //To return a nested struct it needs to be wrapped
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  MexUnpacker<MXVecStruct<userStruct>> my_unpack(nrhs, prhs);

  try {

    auto [d] = my_unpack.unpackMex();

    MexPacker<MXStruct<userStruct2> >my_pack(nlhs, plhs);

    auto d1 = MXVecStruct(d);
    userStruct2 output_struct{1,4.5,"doc",d1};

    my_pack.PackMex(MXStruct(output_struct));

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}