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
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  //Indicate we will receive a struct array we want to unpack into std::vector<userStruct>
  MexUnpacker<MXVecStruct<userStruct>> my_unpack(nrhs, prhs);

  try {

    auto [d] = my_unpack.unpackMex(); //d is std::vector<userStruct>

    userStruct2 output_struct{d[1].a,d[1].b,"doc"};
    std::vector<userStruct2> output_vec;
    output_vec.push_back(output_struct);
    output_vec.push_back(output_struct);

    MexPacker<MXVecStruct<userStruct2> >my_pack(nlhs, plhs);
    my_pack.PackMex(MXVecStruct(output_vec)); //MXVecStruct just stores a reference to the vector of structs we are returning

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}