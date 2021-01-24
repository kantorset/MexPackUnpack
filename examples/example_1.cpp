#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>

using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // We expect two real double matrices (treated as eigen maps), a scalar double, and a string
  MexUnpacker<EDRM,  EDRM,double,std::string> my_unpack(nrhs, prhs);

  try {

    auto [a,b,c,d] = my_unpack.unpackMex(); //a,b are Eigen maps, c is double, d is std::string
    Eigen::MatrixXd e = a+b*c;
    d+=std::string(" appended to the string");
    MexPacker<Eigen::MatrixXd,std::string> my_pack(nlhs, plhs); 
    my_pack.PackMex(e,d);    // We return an eigen Matrix and a string

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}

