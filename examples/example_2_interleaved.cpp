#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>

using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


  //argument 1 : real double precision scalar
  //argument 2 : single precision complex scalar
  //argument 3 : Double precision complex matrix (Eigen Map) 
  //argument 4 : Double precision complex matrix represented as a pair of pointers, number of rows, and number of columns
  //Note: Requires MATLAB 2018a or newer and  the -R2018a compilation flag
  MexUnpacker<double, std::complex<float>, EDCIM,  CDIP> my_unpack(nrhs, prhs);
 
  try {

    auto [a, b, c, d] = my_unpack.unpackMex();
    auto [d_p, d_M, d_N] = d; //dp is std::complex<double>* (pointer to complex data), d_M is number of rows, d_N is number of columns
    
    Eigen::MatrixXcd e_comp = b*a;

    MexPacker<std::complex<double>, int, Eigen::MatrixXcd, CDIP> my_pack(nlhs, plhs);
    my_pack.PackMex(b, 2, e_comp, d);

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
