#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <iostream>


using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


  //argument 1 : double
  //argument 2 : Single precision complex
  //argument 3 : Double precision complex matrix represented as pair of real matrices corresponding to real/imaginary component 
  //argument 4 : Double precision complex matrix represented as a pair of pointers, number of rows, and number of columns
  //Note: on MATLAB 2018b an newer (with the -R2018a compilation flag) you would replace
  //CDSP with CDIP and EDCSM with EDCIM
  MexUnpacker<double, std::complex<float>, EDCSM,  CDSP> my_unpack(nrhs, prhs);

  
  try {

    auto [a, b, c, d] = my_unpack.unpackMex();
    auto [cr, ci] = c; //cr and ci and Eigen::Map<MatrixXd> corresponding to real and imaginary parts of matrix
    auto [d_p, d_M, d_N] = d; //dp is std::pair<double*,double*> (pointers to real and imginary part), d_M is number of rows, d_N is number of columns

    Eigen::MatrixXcd c_comp(cr.rows(), cr.cols());
    c_comp.real() = cr;
    c_comp.imag() = ci;
    c_comp *= a;
    b*=a;
    MexPacker<std::complex<double>, int, Eigen::MatrixXcd, CDSP> my_pack(nlhs, plhs);
    my_pack.PackMex(b, 2, c_comp, d);

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
