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

    uint32_t x[5] = {1,2,3,4,5};

    std::unique_ptr<uint32_t[]> y(new uint32_t[5]);
    for(uint32_t i =0; i<5;i++)
        y[i]=i*i;

    std::shared_ptr<uint32_t[]> v(new uint32_t[5]);
    for(int i =0; i<5;i++)
        v[i]=i*i*i;

    shared_ptr_tuple<uint32_t> u(v,5,1);
    unique_ptr_tuple<uint32_t> z(std::move(y),5,1);
    userStruct w = {1, 3.0, std::move(z) };

    std::vector<uint16_t> t = {10,11,12,13,14};

    MexPacker<uint32_t[5],shared_ptr_tuple<uint32_t>,userStruct,std::vector<uint16_t>> my_pack(nlhs, plhs); //Create a packing object
    my_pack.PackMex(x,u,w,t); 

//    MexPacker<uint32_t[5],shared_ptr_tuple<uint32_t>,userStruct> my_pack(nlhs, plhs); //Create a packing object
//    my_pack.PackMex(x,u,w); 

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
