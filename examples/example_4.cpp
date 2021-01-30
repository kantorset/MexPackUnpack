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


struct userStruct2{
  int a;
  double b;
  std::string c; 
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  //Indicate we will receive a struct array we want to unpack into std::vector<userStruct>
  MexUnpacker<std::vector<userStruct>> my_unpack(nrhs, prhs);

  try {

    auto [d] = my_unpack.unpackMex(); //d is std::vector<userStruct>

    userStruct2 output_struct{d[1].a,d[1].b,"doc"};
    std::vector<userStruct2> output_vec;
    output_vec.push_back(output_struct);
    output_vec.push_back(output_struct);

    std::map<std::type_index,std::vector<std::string> > my_name_map;
    my_name_map[typeid(userStruct2)] = {"my_int","my_double", "my_string"};

    MexPacker<std::vector<userStruct2> >my_pack(nlhs, plhs,my_name_map);
    my_pack.PackMex(output_vec); 

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}