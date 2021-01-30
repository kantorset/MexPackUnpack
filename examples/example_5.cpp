#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <typeindex>
#include <map>

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
  std::vector<userStruct> d; 
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  MexUnpacker<std::vector<userStruct>> my_unpack(nrhs, prhs);

  try {

    auto [d] = my_unpack.unpackMex();


    std::map<std::type_index,std::vector<std::string> > my_name_map;

    my_name_map[typeid(userStruct)] = {"my_int","my_double","my_array", "my_string"};
    my_name_map[typeid(userStruct2)] = {"my_int","my_double", "my_string","my_vector"};

    MexPacker<userStruct2>my_pack(nlhs, plhs,my_name_map);

    userStruct2 output_struct{1,4.5,"doc",d};

    my_pack.PackMex(output_struct);

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}