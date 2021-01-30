#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <variant>
#include <typeindex>
#include <map>

using namespace MexPackUnpackTypes;
struct userStruct1{
  int a;
  double b;
};


struct userStruct2{
  int a;
  std::string c; 
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  try {

    // This example takes no inputs
    userStruct1 output_struct1{1,4.5};
    userStruct2 output_struct2{1,"doc"};

    std::map<std::type_index,std::vector<std::string> > my_name_map;
    my_name_map[typeid(userStruct1)] = {"my_int","my_double"};
    my_name_map[typeid(userStruct2)] = {"my_int", "my_string"};


    std::vector<std::variant<userStruct1,userStruct2,double>>variant_vec;
    variant_vec.push_back(output_struct1);
    variant_vec.push_back(output_struct1);
    variant_vec.push_back(output_struct2);  
    variant_vec.push_back(output_struct2);
    variant_vec.push_back(3.0);

    MexPacker<std::vector<std::variant<userStruct1,userStruct2,double>>>my_pack(nlhs, plhs,my_name_map);
    my_pack.PackMex(variant_vec);

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}