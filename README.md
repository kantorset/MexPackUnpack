# MexPackUnpack

MexPackUnpack is a C++ header library to help automate passing MATLAB/Octave data into and out of user C++ mex code. It uses some lite (depending on your point of view) template metaprogramming with c++ variadic templates. It requires C++17 and is intended to be used in conjunction with C++17 structured binding syntax. 

The current version extracts the MATLAB/Octave numeric arrays to either [Eigen](https://eigen.tuxfamily.org/) Map objects (Eigen Matrices wrapping an external  pointer) or a tuple of the underlying pointers and array sizes. MATLAB/Octave strings can be converted to C++ std::strings. MATLAB/Octave structs and struct arrays can be extracted into user defined C++ structs with a matching layout. These types can then be returned and converted to MATLAB/Octave objects.
  

## Compilation

Eigen and boost::pfr are included as git submodules and need to be in the include path when compiling. For example in Octave.

```matlab
>> mkoctfile --mex -v -std=c++17 -I./eigen/Eigen -I./pfr/include/ ./examples/example_0.cpp  
``` 
The library has been tested and compiles successfully with g++ 7.5 or visual studio 2019 for both MATLAB and Octave.

To compile for MATLAB the appropriate command is for linux:
```matlab
>> mex -R2018a -v CXXFLAGS="-std=c++17 -fPIC" -I./eigen/Eigen -I./pfr/include/  -I.  ./examples/example_0.cpp  
```

For windows:
```matlab
>> mex -R2018a -v COMPFLAGS="/std:c++17" -I./eigen/Eigen -I./pfr/include/ -I.  ./examples/example_0.cpp  
```


## Usage
### Simple Numeric Object Unpacking 

#### Example 0
Note that in all cases what is expected to be passed into and out of the  ```MexUnpacker``` and ```MexPacker```  classes is specified by setting the template parameters as shown below. The [type aliases](#type-aliases) used are shown further down. ```unpackMex``` returns a tuple so we use C++17 structured binding to destructure it.
```cpp
#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <iostream>

using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Create an unpacker with template parameters corresponding to what we expect from matlab
  // argument 1: Scalar double
  // argument 2: EDRM (Eigen map which is basically an Eigen matrix wrapping an existing pointer)
  MexUnpacker<double,  EDRM> my_unpack(nrhs, prhs); //We need the pointers to the inputs from matlab/octave

  try {

    auto [a, b] = my_unpack.unpackMex(); //  a is double, b is Eigen::Map<Eigen::MatrixXd> (which EDRM aliases)
                                         //  Note that in general one must be careful using auto with eigen 
                                         //  objects due to its use of expression templates. However, in 
                                         //  this case auto correctly deduces returned eigen objects types.

    Eigen::MatrixXd c=a*b; //Eigen matrices have arithmetic operations overloaded

    MexPacker<Eigen::MatrixXd> my_pack(nlhs, plhs); //Create a packing object
    my_pack.PackMex(c); //Return this object to matlab

// unpackMex can throw if the matlab objects passed at run time
// are not compatible with what is expected
  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
```

MATLAB/Octave code using the resulting mex function:
```matlab
>> m=[[1,2,3];[4,5,6]];
>> n=example_0(3,m);
>> n
n =

    3    6    9
   12   15   18
```




#### Example 1
Here is a slightly more complicated example:
```cpp
#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>

using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // We expect two real double matrices (treated as eigen maps), a scalar double, and a string
  MexUnpacker<EDRM, EDRM, double, std::string> my_unpack(nrhs, prhs);

  try {

    auto [a,b,c,d] = my_unpack.unpackMex(); //a,b are Eigen maps, c is double, d is std::string
    Eigen::MatrixXd e = a+b*c;
    d+=std::string(" appended to the string");

    // We return an eigen Matrix and a string
    MexPacker<Eigen::MatrixXd,std::string> my_pack(nlhs, plhs);
    my_pack.PackMex(e,d);

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
```

```matlab
>> m=[[1,2,3];[4,5,6]];
>> n=[[6,7,8];[9,10,11]];
>> [out1,out2] = example_1(m,n,2.0,'More stuff')
out1 =

   13   16   19
   22   25   28

out2 = More stuff appended to the string
```
#### Error Handling

In principle the code checks that what is received at runtime from MATLAB/Octave is compatible with the types specified in MexUnpacker instantiation. If not, an exception is raised. The intended use is that ```unpackMex``` is called at the beginning of the try block before anything has been done that would require cleanup so that the catch can just print the error and leave things in a good state. 
For the [previous example](#example-1), the last argument should be a string, if not the program just prints a message about the type being wrong to the Octave/MATLAB terminal (note the argument number is 0 based ...)

```matlab
>> [out1,out2] = example_1(m,n,2.0,1.2)
Argument 3 not an string
```
It's likely that all possible error paths have not been completely tested so some fixes might need to be made. Ideally, it should be hard to crash MATLAB/Octave using the library even if the types passed at runtime are totally wrong. 


#### Example 2 (Complex Matrices)

The library handles complex as well as real matrices. There is some complexity due to the fact that matlab 2017b and earlier as well as Octave use a split real/imaginary representation while MATLAB 2018a and later uses an interleaved representation. See discussion [below](#notes-on-complex-numeric-types) about this. This example show use of split complex types. 

```cpp
//Example with complex numbers and separate real / imaginary pointers
#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <iostream>

using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  //argument 1 : real double precision scalar
  //argument 2 : Single precision complex scalar
  //argument 3 : Double precision complex matrix (Eigen Map) represented as pair of real matrices corresponding to real/imaginary component 
  //argument 4 : Double precision complex matrix represented as a pair of pointers, number of rows, and number of columns 
  //This will only work with Octave or matlab before 2017b (or newer MATLABs compiled without the -R2018a flag)
  //Note: on MATLAB 2018a and newer (with the -R2018a compilation flag) you would replace
  //CDSP with CDIP and EDCSM with EDCIM

  MexUnpacker<double, std::complex<float>, EDCSM,  CDSP> my_unpack(nrhs, prhs);

  try {

    auto [a, b, c, d] = my_unpack.unpackMex();
    auto [cr, ci] = c; //cr and ci and Eigen::Map<MatrixXd> corresponding to real and imaginary parts of matrix
    auto [d_p, d_M, d_N] = d; //dp is std::pair<double*,double*> (pointers to real and imaginary part), d_M is number of rows, d_N is number of columns

    Eigen::MatrixXcd c_comp(cr.rows(), cr.cols());
    c_comp.real() = cr;
    c_comp.imag() = ci;
    c_comp *= a;
    MexPacker<std::complex<double>, int, Eigen::MatrixXcd, CDSP> my_pack(nlhs, plhs);
    my_pack.PackMex(b, 2, c_comp, d);

  } catch (std::string s) {
    mexPrintf(s.data());
  }
}
```
Code using the resulting mex function

```matlab
>> m=([[1,2,3];[4,5,6]]+i*[[7,8,9];[10,11,12.0]]);
>> [a,b,c,d] = example_2(3.0,single(2+1.7*i),m,conj(m))
a =  6.0000 + 5.1000i
b = 2
c =

    3 + 21i    6 + 24i    9 + 27i
   12 + 30i   15 + 33i   18 + 36i

d =

    1 -  7i    2 -  8i    3 -  9i
    4 - 10i    5 - 11i    6 - 12i
```

The analogous example using MATLAB 2018a or later with the -R2018a flag which exposes the newer interleaved complex API would be as follows. Note that the Eigen matrices (maps) and pointer tuples to complex matrices are intrinsically complex in this case and not represented as separate real / imaginary components. 
```cpp
//Example with interleaved real / imaginary numbers
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
    
    Eigen::MatrixXcd e_comp = b*c;

    MexPacker<std::complex<double>, int, Eigen::MatrixXcd, CDIP> my_pack(nlhs, plhs);
    my_pack.PackMex(b, 2, e_comp, d);

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
```


### Type Aliases
```cpp

namespace MexPackUnpackTypes {

using EDRM = Eigen::Map<Eigen::MatrixXd>;                                          // Eigen Double Real Map
using EDCIM = Eigen::Map<Eigen::MatrixXcd>;                                        // Eigen double complex interleaved map
using EDCSM = std::pair<Eigen::Map<Eigen::MatrixXd>, Eigen::Map<Eigen::MatrixXd>>; // Eigen double complex split map
using EDR = Eigen::MatrixXd;                                                       // Eigen Double Real Matrix
using RDP = std::tuple<double *, std::size_t, std::size_t>;                        // Pointer to double of underlying data (real) with dimensions (row,col)
using CDIP = std::tuple<std::complex<double> *, std::size_t, std::size_t>;         // Pointer to complex<double> of underlying data (interleaved complex) with dimensions (row,col)
using CDSP = std::tuple<std::pair<double *, double *>, std::size_t, std::size_t>;  // Pair of pointers to doubles of underlying data (split complex) with dimensions (row,col)

using EFRM = Eigen::Map<Eigen::MatrixXf>;                                          // Eigen Single prec Real Map
using EFCIM = Eigen::Map<Eigen::MatrixXcf>;                                        // Eigen Single prec complex interleaved map
using EFCSM = std::pair<Eigen::Map<Eigen::MatrixXf>, Eigen::Map<Eigen::MatrixXf>>; // Eigen Single prec complex split map
using EFR = Eigen::MatrixXf;                                                       // Eigen Single Prec Real Matrix
using RFP = std::tuple<float *, std::size_t, std::size_t>;                         // Pointer to float of underlying data (real) with dimensions (row,col)
using CFIP = std::tuple<std::complex<float> *, std::size_t, std::size_t>;          // Pointer to complex<double> of underlying data (interleaved complex) with dimensions (row,col)
using CFSP = std::tuple<std::pair<float *, float *>, std::size_t, std::size_t>;    // Pair of pointers to floats of underlying data (interleaved complex) with dimensions (row,col)
using EDC = Eigen::MatrixXcd;
using EFC = Eigen::MatrixXcf;

template<typename T> //tuple with pointer and array dimensions, RDP, RFP,CDIP,CFIP are special cases
using ptr_tuple = std::tuple<T*,std::size_t,std::size_t>;

template<typename T>
using unique_ptr_tuple = std::tuple<std::unique_ptr<T[]>,std::size_t,std::size_t>;

template<typename T>
using shared_ptr_tuple = std::tuple<std::shared_ptr<T[]>,std::size_t,std::size_t>;


template<typename T>
using ptr_tuple_3dim = std::tuple<T*,std::size_t,std::size_t,std::size_t>;

//For 2 dimensional arrays with split complex representation, CDSP, and CFSP are special cases
template<typename T>
using ptr_tuple_CS = std::tuple<std::pair<T*,T*>,std::size_t,std::size_t>;
 
//For 3 dimensional arrays with split complex representation
template<typename T>
using ptr_tuple_3dim_CS = std::tuple<std::pair<T*,T*>,std::size_t,std::size_t,std::size_t>;

} // namespace MexPackUnpackTypes
```


### Structs and Struct Arrays

To deal with matlab structs and struct arrays we can  define a user C++ structure and then use that in the type list to automatically extract matlab structs which have fields matching the user specified C++ struct. This is done using the reflection/introspection capability of  [Boost PFR](https://apolukhin.github.io/magic_get/index.html). Any user specified type that is not one of the  [default](#type-aliases) type is assumed to be a struct intended to ingest or return data to/from a MATLAB/Octave struct/struct array.

#### Simple Struct Example
```cpp
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


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


  MexUnpacker<userStruct> my_unpack(nrhs, prhs);

  try {

    auto [d] = my_unpack.unpackMex();
    d.a=4;

    MexPacker<userStruct >my_pack(nlhs, plhs); 
    my_pack.PackMex(d); 

  } catch (std::string s) {
    mexPrintf(s.data());
  }
}
```
Note that while we can determine the types of the user defined fields in the C++ struct using [Boost PFR](https://apolukhin.github.io/magic_get/index.html), we cannot determine the names of the fields from within C++. So currently when structs are returned they are by default just named field_0, field_1, ..., etc. 

```matlab
>> m=struct('a',int32(1),'b',3.7,'c',[1.0,2.0,5.0],'d','whats up');
>> a=example_3(m)
a =

  scalar structure containing the fields:

    field_0 = 4
    field_1 =  3.7000
    field_2 =

       1   2   5

    field_3 = whats up
```

There is an optional third argument to the MexPacker constructor that takes a user defined map between the user struct type (in the form of a ```std::type_index```) and the field names to be used in MATLAB/Octave as a vector of ```strings```. The use of this is illustrated below. 

```cpp
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


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  MexUnpacker<userStruct> my_unpack(nrhs, prhs);

  try {

    auto [d] = my_unpack.unpackMex();
    d.a=4;

    std::map<std::type_index,std::vector<std::string> > my_name_map;
    my_name_map[typeid(userStruct)] = {"my_int","my_double","my_array", "my_string"};

    MexPacker<userStruct >my_pack(nlhs, plhs,my_name_map); 
    my_pack.PackMex(d); 

  } catch (std::string s) {
    mexPrintf(s.data());
  }
}
```
For this version we have the following results.
```matlab
>> m=struct('a',int32(1),'b',3.7,'c',[1.0,2.0,5.0],'d','whats up');
>> a=example_3a(m)
a =

  scalar structure containing the fields:

    my_int = 4
    my_double =  3.7000
    my_array =

       1   2   5

    my_string = whats up
```
Note how the user field names as used to name the struct fields.
#### Struct Array Example
Struct arrays with more than 1 element can be accomodated by specifying a ```std::vector<S>``` where ```S``` is the user defined struct type.
```cpp
#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include<typeindex>

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
```
Example octave code calling it.
```
>> m=struct('a',int32(1),'b',3.7,'c',[1.0,2.0,5.0],'d','whats up');
>> m=[m,m,m,m]
m =

  1x4 struct array containing the fields:

    a
    b
    c
    d

>> a=example_4(m)
a =

  2x1 struct array containing the fields:

    my_int
    my_double
    my_string

>> a(2)
ans =

  scalar structure containing the fields:

    my_int = 1
    my_double =  3.7000
    my_string = doc
```

#### Slightly Ridiculous (Mis)use of Struct/Struct Array Functionality
It's actually possible to have nested structs in C++ returned as nested structs in MATLAB/Octave.  In principle this should work for arbitrary levels of nesting (until the compiler gets upset at least).
```cpp
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
```
Example octave code using this.
```matlab
>> m=struct('a',int32(1),'b',3.7,'c',[1.0,2.0,5.0],'d','whats up');
>> m=[m;m;m;m]
m =

  4x1 struct array containing the fields:

    a
    b
    c
    d

>> a=example_5(m)
a =

  scalar structure containing the fields:

    my_int = 1
    my_double =  4.5000
    my_string = doc
    my_vector =

      4x1 struct array containing the fields:

        my_int
        my_double
        my_array
        my_string


>> a.my_vector(3)
ans =

  scalar structure containing the fields:

    my_int = 1
    my_double =  3.7000
    my_array =

       1   2   5

    my_string = whats up
```


### Cell Arrays
There is very limited ability to pass cell arrays from MATLAB/Octave into C++ currently. At the moment, the only case handled for going into C++ is that a cell array of strings (really character arrays) is passed to C++ as a vector of strings. A vector of strings is returned as cell array of character vectors. There is, however, a rather general ability to return a C++ vector of  ```std::variant``` of any set of types the library can handle back to MATLAB/Octave as a cell array. See below for an example.
```cpp
#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <variant>
#include <typeindex>
#include <map>
#include <iostream>

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

  MexUnpacker<std::vector<std::string>> my_unpack(nrhs, prhs);
  try {
    auto [string_in] = my_unpack.unpackMex();

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

    MexPacker<std::vector<std::variant<userStruct1,userStruct2,double>>,std::vector<std::string>>my_pack(nlhs, plhs,my_name_map);
    my_pack.PackMex(variant_vec,string_in);

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
```

Example MATLAB/Octave code using this
```matlab
>> [a,b]=example_6({'abcd','efgh'})
abcd
efgh
a =
{
  [1,1] =

    scalar structure containing the fields:

      my_int = 1
      my_double =  4.5000

  [2,1] =

    scalar structure containing the fields:

      my_int = 1
      my_double =  4.5000

  [3,1] =

    scalar structure containing the fields:

      my_int = 1
      my_string = doc

  [4,1] =

    scalar structure containing the fields:

      my_int = 1
      my_string = doc

  [5,1] =  3
}
b =
{
  [1,1] = abcd
  [2,1] = efgh
}
```

### Additional Types

The full set of fixed integer width types (int32,uint32,int16,uint16,int8,uint8) are supported, however, currently they can only be passed into c++ as a scalar or a [ptr_tuple](#type-aliases) (tuple with the pointer, row size, column size).

```cpp
#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>

using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


  MexUnpacker<uint32_t,  ptr_tuple<uint32_t>,ptr_tuple<uint8_t>,ptr_tuple<int16_t>> my_unpack(nrhs, prhs);//We need the pointers to the inputs from matlab/octave

  try {

    auto [a,b,c,d] = my_unpack.unpackMex(); 

    a = 16;

    MexPacker<uint32_t,  ptr_tuple<uint32_t>,ptr_tuple<uint8_t>,ptr_tuple<int16_t> > my_pack(nlhs, plhs); //Create a packing object
    my_pack.PackMex(a,b,c,d); //Return these objects to matlab

// unpackMex can throw if the matlab objects passed at run time
// are not compatible with what is expected
  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
```
Example MATLAB/Octave code.
```matlab
>> [a,b,c,d]=example_7(uint32(1),uint32([1,2,3]),uint8([1,2,3]),int16([5.0,6]))
a = 16
b =

  1  2  3

c =

  1  2  3

d =

  5  6
```

For returning arrays of fixed width types (as well as floats and doubles) back to MATLAB/Octave they can be returned as a  [ptr_tuple](#type-aliases) of raw pointers, shared pointers, unique pointers,  an array, or as a ```std::vector```. Note that for arrays the array size must be included in the template type. Also, note that due to limitations of Boost PFR, while pointer tuples inside structs will correctly be returned to MATLAB/Octave as structs with array members, arrays (of the form ```T[N]```) inside of structs will generate a compilation error that the struct is not a simple aggregate.
```cpp
#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>

using namespace MexPackUnpackTypes;

struct userStruct{
  int a;
  double b;
  unique_ptr_tuple<uint32_t> c;
  //uint32_t c[5]; //Would be a compilation error
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

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
```

Example MATLAB/Octave code
```matlab
>> [a,b,c,d]=example_8
a =

  1
  2
  3
  4
  5

b =

   0
   1
   8
  27
  64

c =

  scalar structure containing the fields:

    field_0 = 1
    field_1 =  3
    field_2 =

       0
       1
       4
       9
      16
d =

   10
   11
   12
   13
   14
```

### Multi-dimensional arrays

There is limited support for 3-dimensional arrays (higher dimensional arrays are not currently supported). 3-dimensional arrays can be passed to c++ as [```ptr_tuple_3dim<S>```](#type-aliases)  (pointer plus the dimensions) where ```S``` is float, int, etc. For complex arrays, when using MATLAB 2018b and later (with -R2018a) the same type can be used as  ```ptr_tuple_3dim<std::complex<S>>``` where ```S``` is double or float. 
For complex arrays (see below) when using Octave and MATLAB before 2017b the ```ptr_tuple_3dim_CS<S>```  with ```S``` as double or float will handle the separate real and imaginary pointers as a pair of pointers in the first component of the tuple. Here is a very simple example that just takes the inputs and passes them back.
```cpp
#include "MexPackUnpack.h"
#include "mex.h"

using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  try {

    MexUnpacker<ptr_tuple_3dim<double>,ptr_tuple_3dim_CS<double> > my_unpack(nrhs, prhs);//We need the pointers to the inputs from matlab/octave
    auto [a,b] = my_unpack.unpackMex();

    MexPacker<ptr_tuple_3dim<double> ,ptr_tuple_3dim_CS<double> > my_pack(nlhs, plhs); //Create a packing object
    my_pack.PackMex(a,b); 

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
```

```matlab
>> m=randn(3,4,5)+i*randn(3,4,5);
>> [a,b]=example_9(real(m),m);
```

#### Eigen Vector Treatment and Least Squares Example

Currently Eigen Vector objects (e.g. ```Eigen::VectorXd```) are not explicitly supported by the interface. But this should not pose any difficulty as they are effectively supported implicitly. If a MATLAB/Octave array is passed in as an Eigen Matrix, you can get an explicit Eigen Vector using the row or column method. Any Eigen Vector that needs to be returned to MATLAB/Octave, can be passed using the corresponding Eigen Matrix type template parameter. This works because Eigen Vectors are basically just Eigen Matrices with one dimension size specialized at compile time to 1. Here is an example passing a matrix of points into C++, doing a least squares fit in C++ (obviously a trivial example that could have just been done in MATLAB/Octave) and returning the parameter vector from the least squares fit to MATLAB/Octave. Note the parameters are an Eigen Vector object but we can return them using the Eigen Matrix template parameter. Also note that even though we have a single input to the Mex function, it is returned from unpackMex as as tuple (with one element) and so needs to be destructured.

```cpp
#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <iostream>
#include <numeric>
#include <cmath>
using namespace MexPackUnpackTypes;

Eigen::MatrixXd LinearFit(const Eigen::MatrixXd & points){
  Eigen::MatrixXd m(points.cols(),2);
  m.col(0) = Eigen::VectorXd::Ones(points.cols());
  m.col(1) = points.row(0);
  Eigen::VectorXd v = m.colPivHouseholderQr().solve(points.row(1).transpose());
  return v;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  MexUnpacker<EDRM> my_unpack(nrhs, prhs);

  try {

    auto [points] = my_unpack.unpackMex(); //points is assumed to be 2 x Number of points matrix
    Eigen::VectorXd e = LinearFit(points);
    MexPacker<Eigen::MatrixXd> my_pack(nlhs, plhs); 
    my_pack.PackMex(e);    // We return an eigen vector through a matrix type. An Eigen vector is just a Matrix with
                           // one dimension (template parameter) fixed to size 1.

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
```

```matlab
>> x=randn(1,300);
>> y=10*x+5+randn;
>> m=[x;y];
>> b=eigen_vector(m);
>> b
b =

   5.2569
  10.0000
```


## Comments on Efficiency
Passing data into C++ is very efficient for pure numeric arrays as the underlying matlab pointer is simply wrapped one way or another. However, for simplicity in handling returning data from C++ to MATLAB/Octave in a uniform fashion, currrently data returned to MATLAB has to be copied into a new MATLAB allocated array and so incurs a copy. In the future we should probably have a separate type that wraps a MATLAB allocated array which could then be returned without forcing a copy. As a workaround for now, one way to avoid this copy is to allocate the output data in MATLAB/Octave and pass that as an additional input which is then written to in C++ (and changes to this array in C++ will be reflected on the MATLAB side without an explicit return). Technically this is not supported by MATLAB, however, as long as only a single reference to the array exists on the MATLAB side it is generally safe. For complicated arrays of structs, many calls to the matlab API will be made to get the subfields out and put them back in when data is passed into or out of C++; so depending on the use case, complicated struct arrays with many entries and fields may have some overhead.

## Notes on Complex Numeric Types

Octave and MATLAB up to R2017b internally store complex matrices as separate arrays (pointers) for the real and imaginary component. MATLAB from 2018a onward internally stores complex arrays with real and imaginary data interleaved and exposes a different API for complex matrices if mex functions are compiled with the -R2018a option. This library supports both in principle. There are #IFDEF statements checking for MX_HAS_INTERLEAVED_COMPLEX which will be set and true for MATLAB R2018a (if compiled with -R2018a) and false or not set at all otherwise. 

If using MATLAB after R2018a with -R2018a, for complex arrays the EDCIM type (eigen double complex interleaved map) which is just an eigen array of ```std::complex<double>```  and CDIP type (complex double interleaved pointer) which exposes a raw  ```std::complex<double>*``` can be used as well as float variants. Note that CDIP is the same as ```ptr_tuple<std::complex<double>> ```. Similarly ```ptr_tuple_3dim<std::complex<S>>``` works fine for ```S``` double or float. 

For MATLAB 2017b and earlier and octave you can use EDSCM (eigen double split complex) which uses a pair of real Eigen maps or CDSP (complex double split pointer) which extracts a pair of real pointers or the single precision variants. There is a   ```ptr_tuple_3dim_CS<S>``` for 3 dimensional complex arrays with split representation where the first components is a pair of pointers. 

In either the split or interleaved case ```std::vector<std::complex<float>>``` and ```std::vector<std::complex<double>>``` can be returned to matlab as 1D arrays.



## Unsupported MATLAB/Octave types
Passing cell arrays from MATLAB/Octave to C++ is not currently supported (with the exception of cell arrays of strings), but creating cell arrays in MATLAB/Octave from C++ is supported as described above. Only 1D, 2D, and 3D arrays are currently supported. Note that in all cases when we say MATLAB/Octave "strings", we mean arrays of characters. Before Matlab 2018 this was the only option. Now there is a separate string type in post-2018 matlab that is not an array of characters. For backwards compatibility this newer style string is not supported.

## Why not MATLAB C++ API
Since 2018a MATLAB has an alternate C++ based interface that is different than the C API used by earlier matlab versions and octave. I have stuck with wrapping the C interface for backwards compatibility with older matlab installations and octave. 


## Contributing
Pull requests, comments, or suggestions are welcome. 

## License
[MIT](LICENSE.txt)