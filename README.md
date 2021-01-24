# MexPackUnpack

MexPackUnpack is a C++ header library to help automate extracting data from MATLAB/Octave data passed into user C++ mex code. It uses some lite (depending on your point of view) template metaprogramming with c++ variadic templates. It requires at least C++14 but as shown below, is intended to be used in conjunction with C++17 structured binding syntax. 

The current version extracts the MATLAB/Octave numeric arrays to either [Eigen](https://eigen.tuxfamily.org/) Map objects (Eigen Matrices wrapping an external  pointer) or a tuple of the underlying pointers and array sizes. MATLAB strings can be converted to C++ std::strings and Matlab structs. Struct arrays can be extracted into user defined C++ structs with a matching layout. These types can then be returned and converted to MATLAB/Octave objects.


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

Code using the resulting mex function:
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
>> [out1,out2] = example_1(m,n,2.0,"More stuff")
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

The library handle complex as well as real matrices. There is some complexity due to the fact that matlab 2017b and earlier as well as octave use a split real/imaginary representation while matlab 2018a and later uses an interleaved representation. See discussion below about this. This example show use of split complex types. 

```cpp
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

//Wrapper for dealing with MATLAB/Octave Structs
template<class T>
  struct MXStruct{
    const T& obj;
    MXStruct(T& in) : obj{in} {}
  };

//Wrapper for dealing with MATLAB/Octave Struct array
template<class T>
  struct MXVecStruct{
    const std::vector<T>& obj;
    MXVecStruct(std::vector<T>& in) : obj{in} {}
  };

}
```


### Structs and Struct Arrays

To deal with matlab structs and struct arrays we can  define a user C++ structure and then use that in the type list to automatically extract matlab structs which have fields matching the user specified C++ struct. This is done using the reflection/introspection capability of  [Boost PFR](https://apolukhin.github.io/magic_get/index.html). To indicate that we should unpack into the user type, the user type needs to be "wrapped" in the  [```MXStruct```](#type-aliases) type.

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


struct userStruct2{
  int a;
  double b;
  std::string c; 
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


  //To indicate that a struct type will be used to receive a matlab struct it needs to be "wrapped"
  //in MXStruct
  MexUnpacker<MXStruct<userStruct>> my_unpack(nrhs, prhs);

  try {

    auto [d] = my_unpack.unpackMex(); //d is a userStruct
    d.a=4;

    MexPacker<MXStruct<userStruct> >my_pack(nlhs, plhs);
    my_pack.PackMex(MXStruct(d)); //MXStruct just stores a reference to the actual struct d to be returned

  } catch (std::string s) {
    mexPrintf(s.data());
  }

}
```
Note that while we can determine the types of the user defined fields in the C++ struct using BOOST PFR, we cannot determine the names of the fields from within C++. So currently when structs are returned they are just named field_0, field_1, ..., etc. If it was needed you could probably add fields to the MXStruct class to indicate what the fieldnames of the returned struct should be.

```matlab
>> m=struct('a',int32(1),'b',3.7,'c',[1.0,2.0,5.0],'d',"what's up");
>> a=example_3(m)
a =

  scalar structure containing the fields:

    field_0 = 4
    field_1 =  3.7000
    field_2 =

       1   2   5

    field_3 = what's up
```

#### Struct Array Example
Struct arrays with more than 1 element can be accomodated by wrapping the user struct in ```MXVecStruct```. 
```cpp
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
```
Example octave code calling it.
```
>> m=struct('a',int32(1),'b',3.7,'c',[1.0,2.0,5.0],'d',"what's up");
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

    field_0
    field_1
    field_2

>> a(2)
ans =

  scalar structure containing the fields:

    field_0 = 1
    field_1 =  3.7000
    field_2 = doc
```

#### Slightly Ridiculous (Mis)use of Struct/Struct Array Functionality
It's actually possible to have nested structs in C++ returned as nested structs in matlab, however, this is probably not a terribly useful feature currently due to the fact that fieldnames cannot be set from the C++ side and the type wrapping with MXStruct becomes a bit cumbersome. In principle this should work for arbitrary levels of nesting (until the compiler gets upset at least).
```cpp
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
```
Example octave code using this.
```matlab
>> m=struct('a',int32(1),'b',3.7,'c',[1.0,2.0,5.0],'d',"whats' up");
>> m=[m,m,m,m]
m =

  1x4 struct array containing the fields:

    a
    b
    c
    d

>> a=example_5(m)
a =

  scalar structure containing the fields:

    field_0 = 1
    field_1 =  4.5000
    field_2 = doc
    field_3 =

      4x1 struct array containing the fields:

        field_0
        field_1
        field_2
        field_3


>> a.field_3(1)
ans =

  scalar structure containing the fields:

    field_0 = 1
    field_1 =  3.7000
    field_2 =

       1   2   5

    field_3 = whats' up
```
However, note that the reverse functionality (ingesting nested MATLAB/octave structs) does not work. This is because a nested struct would need to have a member that was ```MXStruct<S>``` which the code would try to unpack an ```S``` into. If you really wanted this to work, one way to deal with this might be to eliminate the need to wrap struct with ```MXStruct``` by just checking that a user specified type was not on the list of expected numeric types and was POD. 

## Notes on Complex Numeric Types

Octave and Matlab up to R2017b internally store complex matrices as separate arrays (pointers) for the real and imaginary component. Matlab from 2018a onward internally stores complex arrays with real and imaginary data interleaved and exposes a different API for complex matrices if mex functions are compiled with the -R2018a option. This library supports both in principle. There are #IFDEF statements checking for MX_HAS_INTERLEAVED_COMPLEX which will be set and true for MATLAB R2018a (if compiled with -R2018a) and false or not set at all otherwise. If using MATLAB after R2018a with -R2018a, for complex arrays the EDCIM type (eigen double complex interleaved map) which is just an eigen array of ```std::complex<double>```  and CDIP type (complex double interleaved pointer) which exposes a raw  ```std::complex<double>*``` can be used, otherwise you can use EDSCM (eigen double split complex) which uses a pair of real Eigen maps or CDSP (complex double split pointer) which extracts a pair of real pointers. Due to lack of access to more recent MATLAB versions on the systems where this library was developed the interleaved complex functionality is not well tested and may require some tweaks or bug fixes.


## Compilation

Eigen and boost::pfr are included as git submodules and need to be in the include path when compiling.

```matlab
>> mkoctfile --mex -v -std=c++17 ./examples/example_1.cpp  -I./eigen/Eigen -I./pfr/include/
```
The library itself compiles with C++14 but the destructuring of the tuples is less elegant as structured binding syntax is not available. This could be cleaned up a bit with some helper functions. 
The requirement to have Eigen available could be separated out and made optional. Also the library could be extended so that matlab/Octave objects could be extracted into other similar types using the same mechanism.

## Unsupported MATLAB/Octave types
At the moment only 1D and 2D arrays are supported, multidimensional arrays should probably be added. Also, currently cell arrays are not supported, only because I have never found a need to use them in mex code. 

## Why not MATLAB C++ API
Since 2018a MATLAB has an alternate C++ based interface that is different than the C API used by earlier matlab versions and octave. I have stuck with wrapping the C interface for backwards compatibility with older matlab installations and octave. 


## Contributing
Pull requests, comments, or suggestions are welcome. 

## License
[MIT](LICENSE.txt)