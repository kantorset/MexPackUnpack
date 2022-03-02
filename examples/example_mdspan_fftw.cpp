#include "MexPackUnpack.h"
#include "mex.h"
#include <Eigen>
#include <iostream>
#include <complex.h>
#include <fftw3.h>

using namespace MexPackUnpackTypes;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Create an unpacker with template parameters corresponding to what we expect from matlab
  // argument 1: 4 dimensional mdspan

  MexUnpacker<span_4d_dynamic_left<double>  > my_unpack(nrhs, prhs); //We need the pointers to the inputs from matlab/octave

  try {

    auto [a] = my_unpack.unpackMex(); 

    if(a.extent(0)<2||a.extent(1)<6||a.extent(2)<4||a.extent(3)<2){
      throw std::string("Input is too small for slice indices");
    }

    //In this example we take a real input and copy it to a complex array
    //We will do a full complex FFT not the specialized real to complex variants FFTW provides
    std::unique_ptr<std::complex<double>[]>a_complex(new std::complex<double>[a.size()]);
    std::copy(a.data(),a.data()+a.size(),a_complex.get());

    //mdspan view of complex data as multidimensional array
    auto a_complex_span = span_4d_dynamic_left<std::complex<double> > {a_complex.get(), a.extents()};

    auto c = stdex::submdspan(a_complex_span,1,std::pair{2ul,6ul},3,stdex::full_extent); //Create a slice (submdspan)
        
    //2D FFT of slice

    fftw_plan p;
    std::array<int,2> n = {c.extent(1),c.extent(0)}; // FFT size to pass to fftw, extents are swapped because of row vs column major ordering
    int rank = 2; //2D FFT
    int howmany = 1; //Just one

    //space to put FFT output 
    std::unique_ptr<std::complex<double>[]>output(new std::complex<double>[c.extent(1)*c.extent(0)]);
    auto out_span =span_2d_dynamic_left<std::complex<double> > {output.get(),std::array<int,2>{c.extent(0),c.extent(1)}}; //span view of output

    int inembed[2] = {c.stride(0),c.stride(1)/c.stride(0)}; //embedding parameters for fftw, tells it how to stride through data and access the 2D subspan

    int idist = 0 ; //doesn't matter since just doing 1 FFT
    int ostride = 1; //output stride (contiguous)
    int odist = 0; //doesn't matter just 1
    
    //Note that c.stride(0) is the input inner stride
    p = fftw_plan_many_dft(rank, n.data(), howmany,
                            reinterpret_cast<fftw_complex*>( c.data() ) , inembed,
                             c.stride(0), idist,
                             reinterpret_cast<fftw_complex*>(output.get()), n.data(),
                             ostride, odist,
                             FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p); 
    fftw_destroy_plan(p);    
    
    //Now do 1D FFT of columns of subspan

    //allocate new output
    std::unique_ptr<std::complex<double>[]>output_1(new std::complex<double>[c.extent(1)*c.extent(0)]);
    auto out_span_1 =stdex::mdspan<std::complex<double>, stdex::extents<stdex::dynamic_extent,stdex::dynamic_extent>,stdex::layout_left>
    (output_1.get(),std::array<int,2>{c.extent(0),c.extent(1)}); //span view of output

    rank = 1; //1D FFT
    howmany = c.extent(1); //1 for each column
    idist = c.stride(1);  //input spacing stride
    odist = c.extent(0);  //each output is separated by one column
    
    //Now doing 1D FFTs the size of each column (n[1] will be ignored)
    n[0]=c.extent(0);
    //Note that c.stride(0) is the input inner stride
    p = fftw_plan_many_dft(rank, n.data(), howmany,
                            reinterpret_cast<fftw_complex*>( c.data() ) , inembed,
                             c.stride(0), idist,
                             reinterpret_cast<fftw_complex*>(output_1.get()), NULL,
                             ostride, odist,
                             FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p); 
    fftw_destroy_plan(p);    
    
  
    
    //1D FFT rows
    //allocate output
    std::unique_ptr<std::complex<double>[]>output_2(new std::complex<double>[c.extent(1)*c.extent(0)]);

    //view as a span
    auto out_span_2 =stdex::mdspan<std::complex<double>, stdex::extents<stdex::dynamic_extent,stdex::dynamic_extent>,stdex::layout_left>
    (output_2.get(),std::array<int,2>{c.extent(0),c.extent(1)}); //Create a slice (submdspan)

    rank = 1; //1D FFT
    howmany = c.extent(0); //1 for each row
    idist = c.stride(0); //spacing is same is subspan inner stride
    odist = 1; //outputs are contiguous
    ostride = c.extent(0); // But output stride is equal to subspan column size

    //Now doing 1D FFTs the size of each row (n[1] will be ignored)
    n[0]=c.extent(1); 
    p = fftw_plan_many_dft(rank, n.data(), howmany,
                            reinterpret_cast<fftw_complex*>( c.data() ) , inembed,
                             c.stride(1), idist,
                             reinterpret_cast<fftw_complex*>(output_2.get()), NULL,
                             ostride, odist,
                             FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);    
    

    MexPacker<decltype(a_complex_span),decltype(c),decltype(out_span),decltype(out_span_1),decltype(out_span_2) > my_pack(nlhs, plhs); //Create a packing object, we use decltype to get the strides from submdspan correct
    my_pack.PackMex(a_complex_span,c,out_span,out_span_1,out_span_2); //Return this object to matlab    

  } catch (std::string s) {
    mexPrintf(s.data());
  }
  
}
