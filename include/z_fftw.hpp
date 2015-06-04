// Copyright (c) 2015 Zachary Kann
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// ---
// Author: Zachary Kann

// A collection of functions to run fourier transforms using the fftw library.
// TODO(Zak): incorporate these into TCF class.
#include <complex>
#include <fftw3.h>

#ifndef _Z_FFTW_HPP_
#define _Z_FFTW_HPP_

// The standard call to prep and run a fourier transform. Assumes symmetry
// around t = 0.
extern void fourier_plus(std::complex<double> out[], std::complex<double> in[],
                         const int corr, const int nzeros = 0);

// Similar to fourier_plus, but without the assumption of symmetry around t = 0.
extern void fourier_noSymm_plus(std::complex<double> out[],
                                std::complex<double> in[], const int corr,
                                const int nzeros);

// Similar to fourier_plus, but doesn't shift the time correlation function
// to be zero after corr steps.
extern void run_fftw(std::complex<double> out[], std::complex<double> in[],
                     const int corr, const int nzeros);

#endif
