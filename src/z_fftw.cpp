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

#include "z_fftw.hpp"

void fourier_plus(std::complex<double> out[], std::complex<double> in[],
                  const int corr, const int nzeros) {
  fftw_plan p;
  for (int i=0; i<corr; i++)
    in[i] = in[i] - in[corr-1];
  int i_max = corr+nzeros;
  for (int i=corr; i<i_max; i++) {
    in[i] = 0.0;
    in[i+nzeros] = 0.0;
  }
  i_max = 2*nzeros+2*corr-2;
  for (int i = corr+2*nzeros - 1; i<i_max; i++)
    in[i] =  conj(in[2*(nzeros+corr)-2-i]);
  p = fftw_plan_dft_1d(2*(corr+nzeros)-2, reinterpret_cast<fftw_complex*>(in),
                       reinterpret_cast<fftw_complex*>(out), FFTW_FORWARD,
                       FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void fourier_noSymm_plus(std::complex<double> out[], std::complex<double> in[],
                         const int corr, const int nzeros) {
  fftw_plan p;
  for (int i=0; i<corr; i++)
    in[i] = 2.0*(in[i] - in[corr-1]);

  const int i_max = corr+nzeros;
  for (int i=corr; i<i_max; i++)
    in[i] = 0.0;

  p = fftw_plan_dft_1d(2*(corr+nzeros)-2, reinterpret_cast<fftw_complex*>(in),
                       reinterpret_cast<fftw_complex*>(out), FFTW_FORWARD,
                       FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void run_fftw(std::complex<double> out[], std::complex<double> in[],
              const int corr, const int nzeros) {
  fftw_plan p;
  int i_max = corr+nzeros;
  for (int i=corr; i<i_max; i++) {
    in[i] = 0.0;
    in[i+nzeros] = 0.0;
  }
  i_max = 2*nzeros+2*corr-2;
  for (int i = corr+2*nzeros - 1; i<i_max; i++)
    in[i] =  conj(in[2*(nzeros+corr)-2-i]);
  p = fftw_plan_dft_1d(2*(corr+nzeros)-2, reinterpret_cast<fftw_complex*>(in),
                       reinterpret_cast<fftw_complex*>(out), FFTW_FORWARD,
                       FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}
