//Copyright (c) 2015 Zachary Kann
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

// ---
// Author: Zachary Kann

#include <cmath>
#include <fstream>
#include <algorithm>
#include <complex>
#include <fftw3.h>
#include "boost/program_options.hpp"
#include "z_string.hpp"
#include "z_fftw.hpp"
#include "z_constants.hpp"

namespace po = boost::program_options;

int main (int argc, char *argv[]) {

  po::options_description desc("Options");
  desc.add_options()
    ("help,h",  "Print help messages")
    ("input,i", po::value<std::string>()->required(), "Name of input file")
    ("output,o", po::value<std::string>()->required(), "Name of output file");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(EXIT_SUCCESS);
  }

  std::ifstream input_file(vm["input"].as<std::string>().c_str());
  std::ofstream output_file(vm["output"].as<std::string>().c_str());

  int corr = std::count(std::istreambuf_iterator<char>(input_file),
                        std::istreambuf_iterator<char>(), '\n');

  int spec_length = 2*corr-2;
  std::complex<double> *spectraIn = new std::complex<double>[spec_length];
  std::complex<double> *spectraOut = new std::complex<double>[spec_length];
  std::vector<double> time(corr);

  for (int i = 0; i < corr; ++i) {
    input_file >> time[i];
    input_file >> spectraIn[i];
  }

  double deltaT = time[1] - time[0];

  fourier_plus(spectraOut, spectraIn, corr);

  double freq_factor = 1.0/2.0/C_SPEED/(deltaT*corr);

  int max_i = 2*corr-2;
  for (int i = corr; i < max_i; ++i)
    output_file << (i-max_i)*freq_factor << std::real(spectraOut[i]) <<
                   std::imag(spectraOut[i]);
  for (int i = 0; i < corr; ++i)
    output_file << i*freq_factor << std::real(spectraOut[i]) <<
                   std::imag(spectraOut[i]);

  output_file.close();
}
 // main
