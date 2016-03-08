#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>

#define pi 3.14159265

using namespace std;

//---Other header files---//

#include "read_inputs.h"


int main() {

  //------Reading inputs from dft file------//

  fstream dft_file;

  string input = "";

  string ini = "";

  string directory = "";

  //  char ans_samp_time = {0};

  char ans_constrains = {0};

  double freq, dftmag, dftmag_err, dftmag2, dftmag2_err, phi, phi_err;

  double repart, repart_err, impart, impart_err, samp_rate, low_bound, up_bound;

  int first_file_number, last_file_number;

  vector<double> dft_freq;

  vector<double> real;

  vector<double> real_error;

  vector<double> imag;

  vector<double> imag_error;


  //------Data selection------//

  cout << "type in the address of the directory where the data is held (e.g. C:\\data\\ or /home/user/data/)\n> ";

  getline(cin , directory);

  cout << "You entered: " << directory << endl << endl;


  //------Choosing right sampling time------//
  /*
  do

    safe_read_char("Use default sampling time? (cqed-lvis experiment) (y/n)\n> "
		   , ans_samp_time);

  while (ans_samp_time != 'y' && ans_samp_time != 'n');

  if (ans_samp_time == 'n')
  */
  // safe_read_num("Sampling rate for DFT: " , samp_rate);
  /*
  else

    samp_time = apd_time_res * binsize;
  */

  //------Constraining data------//

  do

    safe_read_char("Exclude data ? (y/n)\n> " , ans_constrains);

  while (ans_constrains != 'y' && ans_constrains != 'n');

  if (ans_constrains == 'y') {

    safe_read_num("Lower bound (Hz): " , low_bound);

    safe_read_num("Upper bound (Hz): " , up_bound);

  }


  //------Data files to analyze------//

  cout << "WARNING: the data must be named as data#_dft.asc or data#_dft.txt"
       << endl;

  safe_read_num("begin with what data file?: " , first_file_number);

  safe_read_num("end with what data file?: " , last_file_number);

  for (int counter = first_file_number;
       counter <= last_file_number;
       counter++) {

    stringstream filename;

    stringstream filenameopt;

    stringstream filenamefinaloutput;

    filename << directory << "data" << counter << "_dft.asc";

    if (ans_constrains == 'y')

      filenamefinaloutput << directory << "data" << counter << "_g2ft_"
			  << low_bound << "to" << up_bound << ".asc";

    else

      filenamefinaloutput << directory << "data" << counter << "_g2ft.asc";

    cout << endl << "opening file " << filename.str() << endl;

    ifstream dft_file;

    dft_file.open (filename.str().c_str());

    if (!dft_file.is_open()) {

      cout << endl << "the file " << filename.str() << " could not be opened"
	   << endl;

      filenameopt << directory << "data" << counter << "_dft.txt";

      cout << endl << "opening file " << filenameopt.str() << endl;

      dft_file.open (filenameopt.str().c_str());

      if (!dft_file.is_open()) {

	cout << endl << "nor could file " << filenameopt.str() << endl;

	continue;
      }
    }

    int linenum = 0;

    while (! dft_file.eof()) {

      //------Sorting------//

      dft_file >> ini;

      linenum++;

      if (ini[0] == '#') {

	dft_file.ignore(256 , '\n');

	continue;

      } else

	conv_str_num(ini , freq);

      dft_file >> dftmag;

      dft_file >> dftmag_err;

      dft_file >> dftmag2;

      dft_file >> dftmag2_err;

      dft_file >> phi;

      dft_file >> phi_err;

      dft_file >> repart;

      dft_file >> repart_err;

      dft_file >> impart;

      dft_file >> impart_err;

      dft_freq.push_back(freq);

      real.push_back(repart);

      real_error.push_back(repart_err);

      imag.push_back(impart);

      imag_error.push_back(impart_err);

      if (dft_file.eof())

	break;
    }

    /****Excluding data: Fourier analysis filtering****/

    int N = dft_freq.size();

    double fmax = dft_freq[N-1];

    if (ans_constrains == 'y') {

      for (int ll = 0; ll < N; ll++) {

	if ((dft_freq[ll] > low_bound && dft_freq[ll] < up_bound)
	    || (dft_freq[ll] > fmax - up_bound
		&& dft_freq[ll] < fmax - low_bound)) {

	  real[ll] = 0;

	  real_error[ll] = 0;

	  imag[ll] = 0;

	  imag_error[ll] = 0;

	}

      }

    }

    //------Discrete Fourier Transform calculation------//

    /*Brute force method for calculating the DFT. It brings the benefit of
      a straight forward way of propagating the errors to the frequency
      domain. The method follows the calculations of:
      Giovanni Betta, Consolatina Liguori, Antonio Pietrosanto,
      Propagation of uncertainty in a discrete Fourier transform algorithm,
      Measurement, Volume 27, Issue 4, June 2000,
      Pages 231-239, ISSN 0263-2241,
      DOI: 10.1016/S0263-2241(99)00068-8 */

    samp_rate = dft_freq[1] - dft_freq[0];

    double *realpart = new double [N];

    double *imagpart = new double [N];

    double *mag = new double [N];

    double *mag2 = new double [N];

    double *phase = new double [N];

    double *time = new double [N];

    double *error2_real = new double [N];

    double *error2_imag = new double [N];

    double *error2_cross = new double [N];

    double *error_mag = new double [N];

    double *error_mag2 = new double [N];

    double *error_phase = new double [N];

    for (int ii = 0; ii < N; ii++) {

      realpart[ii] = 0;

      imagpart[ii] = 0;

      error2_real[ii] = 0;

      error2_imag[ii] = 0;

      error2_cross[ii] = 0;

      for (int jj = 0; jj < N; jj++) {

	double kbeta = 2 * pi * ii * jj / N;

	realpart[ii] += real[jj] * cos(kbeta) - imag[jj] * sin(kbeta);

	imagpart[ii] += real[jj] * sin(kbeta) + imag[jj] * cos(kbeta);

	double real_err2_jth = pow(cos(kbeta) * real_error[jj], 2)
	  + pow(sin(kbeta) * imag_error[jj], 2);

	double imag_err2_jth = pow(sin(kbeta) * real_error[jj] , 2)
	  + pow(cos(kbeta) * imag_error[jj] , 2);

	error2_real[ii] += real_err2_jth;

	error2_imag[ii] += imag_err2_jth;

       	error2_cross[ii] += sqrt(real_err2_jth) * sqrt(imag_err2_jth);
      }

      mag2[ii] = (pow(realpart[ii] , 2) + pow(imagpart[ii] , 2)) / N;

      mag[ii] = sqrt(mag2[ii]);

      phase[ii] = atan2(imagpart[ii], realpart[ii]);

      error_mag[ii] = sqrt((pow(realpart[ii], 2) * error2_real[ii] + pow(imagpart[ii], 2) * error2_imag[ii] + 2 * realpart[ii] * imagpart[ii] * error2_cross[ii]) / ((pow(N, 2) * pow(mag[ii], 2))));

      error_mag2[ii] = 2 * mag[ii] * error_mag[ii];

      error_phase[ii] = sqrt((pow(realpart[ii], 2) * error2_imag[ii] + pow(imagpart[ii], 2) * error2_real[ii] - 2 * realpart[ii] * imagpart[ii] * error2_cross[ii]) / ((pow(N, 2) * pow(mag[ii], 4))));

      time[ii] = ii / (N * samp_rate) - 1 / (2 * samp_rate);
    }

    //------Output file operations------//

    ofstream finaloutputfp;

    finaloutputfp.open(filenamefinaloutput.str().c_str() , ios::trunc);

    finaloutputfp << "#data " << counter << " g2 and error using DFT filter\n";

    finaloutputfp << "#time\t";

    finaloutputfp << "mag\t";

    finaloutputfp << "error\t";

    finaloutputfp << "mag2\t";

    finaloutputfp << "error\t";

    finaloutputfp << "phase\t";

    finaloutputfp << "error\t";

    finaloutputfp << "real\t";

    finaloutputfp << "error\t";

    finaloutputfp << "imag\t";

    finaloutputfp << "error\n";


    //------g2 and errors------//

    for (int kk = 1; kk < N; kk++) {

      finaloutputfp << setprecision(8) << time[kk] * 1.0e6 << "  ";

      finaloutputfp << setprecision(8) << mag[kk] << "  ";

      finaloutputfp << setprecision(8) << error_mag[kk] << "  ";

      finaloutputfp << setprecision(8) << mag2[kk] << "  ";

      finaloutputfp << setprecision(8) << error_mag2[kk] << "  ";

      finaloutputfp << setprecision(8) << phase[kk] << "  ";

      finaloutputfp << setprecision(8) << error_phase[kk] << "  ";

      finaloutputfp << setprecision(8) << realpart[kk] / sqrt(N) << "  ";

      finaloutputfp << setprecision(8) << sqrt(error2_real[kk] / N) << "  ";

      finaloutputfp << setprecision(8) << imagpart[kk] / sqrt(N) << "  ";

      finaloutputfp << setprecision(8) << sqrt(error2_imag[kk] / N) << "\n";
    }

    finaloutputfp.close();

    dft_freq.clear();

    real.clear();

    real_error.clear();

    imag.clear();

    imag_error.clear();

    delete [] realpart;

    delete [] imagpart;

    delete [] error2_real;

    delete [] error2_imag;

    delete [] error2_cross;

    delete [] error_mag;

    delete [] error_mag2;

    delete [] error_phase;

    delete [] mag;

    delete [] mag2;

    delete [] phase;

    delete [] time;
  }

  return 0;
}
