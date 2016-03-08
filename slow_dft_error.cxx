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
#define apd_time_res 164.61e-12
#define binsize 100

using namespace std;

//---Other header files---//

#include "read_inputs.h"


int main() {

  //------Reading inputs from g2 file------//

  fstream g2_file;

  string input = "";

  string ini = "";

  string directory = "";

  char ans_samp_time = {0};

  char ans_constrains = {0};

  double time, counts, corr, corr_error, samp_time, low_bound, up_bound;;

  int first_file_number, last_file_number;

  vector<double> g2_time;

  vector<double> g2;

  vector<double> g2_error;


  //------Data selection------//

  cout << "type in the address of the directory where the data is held (e.g. C:\\data\\ or /home/user/data/)\n> ";

  getline(cin , directory);

  cout << "You entered: " << directory << endl << endl;


  //------Choosing right sampling time------//

  do

    safe_read_char("Use default sampling time? (cqed-lvis experiment) (y/n)\n> "
		   , ans_samp_time);

  while (ans_samp_time != 'y' && ans_samp_time != 'n');

  if (ans_samp_time == 'n')

    safe_read_num("Sampling time of data: " , samp_time);

  else

    samp_time = apd_time_res * binsize;


  //------Constraining data------//

  do

    safe_read_char("Constrain data ? (y/n)\n> " , ans_constrains);

  while (ans_constrains != 'y' && ans_constrains != 'n');

  if (ans_constrains == 'y') {

    safe_read_num("Lower bound (microseconds): " , low_bound);

    safe_read_num("Upper bound (microseconds): " , up_bound);

  }


  //------Data files to analyze------//

  cout << "WARNING: the data must be named as data#_g2.asc or data#_g2.txt"
       << endl;

  safe_read_num("begin with what data file?: " , first_file_number);

  safe_read_num("end with what data file?: " , last_file_number);

  for (int counter = first_file_number;
       counter <= last_file_number;
       counter++) {

    stringstream filename;

    stringstream filenameopt;

    stringstream filenamefinaloutput;

    filename << directory << "data" << counter << "_g2.asc";

    if (ans_constrains == 'y')

      filenamefinaloutput << directory << "data" << counter << "_dft_"
			  << low_bound << "to" << up_bound << ".asc";

    else

      filenamefinaloutput << directory << "data" << counter << "_dft.asc";

    cout << endl << "opening file " << filename.str() << endl;

    ifstream g2_file;

    g2_file.open (filename.str().c_str());

    if (!g2_file.is_open()) {

      cout << endl << "the file " << filename.str() << " could not be opened"
	   << endl;

      filenameopt << directory << "data" << counter << "_g2.txt";

      cout << endl << "opening file " << filenameopt.str() << endl;

      g2_file.open (filenameopt.str().c_str());

      if (!g2_file.is_open()) {

	cout << endl << "nor could file " << filenameopt.str() << endl;

	continue;
      }
    }

    int linenum = 0;

    while (! g2_file.eof()) {

      //------Sorting------//

      g2_file >> ini;

      linenum++;

      if (ini[0] == '#') {

	g2_file.ignore(256 , '\n');

	continue;

      } else

	conv_str_num(ini , time);

      g2_file >> counts;

      g2_file >> corr;

      g2_file >> corr_error;

      if (ans_constrains == 'n') {

	g2_time.push_back(time);

	g2.push_back(corr);

	g2_error.push_back(corr_error);

      } else if (time > low_bound && time < up_bound) {

	g2_time.push_back(time);

	g2.push_back(corr);

	g2_error.push_back(corr_error);

      }

      if (g2_file.eof())

	break;
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

    int N = g2.size();

    double *realpart = new double [N];

    double *imagpart = new double [N];

    double *mag = new double [N];

    double *mag2 = new double [N];

    double *phase = new double [N];

    double *freq = new double [N];

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

	realpart[ii] += g2[jj] * cos(kbeta);

	imagpart[ii] += g2[jj] * sin(kbeta);

	error2_real[ii] += pow(cos(kbeta) * g2_error[jj] , 2);

	error2_imag[ii] += pow(sin(kbeta) * g2_error[jj] , 2);

	error2_cross[ii] += sin(kbeta) * cos(kbeta) * pow(g2_error[jj] , 2);
      }

      mag2[ii] = (pow(realpart[ii] , 2) + pow(imagpart[ii] , 2));

      mag[ii] = sqrt(mag2[ii]);

      phase[ii] = atan2(imagpart[ii], realpart[ii]);

      error_mag[ii] = sqrt((pow(realpart[ii], 2) * error2_real[ii] + pow(imagpart[ii], 2) * error2_imag[ii] + 2 * realpart[ii] * imagpart[ii] * error2_cross[ii]) / pow(mag[ii], 2));

      error_mag2[ii] = 2 * mag[ii] * error_mag[ii];

      error_phase[ii] = sqrt((pow(realpart[ii], 2) * error2_imag[ii] + pow(imagpart[ii], 2) * error2_real[ii] - 2 * realpart[ii] * imagpart[ii] * error2_cross[ii]) / pow(mag[ii], 4));

      freq[ii] = ii / (N * samp_time);
    }

    //------Output file operations------//

    ofstream finaloutputfp;

    finaloutputfp.open(filenamefinaloutput.str().c_str() , ios::trunc);

    finaloutputfp << "#data " << counter << " DFT and error\n";

    finaloutputfp << "#freq\t";

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

    //------g2 DFT and errors------//

    for (int kk = 1; kk < N; kk++) {

      finaloutputfp << setprecision(8) << freq[kk] << "  ";

      finaloutputfp << setprecision(8) << mag[kk] << "  ";

      finaloutputfp << setprecision(8) << error_mag[kk] << "  ";

      finaloutputfp << setprecision(8) << mag2[kk] << "  ";

      finaloutputfp << setprecision(8) << error_mag2[kk] << "  ";

      finaloutputfp << setprecision(8) << phase[kk] << "  ";

      finaloutputfp << setprecision(8) << error_phase[kk] << "  ";

      finaloutputfp << setprecision(8) << realpart[kk] << "  ";

      finaloutputfp << setprecision(8) << sqrt(error2_real[kk]) << "  ";

      finaloutputfp << setprecision(8) << imagpart[kk] << "  ";

      finaloutputfp << setprecision(8) << sqrt(error2_imag[kk]) << "\n";
    }

    finaloutputfp.close();

    g2_time.clear();

    g2.clear();

    g2_error.clear();

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

    delete [] freq;
  }

  return 0;
}
