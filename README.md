# slow_dft_error - Calculates (slowly) the DFT of a dataset and propagates the uncertainties

## SHORT DESCRIPTION:

Brute force method for calculating the DFT. It brings the benefit of a straight forward way of propagating the errors to the frequency domain. The method follows the calculations of: Giovanni Betta, Consolatina Liguori, Antonio Pietrosanto, Propagation of uncertainty in a discrete Fourier transform algorithm, Measurement, Volume 27, Issue 4, June 2000, Pages 231-239, ISSN 0263-2241, DOI: 10.1016/S0263-2241(99)00068-8

## LICENSE:

This program is released under the terms of the GNU Public License v3 (GPLv3). For more information please see:
http://www.gnu.org/licenses/gpl-3.0.html

## HOW TO USE:

At present, the program can only handle input ascii text files which contain 4 columns of data. First column displaying the time, the second one showing the raw count, the third containing the normalized values and the fourth displaying the error or uncertainty of the third column (this was the data file specs of my experiment at the time).

It requires several user inputs. Namely the sampling time and the data constraints (optional).

It is designed to accept input of the form dataX_g2.asc. Where X can be any positive integer. In this way, it can be run for multiple datasets.

Its output is also a text file. As a header, it displays the dataset number. Then shows columns for frequency (Hz), Magnitude (Mag), Mag error, Mag squared, Mag2 error, phase and phase error.

To run the program on a Unix-based system, open a terminal and then type:

./slow_dft_error

## COMPILATION:

It is written in plain C++ using only the standard libraries. We have only thoroughly tested g++ (The GNU C++ compiler) as compiler, on a GNU/Linux 64 bits machine. We have not tested Mac OSX, but given g++ is easily available for this platform, it should be easy to produce an executable for it.

### FOR GNU/LINUX (AND POTENTIALLY MACOSX):

In our GNU/Linux machine, to compile we simply do:

g++ -Wall slow_dft_error.cxx read_inputs.cxx -o slow_dft_error
(Only if the read_inputs.cxx and read_inputs.h files are located in the same directory)

or more generally

g++ -Wall -I/path-to-header-file gfncalc.cxx /path-to-file/read_inputs.cxx -o slow_dft_error

This binary has been tested on Debian GNU/Linux systems version 5/6/7 and CentOS 5. However, it should work on any Unix-based OS.
