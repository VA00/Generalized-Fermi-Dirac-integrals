# Generalized-Fermi-Dirac-integrals
Provide reliable, accurate down to several ULP's and fast numerical framework to compute generalized Fermi-Dirac integrals aiming at full floating-point parameter coverage, generality, precision control and speed.


To install library, follow standard
Linux autoconf procedure (see also INSTALL):

./configure
make
sudo make install
sudo ldconfig

To install in e.g. home directory:

./configure --prefix=$HOME
make
make install


If you prefer to use other compilers, e.g:

./configure CC=icc CFLAGS="-qopenmp -xhost -O3"


Use of the library is straightforward. Simplest "example.c" C code:

#include <fermidirac.h>
#include <stdio.h>
#include <math.h>

int main()
{
  
  double result;
  
  
  result = Ffermi(4.0,1.0,1.0);
  
  printf("%.18lf\n",k,result);
  
  return 0;

}


Compile with:

gcc example.c -o example -lm -lfermidirac

Run:

./example

with expected reult:

114.066877991379016066


To use library in Mathematica, first compile MathLink code (replace 12.1 with your Mathematica $VersionNumber),
 which is located in doc subdirectory (Fermi-Dirac.tm and Fermi-Dirac.c files). 

/usr/local/Wolfram/Mathematica/12.1/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions/mprep Fermi-Dirac.tm -o Fermi-Dirac.tm.c

gcc -c Fermi-Dirac.tm.c Fermi-Dirac.c -I /usr/local/Wolfram/Mathematica/12.1/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions/

g++ Fermi-Dirac.tm.o Fermi-Dirac.o -o Fermi-Dirac -L /usr/local/Wolfram/Mathematica/12.1/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions/ -lML64i4 -lpthread -lrt -lstdc++ -ldl -luuid -lfermidirac

You might need to install libuuid developement files. 
On Ubuntu:
sudo apt install uuid-dev
On RedHat:
yum install libuuid-devel

Then, in Mathematica use:

SetDirectory["$HOME/libfermidirac-0.25/doc"]

Install["Fermi-Dirac"]

?FFermi

In:= FFermi[1.0, 1.0, 1.0]

Out = 2.64283

In:= GFermi[1.0, 1.0, 1.0]

Out = 58.676


Use in Fortran


Use in Python
