#    -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([libfermidirac],[0.27],[andrzej.odrzywolek@uj.edu.pl])
AC_CONFIG_SRCDIR([fermidirac.h])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror gnu])
AC_CONFIG_MACRO_DIR([m4])
        
AC_CONFIG_FILES([
Makefile
src/Makefile
doc/Makefile
examples/Makefile
])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AM_PROG_AR
LT_INIT
LT_INIT

# Check for presence of pdfLaTeX
AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
if test -z "$PDFLATEX"; then
  AC_MSG_WARN([Unable to create PDF version of the user manual.])
fi

AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt, nan, tgamma])
AC_CHECK_LIB([quadmath], [sqrtq, tgammaq])
# Checks for library functions.
AC_CHECK_FUNCS([sqrt nan tgamma])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdio.h stdlib.h string.h unistd.h math.h float.h quadmath.h])
AC_CHECK_HEADER([acb_calc.h], [], [AC_MSG_ERROR("Arb integrator missing! Install FLINT and Arb:\nsudo apt install libflint-arb-dev\nsudo apt install libflint-arb2\nsudo apt install libflint-dev")], [])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_HEADER_STDBOOL


AC_OUTPUT
