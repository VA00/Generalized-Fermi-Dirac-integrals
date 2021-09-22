! Compile with gfortran example.f90 -o example -lfermidirac
      implicit none
! interface for C library http://th.if.uj.edu.pl/~odrzywolek/homepage/codes/libfermidirac/libfermidirac.html
! See https://stackoverflow.com/tags/fortran-iso-c-binding/info
      interface
        real(C_DOUBLE) function ofermi(k, eta, theta) bind(C, name='Ffermi')
          use ISO_C_BINDING
          implicit none
          real(C_DOUBLE), value :: k, eta, theta
        end function ofermi
      end interface

      real*8   :: f
! evaluate the fermi function for k=4, eta=1, theta=1
      f = ofermi(4.0d0, 1.0d0, 1.0d0)
	  write(*,*) f
	  
! Expected reference result
! Ffermi(4,1,1)=114.06687799137902589992508824775903265360960251309372394266372713106403373794999046746440912671766590817078780479429221800302036

end program example


