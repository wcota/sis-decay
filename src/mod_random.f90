! ## File: mod_random.f90
! ## - module: random number generator. This is just a module to be used in another program.
! ## See README.md for more information and use
!-----------------------------------------------------------------------------
! SIS epidemic model algorithm based on the article
!           Computer Physics Communications 219C (2017) pp. 303-312
!           "Optimized Gillespie algorithms for the simulation of 
!            Markovian epidemic processes on large and heterogeneous networks"
! Copyright (C) 2017 Wesley Cota, Silvio C. Ferreira
! 
! Please cite the above cited paper (available at <http://dx.doi.org/10.1016/j.cpc.2017.06.007> ) 
! as reference to our code.
! 
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------
! Author    : Wesley Cota
! Email     : wesley@wcota.me
! Date      : 27 Mar 2017
! Version   : 1.0
!-----------------------------------------------------------------------------
! See README.md for more details
! This code is available at <https://github.com/wcota/dynSIS>
! For pure Python, see <https://github.com/wcota/dynSIS-py>
! For NetworkX library, see <https://github.com/wcota/dynSIS-networkx> (NetworkX implementation)
! 
! IMPORTANT:
! THIS CODE WAS MODIFIED FROM http://web.mst.edu/~vojtat/class_5403/kiss05/rkiss05.f90
! ! FORTRAN implementation by Thomas Vojta, vojta@mst.edu
! ! built on a module found at www.fortran.com

! Solution to ifort and gfortran
!DEC$ IF(.FALSE.)
module ifport
end module ifport
!DEC$ ENDIF

module mod_random
    
    integer :: iseed
    
contains

    function random_d() ! KISS
        implicit none

        integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
        integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 
        real(r8b),parameter    :: am=4.656612873077392578d-10       ! multiplier 1/2^31

        real(r8b)              :: random_d
        integer(i4b)           :: kiss
        integer(i4b)           :: x,y,z,w              ! working variables for the four generators
        common /kisscom/x,y,z,w 

        x = 69069 * x + 1327217885
        y = ieor (y, ishft (y, 13)); y= ieor (y, ishft (y, -17)); y= ieor (y, ishft (y, 5))
        z = 18000 * iand (z, 65535) + ishft (z, - 16)
        w = 30903 * iand (w, 65535) + ishft (w, - 16)
        kiss = ishft(x + y + ishft (z, 16) + w , -1)
        random_d=kiss*am
    end function
    
    function random_int(i1,i2)
        integer, intent(in) :: i1, i2
        integer             :: random_int
        
        random_int = min(int(random_d()*(i2+1-i1))+i1,i2)
    end function

    subroutine random_ini()
        use ifport
        implicit none
        integer,parameter     :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
        integer,parameter     :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

        integer(i4b)          :: idum,ia,im,iq,ir
        integer(i4b)          :: k,x,y,z,w,c1
        real(r8b)             :: rdum
        
        ! Generate iseed using time and pid values
        integer :: s, pid
        
        parameter (ia=16807,im=2147483647,iq=127773,ir=2836)
        common /kisscom/x,y,z,w
        
        ! Generate iseed using time and pid values
        ! Modified from <http://stackoverflow.com/questions/8920411/possible-sources-for-random-number-seeds>
        pid = getpid()
        call system_clock(s)
        iseed = abs( mod((s*181)*((pid-83)*359), 101248729) )
        ! \ end generate iseed

        !!! Test integer representation !!!
        c1=-8
        c1=ishftc(c1,-3)
        !     print *,c1
        if (c1.ne.536870911) stop 'Nonstandard integer representation. Stopped.'

        idum=iseed
        idum= abs(1099087573 * idum)               ! 32-bit LCG to shuffle seeds
        if (idum.eq.0) idum=1
        if (idum.ge.IM) idum=IM-1

        k=(idum)/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum = idum + IM
        if (idum.lt.1) then
                x=idum+1 
            else 
                x=idum
        endif
        k=(idum)/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum = idum + IM
        if (idum.lt.1) then 
                y=idum+1 
            else 
                y=idum
        endif
        k=(idum)/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum = idum + IM
        if (idum.lt.1) then
                z=idum+1 
        else 
                z=idum
        endif
        k=(idum)/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum = idum + IM
        if (idum.lt.1) then
                w=idum+1 
            else 
        w=idum
        endif

        rdum=random_d()

        return
    end subroutine
    
end module
