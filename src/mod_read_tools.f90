! ## File: mod_read_tools.f90
! ## - module: Subroutines used to read parameters. This is just a module to be used in another program.
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

module mod_read_tools
    integer                        :: inp_pos
    character*3, private, parameter :: sub_read = '#? '
    
contains

    subroutine read_i(i1,a1)
        integer :: i1
        character(len=*) :: a1
        
        write(*,'(a)') sub_read//trim(adjustl(a1))
        read(*,*) i1;
        write(0,'(a,i0.9)') 'r# '//trim(adjustl(a1))//' = ', i1 
    end subroutine    

    subroutine read_l(l1,a1)
        logical :: l1
        integer :: tmp
        character(len=*) :: a1
        
        write(*,'(a)') sub_read//trim(adjustl(a1))
        read(*,*) tmp
        l1 = .false.
        if (tmp == 1) l1 = .true.
        write(0,'(a,l)') 'r# '//trim(adjustl(a1))//' = ', l1
    end subroutine
    
    subroutine read_f(f1,a1)
        real*8 :: f1
        character(len=*) :: a1
        
        write(*,'(a)') sub_read//trim(adjustl(a1))
        read(*,*) f1;
        write(0,'(a,f7.4)') 'r# '//trim(adjustl(a1))//' = ', f1 
    end subroutine
    
    subroutine read_a(a1,a2)
        character(len=*) :: a1, a2
        
        write(*,'(a)') sub_read//trim(adjustl(a2))
        read(*,*) a1
        write(0,'(a,a)') 'r# '//trim(adjustl(a2))//' = ', trim(adjustl(a1))
    end subroutine
    
    subroutine read_arg(a1)
        character(len=*) :: a1
        
        call getarg(inp_pos,a1)
        inp_pos = inp_pos + 1
    end subroutine

    subroutine print_warning(c1)
        character(*) :: c1
        
        write(0,'(a,a)') "##! Alert !## ", c1
    end subroutine

    subroutine print_error(c1)
        character(*) :: c1
        
        write(0,'(a,a)') "##!!! ERROR !!!## ", c1
        stop ""
    end subroutine

    subroutine print_info(c1)
        character(*)           :: c1
        
        write(*,'(a)') c1
    end subroutine

    subroutine print_progress(c1)
        character(*)           :: c1
        
        write(0,'(a)') c1//'... '
    end subroutine

    subroutine print_done()
        
        write(0,'(a)') char(9)//"done."
    end subroutine
    
    subroutine deal(arr)
        integer,allocatable :: arr(:)
        
        if (allocated(arr)) deallocate(arr)
    end subroutine
    
end module