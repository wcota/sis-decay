! ## File: mod_read_tools.f90
! ## - module: Subroutines used to read parameters. This is just a module to be used in another program.
! ## See README.md for more information and use
!-----------------------------------------------------------------------------
! SIS epidemic model algorithm based on the article "Simulation of Markovian epidemic models on large networks"
! Copyright (C) 2016 Wesley Cota, Silvio C. Ferreira
! 
! Please cite the above cited paper as reference to our code.
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
! Email     : wesley.cota@ufv.br
! Date      : July 2016
!-----------------------------------------------------------------------------
! See README.md for more details
! This code is available at <https://github.com/wcota/dynSIS>

module mod_read_tools
    integer :: inp_pos
    
contains

    subroutine read_i(i1,a1)
        integer :: i1
        character(len=*) :: a1
        
        write(0,'(a)') 'read: '//trim(adjustl(a1))
        read(*,*) i1;
        write(*,'(a,i0.9)') 'r# '//trim(adjustl(a1))//' = ', i1 
    end subroutine    

    subroutine read_l(l1,a1)
        logical :: l1
        integer :: tmp
        character(len=*) :: a1
        
        write(0,'(a)') 'read: '//trim(adjustl(a1))
        read(*,*) tmp
        l1 = .false.
        if (tmp == 1) l1 = .true.
        write(*,'(a,l)') 'r# '//trim(adjustl(a1))//' = ', l1
    end subroutine
    
    subroutine read_f(f1,a1)
        real*8 :: f1
        character(len=*) :: a1
        
        write(0,'(a)') 'read: '//trim(adjustl(a1))
        read(*,*) f1;
        write(*,'(a,f7.4)') 'r# '//trim(adjustl(a1))//' = ', f1 
    end subroutine
    
    subroutine read_a(a1,a2)
        real*8 :: f1
        character(len=*) :: a1, a2
        
        write(0,'(a)') 'read: '//trim(adjustl(a2))
        read(*,*) a1
        write(*,'(a,a)') 'r# '//trim(adjustl(a2))//' = ', trim(adjustl(a1))
    end subroutine
    
    subroutine read_arg(a1)
        character(len=*) :: a1
        
        call getarg(inp_pos,a1)
        inp_pos = inp_pos + 1
    end subroutine
    
end module