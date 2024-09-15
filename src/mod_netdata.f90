! ## File: mod_netdata.f90
! ## - module: network and data read. This is just a module to be used in another program.
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

module mod_netdata
use mod_read_tools
implicit none

    ! file input names and signals
    character*1024              :: f_input
    integer                     :: f_sig ! file read signal
    
    ! total number of vertices and edges
    integer                     :: net_N, net_skk
    
    ! network data and adjacency list
    integer, allocatable        :: net_con(:), net_ini(:), net_k(:) 
    
contains
    subroutine readEdges()
    
        ! input file read variables
        integer                     :: f_input_nl ! number of lines
        
        ! Temporary integers
        integer, allocatable        :: tmp_con(:,:) ! temporary edges lists
        integer                     :: pos_con ! position of edges
        integer, allocatable        :: aux_ini(:)
        
        ! auxiliar integers
        integer                     :: aux_i1, aux_i2
        integer                     :: iaux, itmp
        integer                     :: ver, ver1, ver2
            
        ! input file is opened
        open(1,file=f_input,action='read',status='old')
        
        ! We set net_N to be the largest value read
        net_N = 0
    
        call print_progress('Calculating number of lines to read and allocating matrices')
        f_input_nl = 0
        f_input_loop : do
            read(1,*,IOSTAT=f_sig) aux_i1, aux_i2
            if (f_sig < 0) exit f_input_loop ! if end of file, exit. No "goto" allowed! :)
            
            if (aux_i1 == aux_i2) call print_error("Self-connection found! Verify your data.")
            if (aux_i1 < 1 .or. aux_i2 < 1) call print_error('Vertex id MUST be >= 1. Verify your data!')
            
            net_N = max(net_N, aux_i1, aux_i2) ! will be the largest value found
            
            f_input_nl = f_input_nl + 1
        enddo f_input_loop
    
        close(1); open(1,file=f_input,action='read',status='old') ! input file is opened again

        ! We already have the size of the network net_N. Now, assuming all edges undirected, net_skk = 2*f_input_nl,
        ! twice the number of lines (see README.md to know how to save the data).
        net_skk = 2*f_input_nl
    
        allocate(tmp_con(2,f_input_nl)) ! We need to calculate the degree of the network, and to do this we will have 
                                        ! to read the file again. So, we save the data on the way.
        allocate(net_k(net_N))
    
        ! number of edges added
        pos_con = 0
        tmp_con = 0 ! list of connections read
        net_k = 0 ! degree is zero at the beginning
        call print_done()
    
        ! Let's read the file again, this time saving the connections, since we already know the total number of connections.
        call print_progress('Reading the file again and collecting data')
        input_con_loop : do
            read(1,*,IOSTAT=f_sig) aux_i1, aux_i2
            if (f_sig < 0) exit input_con_loop ! if end of file, exit.

            pos_con = pos_con + 1 ! new connection is saved
        
            ! save line to array
            tmp_con(1,pos_con) = aux_i1
            tmp_con(2,pos_con) = aux_i2
        
            ! degree is added to both vertices
            net_k(aux_i1) = net_k(aux_i1) + 1
            net_k(aux_i2) = net_k(aux_i2) + 1
        enddo input_con_loop
        close(1) ! we do not need the file anymore
    
        ! Is REALLY everything ok? Let's check!
        if (count(tmp_con == 0) > 0) call print_error('Please, verify your data. Something is not right! [count(tmp_con == 0) > 0]')
        if (sum(net_k) .ne. net_skk) call print_error('Please, CHECK! Sum of degrees IS NOT equal the number of connections found & 
& in the file')
        call print_done()
        
        call print_progress('Builing adjacency list')
        ! Everything OK! 
        ! Now we will build the REAL list of adjacency.
        call make_ini() ! This will use the degrees array and build the net_ini matrix.
        allocate(aux_ini(net_N)) ! We will use to add the valid connections. It tell us the last *free* position to add to the net_con list.
    
        ! To begin, we set the initial values.
        aux_ini = net_ini
    
        ! We allocate the real list of adjacency.
        allocate(net_con(net_skk))
    
        ! Now we read again the data and add the connections
        net_con = 0 ! We do not have any connection added.
        do iaux=1,f_input_nl
            ver1 = tmp_con(1,iaux)
            ver2 = tmp_con(2,iaux)
            
            ! Just to check...
            if (net_con(aux_ini(ver1)) .ne. 0) call print_error('Somethings is wrong, con list overflow. Verify your data! &
& [net_con(aux_ini(ver1)) .ne. 0]')
            if (net_con(aux_ini(ver2)) .ne. 0) call print_error('Somethings is wrong, con list overflow. Verify your data! &
& [net_con(aux_ini(ver2)) .ne. 0]')
            
            ! Ok...
            net_con(aux_ini(ver1)) = ver2
            aux_ini(ver1) = aux_ini(ver1) + 1
            
            net_con(aux_ini(ver2)) = ver1
            aux_ini(ver2) = aux_ini(ver2) + 1
        enddo
        
        ! Let's check again...
        if (count(net_con == 0) > 0) call print_error('Verify your data. [count(net_con == 0) > 0]')
    
        ! We do not need tmp_con anymore!
        deallocate(tmp_con)
        
        ! Now, everything is fine and we have the network data ready to be used! ;)
        call print_done()
    
    contains

        subroutine make_ini() 
            
            call deal(net_ini)
            allocate(net_ini(net_N))
            net_ini = 0
    
            iaux = 1
            do ver=1,net_N
                itmp = net_k(ver)
                net_ini(ver) = iaux
                iaux = iaux + itmp
            enddo
        end subroutine 
    
    end subroutine readEdges
    
end module