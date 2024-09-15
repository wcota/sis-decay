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

module netdata_mod
    use read_tools_mod
    implicit none
    private
    
        type :: vertex
            integer :: degree ! degree
            integer, allocatable :: nei(:)
        end type vertex
    
        type :: network
            ! total number of vertices and edges
            integer                     :: N, skk
    
            type(vertex), allocatable :: vertices(:)
        contains
            procedure :: build_from_edges
            procedure :: readEdges
        end type network
    
        public :: network, edges_from_file
        
    contains
        ! to be compatible with old calls
        subroutine readEdges(this, f_input)
            class(network), intent(inout) :: this
            character(len=*) :: f_input

            integer, allocatable :: edges(:,:)

            call build_from_edges(this, edges_from_file(f_input))
            
        end subroutine

        subroutine build_from_edges(this, edges)
            class(network), intent(inout) :: this
            integer, intent(in) :: edges(:,:)

            ! number of edges (undirected)
            integer                     :: n_edges ! number of lines
            
            ! Temporary integers
            integer, allocatable        :: aux_degree(:) ! to check if all neighbors are in the list
            
            ! auxiliar integers
            integer                     :: iaux
            integer                     :: ver, ver1, ver2

            ! The network size is the maximum values. Now, assuming all edges undirected, this%skk = 2*number of edges,
            ! twice the number of lines (see README.md to know how to save the data).

            n_edges = size(edges,2)

            this%N = maxval(edges)
            this%skk = 2*n_edges

            allocate(this%vertices(this%N))
            this%vertices(:)%degree = 0 ! degree is zero at the beginning
        
            ! Is REALLY everything ok? Let's check!
            if (count(edges == 0) > 0) call print_error('Please, verify your data. Something is not right! [count(edges == 0) > 0]')

            call print_progress('Calculating degrees...')

            ! We will use aux_degree to store the position of the next neighbor to be added
            allocate(aux_degree(this%N))
            aux_degree = 0 ! auxiliar degree is zero at the beginning
            
            do iaux=1,size(edges,2)
                ver1 = edges(1,iaux)
                ver2 = edges(2,iaux)
                
                ! degree is added to both vertices
                this%vertices(ver1)%degree = this%vertices(ver1)%degree + 1
                this%vertices(ver2)%degree = this%vertices(ver2)%degree + 1
            enddo
            
            ! Allocate the neighbors list
            do concurrent (ver=1:this%N)
                allocate(this%vertices(ver)%nei(this%vertices(ver)%degree))
                this%vertices(ver)%nei = 0
            end do

            if (sum(this%vertices(:)%degree) .ne. this%skk) call print_error('Please, CHECK! Sum of degrees IS NOT equal the number of connections found & 
    & in the file')
            call print_done()

            call print_progress('Builing network object')
            ! Everything OK! 
            ! Now we will build the REAL list of adjacency.
            
            do iaux=1,size(edges,2)
                ver1 = edges(1,iaux)
                ver2 = edges(2,iaux)
                
                ! Just to check...
                if (this%vertices(ver1)%nei(aux_degree(ver1)+1) .ne. 0) call print_error('Somethings is wrong, con list overflow. Verify your data! &
    & [this%con(aux_ini(ver1)) .ne. 0]')
                if (this%vertices(ver2)%nei(aux_degree(ver2)+1) .ne. 0) call print_error('Somethings is wrong, con list overflow. Verify your data! &
    & [this%con(aux_ini(ver2)) .ne. 0]')
                
                ! Ok...
                this%vertices(ver1)%nei(aux_degree(ver1)+1) = ver2 
                aux_degree(ver1) = aux_degree(ver1) + 1

                this%vertices(ver2)%nei(aux_degree(ver2)+1) = ver1
                aux_degree(ver2) = aux_degree(ver2) + 1
            enddo
            
            ! Let's check again...
            if (any(this%vertices(:)%degree /= aux_degree)) call print_error('Verify your data. [count(this%con == 0) > 0]')
        
            ! We do not need edges and aux_degree anymore!
            deallocate(aux_degree)
            
            ! Now, everything is fine and we have the network data ready to be used! ;)
            call print_done()
            
        end subroutine

        function edges_from_file(f_input) result(edges)
            integer, allocatable        :: edges(:,:)
        
            ! input file read variables
            integer                     :: f_input_nl ! number of lines
            
            ! Temporary integers
            integer                     :: pos_con ! position of edges
            integer, allocatable        :: aux_degree(:) ! to check if all neighbors are in the list
            
            ! auxiliar integers
            integer                     :: aux_i1, aux_i2
            integer                     :: iaux
            integer                     :: ver, ver1, ver2
    
            ! file io
            character(len=*)            :: f_input
            integer                     :: f_sig ! file signal
            integer                     :: f_unit

            ! input file is opened
            open(newunit=f_unit,file=f_input,action='read',status='old')
        
            call print_progress('Calculating number of lines to read and allocating matrices')
            f_input_nl = 0
            f_input_loop : do
                read(f_unit,*,IOSTAT=f_sig) aux_i1, aux_i2
                if (f_sig < 0) exit f_input_loop ! if end of file, exit. No "goto" allowed! :)
                
                if (aux_i1 == aux_i2) call print_error("Self-connection found! Verify your data.")
                if (aux_i1 < 1 .or. aux_i2 < 1) call print_error('Vertex id MUST be >= 1. Verify your data!')
                
                f_input_nl = f_input_nl + 1
            enddo f_input_loop
        
            close(f_unit)
            open(newunit=f_unit,file=f_input,action='read',status='old') ! input file is opened again
            
            allocate(edges(2,f_input_nl)) ! We need to calculate the degree of the network, and to do this we will have 
                                            ! to read the file again. So, we save the data on the way.
                    
            ! number of edges added
            pos_con = 0
            edges = 0 ! list of connections read
            call print_done()
        
            ! Let's read the file again, this time saving the connections, since we already know the total number of connections.
            call print_progress('Reading the file again and collecting data')
            input_con_loop : do
                read(f_unit,*,IOSTAT=f_sig) aux_i1, aux_i2
                if (f_sig < 0) exit input_con_loop ! if end of file, exit.
    
                pos_con = pos_con + 1 ! new connection is saved
            
                ! save line to array
                edges(1,pos_con) = aux_i1
                edges(2,pos_con) = aux_i2
            enddo input_con_loop
            close(f_unit) ! we do not need the file anymore
            call print_done()

        end function edges_from_file
        
    end module