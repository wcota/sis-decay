! ## File: dynamics.f90
! ## - main program: running the SIS dynamics, based on OGA (Optimized Gillespie Algorithm).
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

program dynSIS
use mod_SIS_OGA
implicit none
    
    call print_info('################################################################################')
    call print_info('### Optimized Gillespie algorithms for the simulation of Markovian epidemic  ###')
    call print_info('############ processes on large and heterogeneous networks: SIS-OGA ############')
    call print_info('##============ Copyright (C) 2017 Wesley Cota, Silvio C. Ferreira ============##')
    call print_info('##===== Paper available at <http://dx.doi.org/10.1016/j.cpc.2017.06.007> =====##')
    call print_info('##======= The codes are available at <https://github.com/wcota/dynSIS> =======##')
    call print_info('##======== Please cite the above cited paper as reference to our code ========##')
    call print_info('##=== This code is under GNU General Public License. Please see README.md. ===##')
    call print_info('################################################################################')

    ! Files arguments read
    call read_arg(f_input)
    call read_arg(f_output)

    ! Read network to memory
    call net%readEdges(f_input)

    ! To be used in the SIS-OGA algorithm. Calculate the k_max of the network
    net_kmax = maxval(net%vertices(:)%degree)

    ! Initate the random generator
    call print_info('')
    call print_progress('Starting pseudo random number generator')
    call rgen%init(seed)
    call print_done()

    ! We are ready! All network data is here, now we need to read the dynamical parameters.
    call print_info('')
    call print_info('Now we need to read the dynamical parameters.')
    call read_dyn_parameters()

    ! Let's run the dynamics. But first, we allocate the average matrices
    allocate(avg_rho(dynp_tmax), avg_t(dynp_tmax), avg_sam(dynp_tmax), avg_samSurv(dynp_tmax))
    avg_rho = 0d0
    avg_t = 0d0
    avg_sam = 0
    avg_samSurv = 0
    dyn_dt_pos_max = 0
    call print_info('')
    call print_info('Running dynamics...')

    ! Loop over all the samples
    do dynp_i=1,dynp_sam
        write(f_temp,*) dynp_i
        call print_progress('Sample # '//trim(adjustl(f_temp)))

        ! Run dynamics
        call dyn_run()

        ! Open file and write info
        open(und_output,file=f_output)
        write(und_output,'(a)')         "## ***** Algorithm used: Optimized Gillespie Algorithm for SIS &
                                        & (SIS-OGA, Fortran) *****"
        write(und_output,'(a)')         "#@ Network file: "//trim(adjustl(f_input))
        write(und_output,'(a,i7)')      "#@ Number of nodes: ", net%N
        write(und_output,'(a,i7)')      "#@ Number of edges: ", net%skk
        write(und_output,'(a,i7)')      "#! Samples: ", dynp_i
        write(und_output,'(a,f11.5)')   "#! Infection rate lambda: ", dynp_lb
        write(und_output,'(a,i7)')      "#! Maximum time steps: ", dynp_tmax
        write(und_output,'(a,f11.5)')   "#! Fraction of infected vertices (initial condition): ", dynp_pINI

        do dyn_dt_pos = 1, dyn_dt_pos_max
            write(und_output,*) 1d0*avg_t(dyn_dt_pos)/avg_sam(dyn_dt_pos) , 1d0*avg_rho(dyn_dt_pos)/dynp_i
            ! If you use /avg_samSurv(dyn_dt_pos) instead of /dynp_i to write avg_rho (2nd column), you have 
            ! QS analysis :)
        enddo

        ! Close output file
        close(und_output)
        call print_done()
    enddo

    call print_info('')
    call print_info('Everything ok!')
    call print_info('Input file: '//trim(adjustl(f_input)))
    call print_info('Output file: '//trim(adjustl(f_output)))
    call print_info('')
    call print_info('*****Algorithm used: Optimized Gillespie Algorithm for SIS (SIS-OGA, Fortran)*****')
    call print_info('Codes available at <https://github.com/wcota/dynSIS>.')
    
end program
