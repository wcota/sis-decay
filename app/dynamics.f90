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

module mod_SIS_OGA
use netdata_mod
use read_tools_mod
use rndgen_mod
use kinds_mod
implicit none
    
    ! File read/write arguments and parameter
    character*1024                :: f_output, f_temp
    integer, parameter            :: und_output = 10
    
    
    ! Dynamics samples variables
    integer                       :: dynp_i, dynp_sam       ! samples vars
    real*8                        :: dynp_lb                ! lambda infection rate. mu is defined as = 1
    integer                       :: dynp_tmax              ! Maximum time steps
    real*8                        :: dynp_pINI              ! fraction of network first infected (random)
    
    ! Dynamics    
    real*8                        :: dyn_t, dyn_dt          ! Times and time step variables
    
    ! SIS-OGA - Dynamics Variables
    real*8                        :: dyn_m                  ! m = M/R. 1 - m = w = W/R
    real*8                        :: dyn_R                  ! Total rate
    integer, allocatable          :: dyn_VI(:), dyn_sig(:)  ! Lists V^I and sigma
    integer                       :: dyn_NI, dyn_Nk         ! # of infected vertex N_I and # of infected edges N_k
    
    ! SIS-OGA - Network structure variables
    integer                       :: net_kmax               ! Used in the rejection probability
    
    ! Output variables and average measures
    real*8, allocatable           :: avg_rho(:), avg_t(:)   ! Average for rho at times t, averaged
    integer, allocatable          :: avg_sam(:), avg_samSurv(:) ! # of samples for each time t, and of survivng ones
    integer                       :: dyn_dt_pos, dyn_dt_pos_max ! auxiliar vars

    ! random generator
    integer(kind=i4) :: seed = 20240915
    type(rndgen) :: rgen
    
contains
    
    subroutine read_dyn_parameters()
        call read_i(dynp_sam,"How much dynamics samples? ")
        call read_f(dynp_lb,"Value of infection rate lambda (mu is defined as equal to 1): ")
        call read_i(dynp_tmax,"Maximum time steps (it stops if the absorbing state is reached): ")
        call read_f(dynp_pINI,"Fraction of infected vertices on the network as initial condition (is random for &
 & each sample): ")
        
        ! Allocate the SIS-OGA lists V^I
        allocate(dyn_VI(net_N),dyn_sig(net_N))
    end subroutine
    
    subroutine random_initial_condition()
        integer :: ver, vti
        
        dyn_sig = 0 ! sigma
        !dyn_VI = 0 ! list V^I (not needed)
        dyn_NI = 0 ! N_I
        dyn_Nk = 0  ! N_k
        
        ! Sort vertices and apply the initial condition
        do vti = 1, int(net_N*dynp_pINI)
            vti_ver : do
                ver = rgen%int_i4(1,net_N)
                if (dyn_sig(ver) == 0) then
                    dyn_NI = dyn_NI + 1
                    dyn_VI(dyn_NI) = ver
                    dyn_sig(ver) = 1
                    dyn_Nk = dyn_Nk + net_k(ver)
                    exit vti_ver
                endif
            enddo vti_ver
        enddo
        
        dyn_t = 0d0
        dyn_dt_pos = 1
    end subroutine    

    subroutine dyn_run()
    
        integer :: pos_inf
        integer :: ver, pos_nei
        real*8  :: rnd
    
        call random_initial_condition()
        
        dyn_time_loop : do while (dyn_t <= dynp_tmax)
        
            rnd = rgen%rnd()

            ! Calculate the total rate
            dyn_R = (dyn_NI + 1d0*dynp_lb * dyn_Nk)
            
            ! Select the time step
            dyn_dt = -log(1.0_dp - rnd) / dyn_R
            
            ! Update the time
            dyn_t = dyn_t + dyn_dt
            
            ! Probability m to heal
            dyn_m = 1d0*dyn_NI/ dyn_R
            
            ! Try to heal
            rnd = rgen%rnd()
            if (rnd < dyn_m) then
                ! Select one infected vertex from V^I
                pos_inf = rgen%int_i4(1,dyn_NI) 
                ver = dyn_VI(pos_inf)
                
                ! Then, heal it
                dyn_sig(ver) = 0
                dyn_Nk = dyn_Nk - net_k(ver)
                dyn_VI(pos_inf) = dyn_VI(dyn_NI) ! Swap positions
                dyn_NI = dyn_NI - 1             ! Then, short the list
                
            ! If not, try to infect: w = 1 - m
            else
                ! Select the infected vertex i with prob. proportional to k_i
                select_infec : do
                    pos_inf = rgen%int_i4(1,dyn_NI)
                    ver = dyn_VI(pos_inf)
                    if (rgen%rnd() < 1d0*net_k(ver)/(1d0*net_kmax)) exit select_infec
                enddo select_infec
                
                ! select one of its neighbors
                pos_nei = rgen%int_i4(net_ini(ver) , net_ini(ver) + net_k(ver) - 1)
                ver = net_con(pos_nei)
                
                ! if not a phantom process, infect
                if (dyn_sig(ver) == 0) then
                    dyn_sig(ver) = 1
                    dyn_Nk = dyn_Nk + net_k(ver)
                    dyn_NI = dyn_NI + 1     ! Increase by 1 the list
                    dyn_VI(dyn_NI) = ver    ! Add one element to list
                endif
            endif
            
            ! Try to save the dynamics by time unit
            do while (dyn_t >= dyn_dt_pos)
                avg_rho(dyn_dt_pos) = avg_rho(dyn_dt_pos) + 1d0*dyn_NI/net_N
                avg_t(dyn_dt_pos) = avg_t(dyn_dt_pos) + dyn_t
                avg_sam(dyn_dt_pos) = avg_sam(dyn_dt_pos) + 1 
                if (dyn_NI .ne. 0) then
                    avg_samSurv(dyn_dt_pos) = avg_samSurv(dyn_dt_pos) + 1
                    dyn_dt_pos_max = max(dyn_dt_pos,dyn_dt_pos_max)         ! The maximum t with non-null rho
                endif
                dyn_dt_pos = dyn_dt_pos + 1
            enddo
            
            ! if a absorbing state is reached, exit
            if (dyn_NI == 0) exit dyn_time_loop
            
        enddo dyn_time_loop
        
    end subroutine

    subroutine run_prog()
        call print_info('################################################################################')
        call print_info('### Optimized Gillespie algorithms for the simulation of Markovian epidemic  ###')
        call print_info('############ processes on large and heterogeneous networks: SIS-OGA ############')
        call print_info('##============ Copyright (C) 2017 Wesley Cota, Silvio C. Ferreira ============##')
        call print_info('##===== Paper available at <http://dx.doi.org/10.1016/j.cpc.2017.06.007> =====##')
        call print_info('##======= The codes are available at <https://github.com/wcota/dynSIS> =======##')
        call print_info('##======== Please cite the above cited paper as reference to our code ========##')
        call print_info('##=== This code is under GNU General Public License. Please see README.md. ===##')
        call print_info('################################################################################')
    
        ! initial value of input counter to read command arguments, used by mod_read_tools
        inp_pos = 1 
    
        ! Files arguments read
        call read_arg(f_input)
        call read_arg(f_output)
    
        ! Read network to memory
        call readEdges()
    
        ! To be used in the SIS-OGA algorithm. Calculate the k_max of the network
        net_kmax = maxval(net_k)
    
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
            write(und_output,'(a,i7)')      "#@ Number of nodes: ", net_N
            write(und_output,'(a,i7)')      "#@ Number of edges: ", net_skk
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
        
    end subroutine

end module

program dynSIS
use mod_SIS_OGA
implicit none
    
    call run_prog()
    
end program
