! PHY 480 Project 2
! 2/3/15

Program lattice

use PhysicalProperties
use ObtainStep
use ForceCalculator

  Implicit None
  real(8):: temp
  real(8):: mass
  integer ::cell_dim,N,counter
  real(8):: lattice_constant, box_length, maxForceDistance, rm, bondingEnergy
  real(8):: timeStep, meanMomentumSq, kineticEnergyAfterStep, dR
  real(8):: potentialEnergyAfterStep
  real(8), allocatable :: pos(:,:)
  real(8), allocatable :: momenta(:,:)
  real(8), allocatable :: force(:,:)
  real(8), allocatable :: potential(:)
  real(8), allocatable :: particleKineticEnergy(:)
  real(8), allocatable :: particlePotential(:)
  real(8), allocatable :: correlation(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A short explanation of the units: the goal is to have rm = mass = bondingEnergy
!   = 1 (i.e. unitless). This means:
! kg -> argon_mass (6.6335209d-26 kg);   m -> rm (3.81637096425d-10 m);
! J = kg m²/s² -> bondingEnergy (1.65d-21 J);
! s = sqrt(kg m²/J) -> sqrt(argon_mass rm²/bondingEnergy)
! 
! As temperature is not constant and there is no reference temperature, K survives
! as a unit but is accompanied by kB, which changes:
!
! kB = 1 J/K = 1 (bondingEnergy/J)^-1 bondingEnergy/K = 1/1.65d-21 K^-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  temp = 1d0		!Temperature is in K (apart from kB, this is the sole exception)

  mass = 1d0		!Mass is in argon_mass
  rm = 1d0		!Length is in rm
  bondingEnergy = 1d0   !Energy is in bondingEnergy
  lattice_constant = 1*rm    !Length is in rm
  maxForceDistance = 5*rm
  dR = 1d-1*rm
  timeStep = 0.0004 !Time is in sqrt(argon_mass rm²/bondingEnergy) = 2.41980663d-12 s
  cell_dim = 5 !This is the number of lattice spacings within the box

  
  N = 4*cell_dim**3
  box_length = cell_dim*lattice_constant
  allocate( pos(N,3) )
  allocate( momenta(N,3) )
  allocate( force(N,3) )
  allocate( potential(N) )
  allocate( particleKineticEnergy(N) )
  allocate( particlePotential(N) )
  allocate( correlation(ceiling(dsqrt(3.0d0)*box_length/(2*dR))+1) )

  call init_random_seed()
  call position_initializer(N,cell_dim,pos)
  pos = pos*lattice_constant
  do counter=1,N
     print *,counter, pos(counter,1),pos(counter,2),pos(counter,3)
  end do
  call momentum_initializer(temp,N,mass,momenta)
  do counter=1,N
     print *,counter, momenta(counter,1),momenta(counter,2),momenta(counter,3)
  end do

  call force_potential_calculator_correlation(rm,bondingEnergy,maxForceDistance,box_length,pos,force,potential,dR,correlation)
  do counter=1,100
    call obtain_step_correlation(pos,momenta,force,box_length,maxForceDistance,rm,bondingEnergy,mass,timeStep,meanMomentumSq, &
    particleKineticEnergy,particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,dR,correlation)
    print *, counter, kineticEnergyAfterStep, potentialEnergyAfterStep
  end do

end program

subroutine init_random_seed()
  implicit none

  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid, t(2), s
  integer(8) :: count, tms
  call random_seed(size = n)
  allocate(seed(n))
  open(newunit=un, file="/dev/urandom", access="stream",&
       form="unformatted", action="read", status="old", &
       iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     call system_clock(count)
     if (count /= 0) then
        t = transfer(count, t)
     else
        call date_and_time(values=dt)
        tms = (dt(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24 * 60 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
        t = transfer(tms, t)
     end if
     s = ieor(t(1), t(2))
     pid = getpid() + 1099279 ! Add a prime
     s = ieor(s, pid)
     if (n >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (n > 3) then
           seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
     else
        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
     end if
  end if
  call random_seed(put=seed)
end subroutine init_random_seed

subroutine position_initializer(N,dim_cells,position)
  integer :: N,dim_cells
  real(8),intent(out), dimension(N,3) :: position
  integer :: counter1,x_counter,y_counter,z_counter,particle_counter
  real(8) :: cell_width
  cell_width = 1d0
  !--------------------------------------------------------------
  !initializing positions of particles in room as FCC lattice
  particle_counter =1
  do z_counter = 0, dim_cells-1
     do y_counter =0,dim_cells-1
        do x_counter = 0, dim_cells-1
           position(particle_counter,:) = (/ 0d0+ x_counter,0d0+ y_counter,0d0 + z_counter /)
           position(particle_counter + 1,:) = (/ 0.5d0 + x_counter, 0.5d0 + y_counter,0d0+ z_counter /)
           position(particle_counter + 2,:) = (/ 0.5d0 + x_counter,0d0+ y_counter, 0.5d0 + z_counter /)
           if (particle_counter+3 > N) then
              print *, "Out of bounds in position array"
              exit
           end if
           position(particle_counter + 3,:) = (/ 0d0+x_counter, 0.5d0 + y_counter, 0.5d0 + z_counter /)
           particle_counter = particle_counter + 4
        end do
     end do
  end do
  !---------------------------------------------------------------
  return
end subroutine position_initializer

subroutine momentum_initializer(temp, N, mass,momentum)
  real(8) :: temp, mass, stand_dev, kb
  integer :: N, counter
  real(8) :: p_x, p_y, p_z
  real(8), intent(out), dimension(N,3) :: momentum
  kb = 1/1.65d-21
  p_x = 0d0
  p_y = 0d0
  p_z = 0d0
  stand_dev = dsqrt(kb*mass*temp)
  !---------------------------------------------------------------
  !initializing velocities of particles according to a Gaussian distribution
  do counter =1,N
     call get_gaussian_momentum(mass,stand_dev,momentum(counter,1))
     call get_gaussian_momentum(mass,stand_dev,momentum(counter,2))
     call get_gaussian_momentum(mass,stand_dev,momentum(counter,3))
  end do
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  !finds total momentum in each direction and reduces to get closer to 0
  do counter=1,N
     p_x = p_x + momentum(counter,1)
     p_y = p_y + momentum(counter,2)
     p_z = p_z + momentum(counter,3)
  end do
  momentum(:,1) = momentum(:,1) - p_x/N
  momentum(:,2) = momentum(:,2) - p_y/N
  momentum(:,3) = momentum(:,3) - p_z/N
  p_x =0d0
  p_y = 0d0
  p_z = 0d0
  do counter=1,N
     p_x = p_x + momentum(counter,1)
     p_y = p_y + momentum(counter,2)
     p_z = p_z + momentum(counter,3)
  end do
  print *,p_x,p_y,p_z
  !---------------------------------------------------------------
  return
end subroutine momentum_initializer

subroutine get_gaussian_momentum(mass,stand_dev,x)

  real(8), intent(in) :: stand_dev, mass
  real(8), intent(out) :: x


  real(8) :: random_1,random_y, I, x_integ, dx, pi, prefactor, y, test
  integer :: counter, check
  x = 3.0d0  
  check = 0
  pi = 4d0*atan(1.0d0)
  dx = stand_dev/10d0
  x= 0d0
  prefactor = 1d0/(stand_dev)*dsqrt(2d0*pi)
  do while (check == 0)
     call random_number(random_1)
     random_1 = random_1
     x_integ = -10d0*stand_dev
     I = 0d0
     do while (x_integ < 10d0*stand_dev)
        I = I + prefactor*exp(-x_integ*x_integ/(2d0*stand_dev*stand_dev))*dx
        if (I > random_1) then
           x = x_integ*mass
           exit
        else
           x_integ = x_integ + dx
        end if
     end do
     call random_number(random_y)
     random_y = 10d0*random_y
     test = prefactor*exp(-x*x/(2d0*stand_dev*stand_dev))
     if (random_y < test) then
        return
        check = 1
     end if
  end do
end subroutine get_gaussian_momentum
