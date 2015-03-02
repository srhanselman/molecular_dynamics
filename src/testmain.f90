! PHY 480 Project 2
! 2/3/15

Program lattice

use PhysicalProperties
use ObtainStep
use ForceCalculator

  Implicit None
  real(8):: temp, tempCalc, tempPre
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
  lattice_constant = dsqrt(2d0)*rm    !Length is in rm
  maxForceDistance = 5*rm
  dR = 1d-1*rm
  timeStep = 0.0004 !Time is in sqrt(argon_mass rm²/bondingEnergy) = 2.41980663d-12 s
  cell_dim = 6 !This is the number of lattice spacings within the box;
               !the number of particles (N) is 4/uc times cell_dim translational copies

  
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
 ! do counter=1,N
  !   print *,counter, pos(counter,1),pos(counter,2),pos(counter,3)
  !end do
  call momentum_initializer(temp,N,mass,momenta)
  call momentum_balancer(momenta,N)
  
  do counter=1,N
    particleKineticEnergy(counter) = dot_product(momenta(counter,:),momenta(counter,:))/(2*mass)
  end do
  print *, sum(particleKineticEnergy)*mass/N
  meanMomentumSq = dot_product(sum(momenta,1),sum(momenta,1))
  print *, meanMomentumSq
  call calculate_temperature(particleKineticEnergy,meanMomentumSq,mass,tempCalc)

!  do counter=1,N
!     print *, momenta(counter,1),";",momenta(counter,2),";",momenta(counter,3)
!  end do


  call force_potential_calculator_correlation(rm,bondingEnergy,maxForceDistance,box_length,pos,force,potential,dR,correlation)
  print *, "Start of simulation!"
  print *, "Initial temperature: ", tempCalc
  print *, "Step;KineticEnergy", "PotentialEnergy", "TotalEnergy", "CalculatedTemperature"
  do counter=1,1000
    
    call verlet_eqs_of_motion_correlation(pos,momenta,force,box_length,maxForceDistance,rm,bondingEnergy,mass,timeStep, &
    meanMomentumSq, particleKineticEnergy,particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,dR,correlation)
    call calculate_temperature(particleKineticEnergy,meanMomentumSq,mass,tempCalc)
    call momentum_balancer(momenta,N)
    print *, counter, kineticEnergyAfterStep, potentialEnergyAfterStep, &
 & kineticEnergyAfterStep+potentialEnergyAfterStep, tempCalc, meanMomentumSq
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
             + dt(3) * 24 * 60 * 60 * 1000 &
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
  integer :: counter
  integer, intent(in) :: N
  real(8) :: p_x, p_y, p_z, pxvar(N), pyvar(N), pzvar(N), conversion
  real(8), intent(out), dimension(N,3) :: momentum

  kb = 1.3806488d-2/1.65d0
  p_x = 0d0
  p_y = 0d0
  p_z = 0d0
  stand_dev = dsqrt(kb*mass*temp)
  !---------------------------------------------------------------
  !initializing velocities of particles according to a Gaussian distribution
  do counter =1,N
     call get_gaussian_momentum(mass,1d0,momentum(counter,1))
     call get_gaussian_momentum(mass,1d0,momentum(counter,2))
     call get_gaussian_momentum(mass,1d0,momentum(counter,3))
  end do
  !---------------------------------------------------------------
  !SRH: applying the standard deviation to the momenta - for some reason
  !the Gaussian RNG algorithm does not give any useful standard deviations!
  !---------------------------------------------------------------
  !finds total momentum in each direction and reduces to get closer to 0
  
  do counter=1,N
     print *,momentum(counter,1), "*", momentum(counter,2), "*", momentum(counter,3)
     p_x = p_x + momentum(counter,1)
     p_y = p_y + momentum(counter,2)
     p_z = p_z + momentum(counter,3)
  end do
  momentum(:,1) = momentum(:,1) - p_x/N
  momentum(:,2) = momentum(:,2) - p_y/N
  momentum(:,3) = momentum(:,3) - p_z/N

  pxvar(:) = momentum(:,1)*momentum(:,1)
  pyvar(:) = momentum(:,2)*momentum(:,2)
  pzvar(:) = momentum(:,3)*momentum(:,3)
  conversion = (sum(pxvar)+sum(pyvar)+sum(pzvar))/N
  print *, "conversion: ", conversion
  conversion = sqrt(conversion)/stand_dev
  print *, "conversion: ", conversion
  momentum = momentum*conversion

  do counter=1,N
     p_x = p_x + momentum(counter,1)
     p_y = p_y + momentum(counter,2)
     p_z = p_z + momentum(counter,3)
  end do


  print *,p_x,p_y,p_z

  !---------------------------------------------------------------
  return
end subroutine momentum_initializer

subroutine momentum_balancer(p,N)

  integer, intent(in) :: N
  real(8), intent(inout) :: p(N,3)
  real(8)                :: meanP(3)
  integer(8)             :: i

  meanP = sum(p,1)
  meanP = meanP/size(p,1)
  do i=1,size(p,1)
      p(i,:) = p(i,:) - meanP
  end do

end subroutine momentum_balancer
  
	

subroutine get_gaussian_momentum(mass,stand_dev,p)

  real(8), intent(in) :: stand_dev, mass
  real(8), intent(out) :: p


  real(8) :: random_1,random_y, I, p_integ, dp, pi, prefactor, y, test, povers, dI
  integer :: counter, check


  check = 0
  pi = 4d0*atan(1.0d0)
  dp = stand_dev/10d0
  p= 0d0
  prefactor = 1d0/(stand_dev)*dsqrt(2d0*pi)
  do while (check == 0)
     call random_number(random_1)
     random_1 = random_1
     p_integ = -10d0*stand_dev
     I = 0d0
     do while (p_integ < 10d0*stand_dev)
!	The cumulative Gaussian contributions appear to be too shaky,
!	so instead a (somewhat more expensive) direct erf calculation.
!	P_cumul = (1+erf(p/(sqrt(2)s)))/2
       if (p_integ > -3d0*stand_dev .AND. p_integ < -3d0*stand_dev) then
         povers = p_integ/stand_dev
         dI = povers/dsqrt(2d0)
         povers = -povers*povers
         counter = 1
         do while (abs(dI) > 1d-29)
	  dI = dI*povers/counter
          print *,dI
          I = I + dI/(2*counter+1)
	  print *,dI/(2*counter+1),I
          counter = counter + 1
        end do
      else
         I = I + prefactor*exp(-p_integ*p_integ/(2d0*stand_dev*stand_dev))*dp
         if (I > random_1) then
            p = p_integ
            exit
         else
            p_integ = p_integ + dp
         end if
     end if
     end do
     call random_number(random_y)
     random_y = 10d0*random_y
     test = prefactor*exp(-p*p/(2d0*stand_dev*stand_dev))
     if (random_y < test) then
        return
        check = 1
     end if
  end do
end subroutine get_gaussian_momentum
