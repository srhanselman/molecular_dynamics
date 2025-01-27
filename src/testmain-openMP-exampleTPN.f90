! PHY 480 Project 2
! 2/3/15

Program lattice

use PhysicalProperties
use ObtainStep
use ForceCalculator
use md_plot

  Implicit None
  real(8):: temp, finalTemp, tempCalc, tempPre
  real(8):: mass
  integer ::cell_dim,N,counter,counterTemp, nTimeStepsEquilibration, &
 nTimeStepsRun, nTimeStepsTempStep, nTimeStepsTempStepRun, &
 nChemicalPotentialTemperatures, innerCounter,nThreads
  real(8):: lattice_constant, boxLength, maxForceDistance, rm, bondingEnergy
  real(8):: timeStep, meanMomentumSq, kineticEnergyAfterStep, dR, pressure
  real(8):: potentialEnergyAfterStep, tempStep, kB
  real(8), allocatable :: pos(:,:)
  real(8), allocatable :: momenta(:,:)
  real(8), allocatable :: force(:,:)
  real(8), allocatable :: potential(:)
  real(8), allocatable :: particleKineticEnergy(:)
  real(8), allocatable :: particlePotential(:)
  real(8), allocatable :: correlation(:)
  real(8), allocatable :: chemicalPotentialTemperatures(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The program can be controlled using these switches:
!  - correctTemperature contains either 0 or the temperature tolerance.
!  - constantPressure is 0 for an TVN system, or /=0 for a TPN system
!  - check_heat_capacity and check_chemical_potential switch on tests after
!      equilibration.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real:: correctTemperature, constantPressure
  logical:: checkHeatCapacity, checkChemicalPotential

  constantPressure = 1d5/29684667.9463 !Pressure is in eps/rm^3 = 29684667.9463 Pa
  correctTemperature = 1d-2
  checkHeatCapacity = .FALSE.
  checkChemicalPotential = .FALSE.

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
  kB = 1.3806488d-2/1.65d0
  temp = 5d0		!Temperature is in K (apart from kB, this is the sole exception)
  finalTemp = 1d2
  tempStep = 1d0
  mass = 1d0		!Mass is in argon_mass
  rm = 1d0		!Length is in rm
  bondingEnergy = 1d0   !Energy is in bondingEnergy
  lattice_constant = dsqrt(2d0)*rm    !Length is in rm
  maxForceDistance = 5*rm
  dR = 2d-2*rm
  timeStep = 0.0004 !Time is in sqrt(argon_mass rm²/bondingEnergy) = 2.41980663d-12 s
  cell_dim = 9 !This is the number of lattice spacings within the box;
               !the number of particles (N) is 4/uc times cell_dim translational copies
  nTimeStepsEquilibration = 1000
  nTimeStepsRun = 500
  nTimeStepsTempStep = 800
  nTimeStepsTempStepRun = 500
  
  N = 4*cell_dim**3
  boxLength = cell_dim*lattice_constant
  allocate( pos(N,3) )
  allocate( momenta(N,3) )
  allocate( force(N,3) )
  allocate( potential(N) )
  allocate( particleKineticEnergy(N) )
  allocate( particlePotential(N) )
  allocate( correlation(ceiling(dsqrt(3.0d0)*boxLength/(2*dR))+1) )

  nThreads = 2

!!!
!!! Opening the files needed to start:
!!!

  open(unit = 1, file="equilibration.csv")
  open(unit = 2, file="run.csv")
  if (checkHeatCapacity) then
    open(unit=10, file="heatcapacity.csv")
  end if
  if (checkChemicalPotential) then
    open(unit=11, file="chemicalpotential.csv")
  end if

!!!
!!! First, the system's particle positions and momenta will be initialized.
!!!

  call init_random_seed()
  call position_initializer(N,cell_dim,pos)
  pos = pos*lattice_constant
 
  call momentum_initializer(temp,N,mass,momenta)
  call momentum_balancer(momenta,N)
  
  do counter=1,N
    particleKineticEnergy(counter) = dot_product(momenta(counter,:),momenta(counter,:))/(2*mass)
  end do
  print *, "Initial total kinetic energy: ", sum(particleKineticEnergy)
  meanMomentumSq = dot_product(sum(momenta,1),sum(momenta,1))
  print *, "Initial mean momentum squared: ", meanMomentumSq
  write (1,*) "The system contains ", N, " particles."
  write (2,*) "The system contains ", N, " particles."
  call calculate_temperature(particleKineticEnergy,meanMomentumSq,mass,tempCalc)
  if(correctTemperature /= 0 .AND. correctTemperature < abs(tempCalc/temp-1)) then
    momenta = momenta * sqrt(temp/tempCalc)
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!   !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! START OF SIMULATION --- START OF SIMULATION --- START OF SIMULATION --- START OF SIMULATION --- START OF SIMULATION !                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!   !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (constantPressure == 0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The subsequent procedure is ONLY MEANT FOR TVN systems - it does not allow for boxLength resizing. !!!!!!!!!
!!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!   !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!
!! TVN ONLY --- TVN ONLY --- TVN ONLY --- TVN ONLY --- TVN  ONLY --- TVN ONLY --- TVN ONLY --- TVN ONLY --- TVN ONLY !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!
!!! The simulation starts by calculating the initial forces and starting the stepping procedure.
!!!

  call force_potential_calculator_correlation(rm,bondingEnergy,maxForceDistance,boxLength,pos,force,potential,dR,correlation)
  print *, "Start of simulation!"
  print *, "Initial temperature: ", tempCalc
  write (1,*) "Step; ", "Kinetic energy; ", "Potential energy; ", "Total energy; ", "Calculated temperature; ", &
 "Mean momentum squared; ", "Pressure; "
  write (2,*) "Step; ", "Kinetic energy; ", "Potential energy; ", "Total energy; ", "Calculated temperature; ", &
 "Mean momentum squared; ", "Pressure; "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   This is the initial procedure, in which the temperature may be equilibrated.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  do counter=1,nTimeStepsEquilibration
    
    call verlet_eqs_of_motion_correlation(pos,momenta,force,boxLength,maxForceDistance,rm,bondingEnergy,mass,timeStep, &
    meanMomentumSq, particleKineticEnergy,particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,dR,correlation)
    call calculate_temperature(particleKineticEnergy,meanMomentumSq,mass,tempCalc)
    call momentum_balancer(momenta,N)
    call calculate_pressure(pos,force,boxLength,N,tempCalc,pressure)
    write (1,*) counter,";",kineticEnergyAfterStep,";",potentialEnergyAfterStep, &
 ";",kineticEnergyAfterStep+potentialEnergyAfterStep,";",tempCalc,";",meanMomentumSq, &
 ";",pressure
    if(correctTemperature /= 0 .AND. correctTemperature < abs(tempCalc/temp)-1) then
      momenta = momenta * sqrt(temp/tempCalc)
    end if
  end do

correlation = 0

if(nTimeStepsTempStep == 0) then
  do counter=1,nTimeStepsRun
    call verlet_eqs_of_motion_correlation(pos,momenta,force,boxLength,maxForceDistance,rm,bondingEnergy,mass,timeStep, &
    meanMomentumSq, particleKineticEnergy,particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,dR,correlation)
    call calculate_temperature(particleKineticEnergy,meanMomentumSq,mass,tempCalc)
    call momentum_balancer(momenta,N)
    call calculate_pressure(pos,force,boxLength,N,tempCalc,pressure)
    write (2,*) counter,";",kineticEnergyAfterStep,";",potentialEnergyAfterStep, &
 ";",kineticEnergyAfterStep+potentialEnergyAfterStep,";",tempCalc,";",meanMomentumSq, &
 ";",pressure
  end do
  write (2,*) "Correlation function: ", correlation
else
 do while (temp < FinalTemp+5d-1)
!$omp parallel do schedule(static,1) shared(nTimeStepsTempStep,boxLength,maxForceDistance, &
!$omp&  rm,bondingEnergy,mass,timeStep,dR,N) private(meanMomentumSq,particleKineticEnergy, &
!$omp&  particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,correlation, &
!$omp&  tempCalc,pressure,innerCounter,counter) firstprivate(pos,momenta,force,temp) &
!$omp&  lastprivate(pos,momenta,force,temp)
   do innerCounter=1,nThreads
!!!!!!!!!!!!!! BEGIN OF EQUILIBRATION !!!!!!!!!!!!!
!$omp flush
    temp=temp+innerCounter*tempStep
    do counter=1,nTimeStepsTempStep
     call verlet_eqs_of_motion_correlation(pos,momenta,force,boxLength,maxForceDistance, &
      rm,bondingEnergy,mass,timeStep,meanMomentumSq,particleKineticEnergy, &
      particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,dR,correlation)
     call calculate_temperature(particleKineticEnergy,meanMomentumSq,mass,tempCalc)
     call momentum_balancer(momenta,N)
     call calculate_pressure(pos,force,boxLength,N,tempCalc,pressure)
!$omp CRITICAL
     write (1,*) counter,";",kineticEnergyAfterStep,";",potentialEnergyAfterStep, &
   ";",kineticEnergyAfterStep+potentialEnergyAfterStep,";",tempCalc,";",temp,";", &
   meanMomentumSq,";",pressure
!$omp END CRITICAL
     if(correctTemperature /= 0 .AND. correctTemperature < abs(tempCalc/temp)-1) then
       momenta = momenta * sqrt(temp/tempCalc)
     end if
    end do
!!!!!!!!!!!!!!!!! BEGIN OF RUN !!!!!!!!!!!!!!!!!!!!
    correlation = 0
   do counter=1,nTimeStepsTempStepRun
     call verlet_eqs_of_motion_correlation(pos,momenta,force,boxLength,maxForceDistance, &
      rm,bondingEnergy,mass,timeStep,meanMomentumSq,particleKineticEnergy, &
      particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,dR,correlation)
     call calculate_temperature(particleKineticEnergy,meanMomentumSq,mass,tempCalc)
     call momentum_balancer(momenta,N)
     call calculate_pressure(pos,force,boxLength,N,tempCalc,pressure)
!$omp CRITICAL
     write (2,*) counter,";",kineticEnergyAfterStep,";",potentialEnergyAfterStep, &
   ";",kineticEnergyAfterStep+potentialEnergyAfterStep,";",tempCalc,";",temp,";", &
   meanMomentumSq,";",pressure
!$omp END CRITICAL
   end do
!$omp CRITICAL
   write (2,*) "Correlation function at temperature ", temp, ": ", correlation
!$omp END CRITICAL
  end do
!$omp end parallel do
 end do
 end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else if (constantPressure > 0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The subsequent procedure is ONLY MEANT FOR TPN systems - it resizes boxLength before EACH time step !!!!!!!!!
!!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!   !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!  !!
!! TPN ONLY --- TPN ONLY --- TPN ONLY --- TPN ONLY --- TPN  ONLY --- TPN ONLY --- TPN ONLY --- TPN ONLY --- TPN ONLY !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!
!!! The simulation starts by calculating the initial forces and starting the stepping procedure.
!!! Note that the boxLength is first resized to fit the ideal gas law (NkBT/p = V);
!!! moreover, the correlation box spacing is decreased/increased correspondingly.
!!!
  dR = dR*((N*kB*temp)/constantPressure)**(1d0/3d0)/boxLength
  boxLength = ((N*kB*temp)/constantPressure)**(1d0/3d0)

  call force_potential_calculator_correlation(rm,bondingEnergy,maxForceDistance,boxLength,pos,force,potential,dR,correlation)
  print *, "Start of simulation!"
  print *, "Initial temperature: ", tempCalc
  write (1,*) "Initial dR is ", dR
  write (2,*) "Initial dR is ", dR

  write (1,*) "Step; ", "Kinetic energy; ", "Potential energy; ", "Total energy; ", "Calculated temperature; ", &
 "Mean momentum squared; ", "Pressure; ", "Box length; ", "Volume"
  write (2,*) "Step; ", "Kinetic energy; ", "Potential energy; ", "Total energy; ", "Calculated temperature; ", &
 "Mean momentum squared; ", "Pressure; ", "Box length; ", "Volume"




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   This is the initial procedure, in which the temperature may be equilibrated.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  do counter=1,nTimeStepsEquilibration

    call verlet_eqs_of_motion_correlation(pos,momenta,force,boxLength,maxForceDistance,rm,bondingEnergy,mass,timeStep, &
    meanMomentumSq, particleKineticEnergy,particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,dR,correlation)
    call calculate_temperature(particleKineticEnergy,meanMomentumSq,mass,tempCalc)
    call momentum_balancer(momenta,N)
    call calculate_pressure(pos,force,boxLength,N,tempCalc,pressure)
    write (1,*) counter,";",kineticEnergyAfterStep,";",potentialEnergyAfterStep, &
 ";",kineticEnergyAfterStep+potentialEnergyAfterStep,";",tempCalc,";",meanMomentumSq, &
 ";",pressure,";",boxLength,";",boxLength**(3d0)
    if(correctTemperature /= 0 .AND. correctTemperature < abs(tempCalc/temp)-1) then
      momenta = momenta * sqrt(temp/tempCalc)
    end if
  end do

correlation = 0

if(nTimeStepsTempStep == 0) then
  do counter=1,nTimeStepsRun
    !  Volume correction (PV)  !
    boxLength = boxLength*((pressure/constantPressure)**(1d0/3d0))
    pos = pos*((pressure/constantPressure)**(1d0/3d0))
    ! End of volume correction !
    call verlet_eqs_of_motion_correlation(pos,momenta,force,boxLength,maxForceDistance,rm,bondingEnergy,mass,timeStep, &
    meanMomentumSq, particleKineticEnergy,particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,dR,correlation)
    call calculate_temperature(particleKineticEnergy,meanMomentumSq,mass,tempCalc)
    call momentum_balancer(momenta,N)
    call calculate_pressure(pos,force,boxLength,N,tempCalc,pressure)
    write (2,*) counter,";",kineticEnergyAfterStep,";",potentialEnergyAfterStep, &
 ";",kineticEnergyAfterStep+potentialEnergyAfterStep,";",tempCalc,";",meanMomentumSq, &
 ";",pressure,";",boxLength,";",boxLength**(3d0)    
  end do
  write (2,*) "Correlation function: ", correlation
else
 do while (temp < FinalTemp+5d-1)
!$omp parallel do schedule(static,1) shared(nTimeStepsTempStep,maxForceDistance, &
!$omp&  rm,bondingEnergy,mass,timeStep,dR,N) private(meanMomentumSq,particleKineticEnergy, &
!$omp&  particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,correlation, &
!$omp&  tempCalc,pressure,innerCounter,counter) firstprivate(pos,momenta,force,temp,boxLength) &
!$omp&  lastprivate(pos,momenta,force,temp,boxLength)
   do innerCounter=1,nThreads
!!!!!!!!!!!!!! BEGIN OF EQUILIBRATION !!!!!!!!!!!!!
!$omp flush
    boxLength = boxLength*(((temp+innerCounter*tempStep)/temp)**(1d0/3d0))
    temp=temp+innerCounter*tempStep
    
    do counter=1,nTimeStepsTempStep
     !  Volume correction (PV)  !
     boxLength = boxLength*((pressure/constantPressure)**(1d0/3d0))
     pos = pos*((pressure/constantPressure)**(1d0/3d0))
     ! End of volume correction !
     call verlet_eqs_of_motion_correlation(pos,momenta,force,boxLength,maxForceDistance, &
      rm,bondingEnergy,mass,timeStep,meanMomentumSq,particleKineticEnergy, &
      particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,dR,correlation)
     call calculate_temperature(particleKineticEnergy,meanMomentumSq,mass,tempCalc)
     call momentum_balancer(momenta,N)
     call calculate_pressure(pos,force,boxLength,N,tempCalc,pressure)
!$omp CRITICAL
    write (1,*) counter,";",kineticEnergyAfterStep,";",potentialEnergyAfterStep, &
 ";",kineticEnergyAfterStep+potentialEnergyAfterStep,";",tempCalc,";",meanMomentumSq, &
 ";",pressure,";",boxLength,";",boxLength**(3d0)
!$omp END CRITICAL
     if(correctTemperature /= 0 .AND. correctTemperature < abs(tempCalc/temp)-1) then
       momenta = momenta * sqrt(temp/tempCalc)
     end if
    end do
!!!!!!!!!!!!!!!!! BEGIN OF RUN !!!!!!!!!!!!!!!!!!!!
    correlation = 0
   do counter=1,nTimeStepsTempStepRun
     call verlet_eqs_of_motion_correlation(pos,momenta,force,boxLength,maxForceDistance, &
      rm,bondingEnergy,mass,timeStep,meanMomentumSq,particleKineticEnergy, &
      particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,dR,correlation)
     call calculate_temperature(particleKineticEnergy,meanMomentumSq,mass,tempCalc)
     call momentum_balancer(momenta,N)
     call calculate_pressure(pos,force,boxLength,N,tempCalc,pressure)
!$omp CRITICAL
    write (2,*) counter,";",kineticEnergyAfterStep,";",potentialEnergyAfterStep, &
 ";",kineticEnergyAfterStep+potentialEnergyAfterStep,";",tempCalc,";",meanMomentumSq, &
 ";",pressure,";",boxLength,";",boxLength**(3d0)
!$omp END CRITICAL
   end do
!$omp CRITICAL
   write (2,*) "Correlation function at temperature ", temp, ": ", correlation
!$omp END CRITICAL
  end do
!$omp end parallel do
 end do
 end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

else
  print *, "No negative pressures allowed!"

end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  print *, "Simulation completed!"

  close(1)
  close(2)
  if (checkHeatCapacity) then
    close(10)
  end if
  if (checkChemicalPotential) then
    close(11)
  end if


  deallocate( pos )
  deallocate( momenta )
  deallocate( force )
  deallocate( potential )
  deallocate( particleKineticEnergy )
  deallocate( particlePotential )
  deallocate( correlation )

end program


!!! Internal (initialization) subroutines below this line. !!!


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
  real(8) :: temp, mass, stand_dev, kb, dummy
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
     call get_gaussian_momentum(mass,stand_dev,momentum(counter,1),momentum(counter,2))
     call get_gaussian_momentum(mass,stand_dev,momentum(counter,3),dummy)
  end do
  !---------------------------------------------------------------
  !SRH: applying the standard deviation to the momenta - for some reason
  !the Gaussian RNG algorithm does not give any useful standard deviations!
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

  pxvar(:) = momentum(:,1)*momentum(:,1)
  pyvar(:) = momentum(:,2)*momentum(:,2)
  pzvar(:) = momentum(:,3)*momentum(:,3)
  conversion = (sum(pxvar)+sum(pyvar)+sum(pzvar))/N
  print *, "conversion: ", conversion, stand_dev

  do counter=1,N
     p_x = p_x + momentum(counter,1)
     p_y = p_y + momentum(counter,2)
     p_z = p_z + momentum(counter,3)
  end do


  print *, "Initial average momentum (before correction): ",p_x,p_y,p_z

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
  
	

subroutine get_gaussian_momentum(mass,stand_dev,p1,p2)

  real(8), intent(in) :: stand_dev, mass
  real(8), intent(out) :: p1, p2

  real(8) :: random_1,random_2, I, twopi, test, povers, dI
  integer :: counter, check

  check = 0
  twopi = 8d0*atan(1.0d0)
  do while (check == 0)
     call random_number(random_1)
     call random_number(random_2)
     if(random_1 > 1d-64) then
       p1 = dsqrt(-2d0*log(random_1))
       p2 = p1*sin(twopi*random_2)*stand_dev
       p1 = p1*cos(twopi*random_2)*stand_dev
           ! The standard deviation is most efficiently applied
           ! by simple multiplication of p/stand_dev by stand_dev.
       check = 1
     end if
  end do  

end subroutine get_gaussian_momentum
