!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ObtainStep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use ForceCalculator

! Changelog 2/4/2015 11:22:
!      meanMomentumSq is now defined as being a sum over all momenta.
! Changelog 2/5/2015:
!      included velocity Verlet algorithm


implicit none
private

public verlet_eqs_of_motion                          ! Call this subroutine when energy is presumed to be stable.
public verlet_eqs_of_motion_checks                   ! A lot slower, but automatically sets warning flags for clear transgressions.
public verlet_eqs_of_motion_correlation			     ! Integrates correlation function calculation.


contains

subroutine verlet_eqs_of_motion(x,p,f,L,maxForceDistance,rm,bondingEnergy,mass,timeStep,meanMomentumSq, &
 particleKineticEnergy,particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep, &
 particleDistanceSq)
	! This function uses a simple velocity Verlet algorithm.
	! Its arguments are as listed above; note that x, p, f are position, momentum
	!    and force respectively, and L is the box length. outOfVelocityBounds 
	!    indicates whether the program should warn/exit due to exceeding the initially
	!    defined velocities, whereas meanMomentumDiverged does the same for a
	!    out-of-bounds shift of mean momentum. rm indicates the individual interaction
	!    potential minimum while bondingEnergy indicates the potential depth.

	integer :: i                                                                  ! a plain old iterator
	real*8, intent(in) ::    mass, timeStep, rm, bondingEnergy, L, maxForceDistance ! maxBeta is dimensionless, mass in kg, L in m, timeStep in s, momentumTolerance is dimensionless (# stddevs)
	real*8, intent(inout) :: x(:,:), p(:,:), f(:,:)                              ! position, momentum and force vectors
	real*8, intent(out) ::   kineticEnergyAfterStep, potentialEnergyAfterStep, &
 particleKineticEnergy(size(x,1)), particlePotential(size(x,1)), &
 particleDistanceSq(size(x,1),size(x,1))
		! these flags signal clear divergences in individual and mean momentum
	real*8 ::                meanMomentumSq

!$omp parallel
        p = p + f*timeStep/2
	x = modulo(x+p*timeStep/mass,L)
	call force_potential_calculator(rm,bondingEnergy,maxForceDistance,L,x,f,particlePotential, &
	 particleDistanceSq)
	p = p + f * timeStep/(2*mass)
!$omp end parallel	

	meanMomentumSq = dot_product(sum(p(:,:),1),sum(p(:,:),1))                     ! sums over all particles BEFORE obtaining the modulus squared - see Changelog 2/4/2015
	particleKineticEnergy = dot_product(p(i,:),p(i,:))/(2*mass)
	kineticEnergyAfterStep = sum(particleKineticEnergy)
	potentialEnergyAfterStep = sum(particlePotential)

end subroutine





subroutine verlet_eqs_of_motion_checks(x,p,f,L,maxBeta,momentumTolerance,mass,timeStep,outOfVelocityBounds, &
 kineticEnergyAfterStep,particleMomentumSq,initialKineticEnergy,particlePotential,potentialEnergyAfterStep, &
 energyTolerance,meanMomentumDiverged,rm,bondingEnergy,energyNotPreserved,maxForceDistance, particleDistanceSq)

	! To automatically check for issues, use the flags displayed under logical, intent(out).

	integer :: i                                                               ! a plain old iterator
	real*8, intent(in) ::    maxBeta, mass, timeStep, momentumTolerance, rm, bondingEnergy, energyTolerance, &
 initialKineticEnergy, L, maxForceDistance
           ! maxBeta is dimensionless, mass in kg, timeStep in s, momentumTolerance/energyTolerance are dimensionless (# stddevs)
	real*8, intent(inout) :: x(:,:), p(:,:), f(:,:)                              ! position, momentum and force vectors
	real*8, intent(out) ::   kineticEnergyAfterStep, particlePotential(size(x,1)), potentialEnergyAfterStep, &
 particleMomentumSq(size(x,1)),  particleDistanceSq(size(x,1),size(x,1))
	logical, intent(out) ::   outOfVelocityBounds(size(x,1)), meanMomentumDiverged, energyNotPreserved ! these flags signal clear divergences in individual and mean momentum
				
	real*8 ::  lightSpeed, temperature, stddevsq, kB, maxMomentum
        real*8 ::  particleKineticEnergy(size(x,1))
        real*8 ::  meanMomentumSq  ! lightSpeed is in m/s, 
!temperature in K, stddevsq in kg^2m^2s^-2, kB in J/K, maxMomentum in kgms^-1
	kB = 1/1.65d-21
	temperature = 273
	lightSpeed = 299792458

	stddevsq = 2*mass*kB*temperature
	maxMomentum = maxBeta*lightSpeed*mass
 

        p = p + 0.5*f*timeStep

	do i=1,size(p,1)
		particleMomentumSq(i) = dot_product(p(i,:),p(i,:))
		if (particleMomentumSq(i) > maxMomentum*maxMomentum) then
			outOfVelocityBounds(i) = .true.
		else
			x(i,:) = modulo(x(i,:)+p(i,:)*timeStep/mass,L)
			outOfVelocityBounds(i) = .false.
		end if
	end do


	call force_potential_calculator(rm,bondingEnergy,maxForceDistance,L,x,f,particlePotential,particleDistanceSq)
	p = p + f * timeStep/2
	

	meanMomentumSq = dot_product(sum(p(:,:),1),sum(p(:,:),1))                     ! sums over all particles BEFORE obtaining the modulus squared - see Changelog 2/4/2015
	meanMomentumDiverged = (meanMomentumSq > momentumTolerance*stddevsq)

	particleKineticEnergy = particleMomentumSq/(2*mass)
	kineticEnergyAfterStep = sum(particleKineticEnergy)
	potentialEnergyAfterStep = sum(particlePotential)/2

	energyNotPreserved = (abs(kineticEnergyAfterStep-potentialEnergyAfterStep) > kB*temperature)


end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine verlet_eqs_of_motion_correlation(x,p,f,L,maxForceDistance,rm,bondingEnergy,mass,timeStep,meanMomentumSq, &
 particleKineticEnergy,particlePotential,kineticEnergyAfterStep,potentialEnergyAfterStep,dR,correlation)
	! This function uses a simple velocity Verlet algorithm.
	! Its arguments are as listed above; note that x, p, f are position, momentum
	!    and force respectively, and L is the box length. outOfVelocityBounds 
	!    indicates whether the program should warn/exit due to exceeding the initially
	!    defined velocities, whereas meanMomentumDiverged does the same for a
	!    out-of-bounds shift of mean momentum. rm indicates the individual interaction
	!    potential minimum while bondingEnergy indicates the potential depth.

	integer :: i, N                                                                ! a plain old iterator
	real*8, intent(in) ::    mass, timeStep, rm, bondingEnergy, L, maxForceDistance, dR
             ! maxBeta is dimensionless, mass in kg, L in m, timeStep in s, momentumTolerance is dimensionless (# stddevs)
	real*8, intent(inout) :: x(:,:), p(:,:), f(:,:), correlation(:)
           ! position, momentum and force vectors
	real*8, intent(out) ::   kineticEnergyAfterStep, potentialEnergyAfterStep, particlePotential(size(x,1)), &
 particleKineticEnergy(size(x,1)) 
		! these flags signal clear divergences in individual and mean momentum
	real*8 ::                meanMomentumSq, psq(size(x,1))

	psq = 0
	N = size(x,1)

        p = p + f*timeStep/(2*mass)
	x = modulo(x+p*timeStep/mass,L)
	call force_potential_calculator_correlation(rm,bondingEnergy,maxForceDistance,L,x,f,particlePotential,dR,correlation)
	p = p + f*timeStep/(2*mass)
	

	meanMomentumSq = dot_product(sum(p(:,:),1),sum(p(:,:),1))                     ! sums over all particles BEFORE obtaining the modulus squared - see Changelog 2/4/2015
	do i=1,size(p,1)
		psq(i) = dot_product(p(i,:),p(i,:))
	end do
	particleKineticEnergy = psq/(2*mass)
	meanMomentumSq = meanMomentumSq/N
	kineticEnergyAfterStep = sum(particleKineticEnergy)
	potentialEnergyAfterStep = sum(particlePotential)/2

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ObtainStep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

