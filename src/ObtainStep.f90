!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ObtainStep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use ForceCalculator

! Changelog 2/4/2015 11:22: 
!      meanMomentumSq is now defined as being a sum over all momenta.


implicit none
private

public obtain_step

contains

subroutine obtain_step(x,p,f,maxBeta,momentumTolerance,mass,timeStep,outOfVelocityBounds,kineticEnergyAfterStep, &
meanMomentumDiverged,rm,bondingEnergy)
	! This function uses a simple velocity Verlet algorithm.
	! Its arguments are as listed above; note that x, p, f are position, momentum
	!    and force respectively; outOfVelocityBounds indicates whether the program
	!    should warn/exit due to exceeding the initially defined velocities,
	!    whereas meanMomentumDiverged does the same for a out-of-bounds shift of
	!    mean momentum. rm indicates the individual interaction potential minimum
	!    whereas 

	integer :: i                                                                  ! a plain old iterator
	real*16, intent(in) ::    maxBeta, mass, timeStep, momentumTolerance, rm, bondingEnergy ! maxBeta is dimensionless, mass in kg, timeStep in s, momentumTolerance is dimensionless (# stddevs)
	real*16, intent(inout) :: x(:,:), p(:,:), f(:,:)                              ! position, momentum and force vectors
	real*16, intent(out) ::   kineticEnergyAfterStep
	logical, intent(out) ::   outOfVelocityBounds(size(x,1)), meanMomentumDiverged 
		! these flags signal clear divergences in individual and mean momentum
	real*16 ::                lightSpeed, temperature, stddevsq, kB, maxMomentum, meanMomentumSq, particleMomentumSq(size(x,1)) 
	real*16 ::                fNext(size(x,1),3)  ! lightSpeed is in m/s, temperature in K, stddevsq in kg^2m^2s^-2, kB in J/K, maxMomentum in kgms^-1
	kB = 1.3806488d-23
	temperature = 273
	lightSpeed = 299792458


	stddevsq = 2*mass*kB*temperature
	maxMomentum = maxBeta*lightSpeed*mass

	p = p + f * timeStep

	if (meanMomentumSq > momentumTolerance*stddevsq) then
		meanMomentumDiverged = .true.
	else
		do i=1,size(p,1)
			particleMomentumSq(i) = dot_product(p(i,:),p(i,:))
			if (particleMomentumSq(i) > maxMomentum*maxMomentum) then
				outOfVelocityBounds(i) = .true.
			else
				x(i,:) = x(i,:) + p(i,:)*timeStep/mass + f(i,:)*timeStep*timeStep/(2*mass)
				outOfVelocityBounds(i) = .false.
			end if
		end do
	end if
	
	fNext = f
	
	call force_calculator(rm,bondingEnergy,x,fNext)

	p = p + (f+fNext) * timeStep/2

	f = fNext
	
	meanMomentumSq = dot_product(sum(p(:,:),1),sum(p(:,:),1))                     ! sums over all particles BEFORE obtaining the modulus squared - see Changelog 2/4/2015
	meanMomentumDiverged = (meanMomentumSq > momentumTolerance*stddevsq)

end subroutine





subroutine obtain_step_energy_check(x,p,f,maxBeta,momentumTolerance,mass,timeStep,outOfVelocityBounds,kineticEnergyAfterStep, &
 initialKineticEnergy, potentialEnergyAfterStep,energyTolerance,meanMomentumDiverged,rm,bondingEnergy,energyNotPreserved)


	integer :: i                                                               ! a plain old iterator
	real*16, intent(in) ::    maxBeta, mass, timeStep, momentumTolerance, rm, bondingEnergy, energyTolerance, initialKineticEnergy ! maxBeta is dimensionless, mass in kg, timeStep in s, momentumTolerance is dimensionless (# stddevs)
	real*16, intent(inout) :: x(:,:), p(:,:), f(:,:)                              ! position, momentum and force vectors
	real*16, intent(out) ::  kineticEnergyAfterStep, potentialEnergyAfterStep
	logical, intent(out) ::  outOfVelocityBounds(size(x,1)), meanMomentumDiverged, energyNotPreserved ! these flags signal clear divergences in individual and mean momentum
	
	real*16 ::  lightSpeed, temperature, stddevsq, kB, maxMomentum
        real*16 ::  fNext(size(x,1),3), potential(size(x,1)), particleMomentumSq(size(x,1))
        real*16 ::  meanMomentumSq  ! lightSpeed is in m/s, 
!temperature in K, stddevsq in kg^2m^2s^-2, kB in J/K, maxMomentum in kgms^-1
	kB = 1.3806488d-23
	temperature = 273
	lightSpeed = 299792458


	stddevsq = 2*mass*kB*temperature
	maxMomentum = maxBeta*lightSpeed*mass


	do i=1,size(p,1)
		particleMomentumSq(i) = dot_product(p(i,:),p(i,:))
		if (particleMomentumSq(i) > maxMomentum*maxMomentum) then
			outOfVelocityBounds(i) = .true.
		else
			x(i,:) = x(i,:) + p(i,:)*timeStep/mass + f(i,:)*timeStep*timeStep/(2*mass)
			outOfVelocityBounds(i) = .false.
		end if
	end do

	
	fNext = f	
	call force_potential_calculator(rm,bondingEnergy,x,fNext,potential)
	p = p + (f+fNext) * timeStep/2
	f = fNext

	meanMomentumSq = dot_product(sum(p(:,:),1),sum(p(:,:),1))                     ! sums over all particles BEFORE obtaining the modulus squared - see Changelog 2/4/2015
	meanMomentumDiverged = (meanMomentumSq > momentumTolerance*stddevsq)

	kineticEnergyAfterStep = sum(particleMomentumSq)/(2*mass)
	potentialEnergyAfterStep = sum(potential)

	energyNotPreserved = (abs(kineticEnergyAfterStep-potentialEnergyAfterStep) > kB*temperature)


end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ObtainStep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

