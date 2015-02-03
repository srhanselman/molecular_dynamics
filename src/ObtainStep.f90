!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ObtainStep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
private

public obtain_step

contains

subroutine obtain_step(x,p,f,L,maxBeta,tolerance,mass,timeStep,outOfVelocityBounds,meanMomentumSq,meanMomentumDiverged)
	! This function uses Euler's method.
	! The subroutine changes x, p to their next step values.

	integer :: i, j                                                                  ! a plain old iterator
	real*16, intent(in) ::    maxBeta, mass, timeStep, tolerance, L               ! maxBeta is dimensionless, mass in kg, timeStep in s, tolerance is dimensionless (# stddevs)
	real*16, intent(inout) :: x(:,:), p(:,:), f(:,:)                              ! position, momentum and force vectors
	real*16, intent(out) ::   meanMomentumSq
	logical, intent(out) ::   outOfVelocityBounds, meanMomentumDiverged           ! these flags signal clear divergences in individual and mean momentum
	real*16 ::                lightSpeed, temperature, stddevsq, kB, maxMomentum  ! lightSpeed is in m/s, temperature in K, stddevsq in kg^2m^2s^-2, kB in J/K, maxMomentum in kgms^-1
	kB = 1.3806488d-23
	temperature = 273
	lightSpeed = 299792458


	stddevsq = 2*mass*kB*temperature
	maxMomentum = maxBeta*lightSpeed*mass


	p = p + f * timeStep

	meanMomentumSq = dot_product(p(:,1),p(:,1))+dot_product(p(:,2),p(:,2))+dot_product(p(:,3),p(:,3))

	if (meanMomentumSq > tolerance*stddevsq) then
		meanMomentumDiverged = .true.
	else
		do i=1,size(p,1)
			if (dot_product(p(i,:),p(i,:)) > maxMomentum*maxMomentum) then
				outOfVelocityBounds = .true.
			else
				x = x + p*timeStep/mass + f*timeStep*timeStep/(2*mass)
                                do j=1,3
					if (x(i,j)<0 .OR. x(i,j)>L) then
						x(i,j) = mod(x(i,j),L)
					end if                
				end do                        
				meanMomentumDiverged = .false.
				outOfVelocityBounds = .false.
			end if
		end do
	end if


end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ObtainStep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

