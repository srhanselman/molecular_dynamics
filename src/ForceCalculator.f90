!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ForceCalculator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! COMMENTS 02/12/2015
!	In calculating the correlation function, it may be convenient to store the r_{ij}^2
!	and simply reuse this array in the correlation function. However, this requires a major
!	rewrite of the force_(potential_)calculator subroutine, involving a second particleDistanceSq
!	array. I'll do just that, but in case we find a more efficient way of coding, we might
!	revert to the previous commit.

use PhysicalProperties

implicit none
private

public force_calculator
public force_potential_calculator
public force_calculator_correlation
public force_potential_calculator_correlation


contains


subroutine force_calculator(rm,bondingEnergy,maxForceDistance,L,x,f,particleDistanceSq)
	! One generally would not need to directly call this subroutine.
	! rm is the potential minimum distance, bondingEnergy is the interaction depth,
	! L is the box length, x contains the particle positions while f contains the
	! force output.

	integer ::		  i, j, N
	real*8, intent(in) :: 	  rm, bondingEnergy, x(:,:), L, maxForceDistance
	real*8, intent(out) ::    f(:,:), particleDistanceSq(size(x,1),size(x,1))
	real*8 ::		  dx(3), distancesq, fPair(3)

	N = size(x,1)
	f = 0
	particleDistanceSq = 0

	do i=1,N
		do j=i+1,N
			dx = x(i,:) - x(j,:)
			dx = dx-nint(dx/L)*L
			distancesq = dot_product(dx,dx)
			particleDistanceSq(i,j) = distancesq
			if (distancesq < maxForceDistance*maxForceDistance) then
				distancesq = 1/distancesq
				fPair =  12*bondingEnergy*dx*(rm**6*distancesq**4-rm**12*distancesq**7)
				f(i,:) = f(i,:) + fPair
				f(j,:) = f(j,:) - fPair
			end if
		end do
	end do


end subroutine


subroutine force_potential_calculator(rm,bondingEnergy,maxForceDistance,L,x,f,potential,particleDistanceSq)

	real*8, intent(in) :: 	  rm, bondingEnergy, x(:,:), L, maxForceDistance
	real*8, intent(inout) ::  f(:,:)
	real*8, intent(out) ::    potential(size(x,1)), particleDistanceSq(size(x,1),size(x,1))
	real*8 ::		  dx(3), distancesq
	real*8 ::                 fPair(3), potentialPair
	integer ::		  i, j, N
	
	N = size(x,1)
	particleDistanceSq = 0

	do i=1,N
		do j=i+1,N
			dx = x(i,:) - x(j,:)
			dx = dx-nint(dx/L)*L
			distancesq = dot_product(dx,dx)
			particleDistanceSq(i,j) = distancesq
			if (distancesq < maxForceDistance*maxForceDistance) then
				distancesq = 1/distancesq
				fPair =  12*bondingEnergy*dx*(rm**6*distancesq**4-rm**12*distancesq**7)
				potentialPair = bondingEnergy*(2*rm**6*distancesq**3-rm**12*distancesq**6)
				f(i,:) = f(i,:) + fPair
				f(j,:) = f(j,:) - fPair
				potential(i) = potential(i) + potentialPair
				potential(j) = potential(j) + potentialPair
			end if
		end do
	end do


end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS AFTER THIS LINE INCLUDE HISTOGRAM CREATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine force_calculator_correlation(rm,bondingEnergy,maxForceDistance,L,x,f,dR,correlation)
	! One generally would not need to directly call this subroutine.
	! rm is the potential minimum distance, bondingEnergy is the interaction depth,
	! L is the box length, x contains the particle positions while f contains the
	! force output.

	integer ::		  i, j, N
	real*8, intent(in) :: 	  rm, bondingEnergy, x(:,:), L, maxForceDistance, dR
	real*8, intent(inout) ::  f(:,:)
	real*8, intent(out) ::    correlation(ceiling(sqrt(3.)*L/(2*dR))+1)
	real*8 ::		  dx(3), distancesq, fPair(3)

	N = size(x,1)
	f = 0
	correlation = 0

	do i=1,N
		do j=i+1,N
			dx = x(i,:) - x(j,:)
			dx = dx-nint(dx/L)*L
			distancesq = dot_product(dx,dx)
			call correlation_function_scalar(distancesq,dR,L,correlation)
			if (distancesq < maxForceDistance*maxForceDistance) then
				distancesq = 1/distancesq
				fPair =  12*bondingEnergy*dx*(rm**6*distancesq**4-rm**12*distancesq**7)
				f(i,:) = f(i,:) + fPair
				f(j,:) = f(j,:) - fPair
			end if
		end do
	end do


end subroutine


subroutine force_potential_calculator_correlation(rm,bondingEnergy,maxForceDistance,L,x,f,potential,dR,correlation)

	real*8, intent(in) :: 	  rm, bondingEnergy, x(:,:), L, maxForceDistance, dR
	real*8, intent(inout) ::  f(:,:)
	real*8, intent(out) ::    potential(size(x,1)), correlation(ceiling(sqrt(3.)*L/(2*dr))+1)
	real*8 ::		  dx(3), distancesq
	real*8 ::                 fPair(3), potentialPair
	integer ::		  i, j, N
	
	N = size(x,1)
	correlation = 0

	do i=1,N
		do j=i+1,N
			dx = x(i,:) - x(j,:)
			dx = dx-nint(dx/L)*L
			distancesq = dot_product(dx,dx)
			call correlation_function_scalar(distancesq,dR,L,correlation)
			if (distancesq < maxForceDistance*maxForceDistance) then
				distancesq = 1/distancesq
				fPair =  12*bondingEnergy*dx*(rm**6*distancesq**4-rm**12*distancesq**7)
				potentialPair = bondingEnergy*(2*rm**6*distancesq**3-rm**12*distancesq**6)
				f(i,:) = f(i,:) + fPair
				f(j,:) = f(j,:) - fPair
				potential(i) = potential(i) + potentialPair
				potential(j) = potential(j) + potentialPair
			end if
		end do
	end do


end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ForceCalculator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
