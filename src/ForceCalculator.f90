!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ForceCalculator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
private

public force_calculator
public force_potential_calculator



contains


subroutine force_calculator(rm,bondingEnergy,maxForceDistance,L,x,f)
	! One generally would not need to directly call this subroutine.
	! rm is the potential minimum distance, bondingEnergy is the interaction depth,
	! L is the box length, x contains the particle positions while f contains the
	! force output.

	integer ::		  i, j, N
	real*8, intent(in) :: 	  rm, bondingEnergy, x(:,:), L, maxForceDistance
	real*8, intent(out) ::    f(:,:)
	real*8 ::		  dx(3), distancesq, fPair(3)

	N = size(x,1)
	f = 0
	do i=1,N
		do j=i+1,N
			dx = x(i,:) - x(j,:)
			dx = dx-nint(dx/L)*L
			distancesq = dot_product(dx,dx)
			if (distancesq < maxForceDistance*maxForceDistance) then
				distancesq = 1/distancesq
				fPair =  12*bondingEnergy*dx*(rm**6*distancesq**4-rm**12*distancesq**7)
				f(i,:) = f(i,:) + fPair
				f(j,:) = f(j,:) - fPair
			end if
		end do
	end do


end subroutine


subroutine force_potential_calculator(rm,bondingEnergy,maxForceDistance,L,x,f,potential)

	real*8, intent(in) :: 	  rm, bondingEnergy, x(:,:), L, maxForceDistance
	real*8, intent(inout) ::  f(:,:)
	real*8, intent(out) ::    potential(size(x,1))
	real*8 ::		  dx(3), distancesq
	real*8 ::                 fPair(3), potentialPair
	integer ::		  i, j, N
	
	N = size(x,1)

	do i=1,N
		do j=i+1,N
			dx = x(i,:) - x(j,:)
			dx = dx-nint(dx/L)*L
			distancesq = dot_product(dx,dx)
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
