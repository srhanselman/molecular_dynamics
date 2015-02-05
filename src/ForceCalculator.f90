!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ForceCalculator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
private

public force_calculator
public force_potential_calculator



contains


subroutine force_calculator(rm,bondingEnergy,x,f)

	integer ::		  i, j, N
	real*16, intent(in) :: 	  rm, bondingEnergy, x(:,:)
	real*16, intent(inout) :: f(:,:)
	real*16 ::		  dx(size(x,1),size(x,1),3), distancesq(size(x,1),size(x,1))

	N = size(x,1)

	do i=1,N
		do j=1,3
			dx(i,i,j) = 0
			distancesq(i,i) = 1
		end do
		do j=i+1,N
			dx(i,j,:) = x(i,:) - x(j,:)
			dx(j,i,:) = -dx(i,j,:)
			distancesq(i,j) = dot_product(dx(i,j,:),dx(i,j,:))
			distancesq(j,i) = distancesq(j,i)
		end do
	end do
	f = 12*bondingEnergy*sum( dx*(rm**6/spread(distancesq**4,1,N)-rm**12/spread(distancesq**7,1,N)) , 2 )

end subroutine


subroutine force_potential_calculator(rm,bondingEnergy,x,f,potential)

	real*16, intent(in) :: 	  rm, bondingEnergy, x(:,:)
	real*16, intent(inout) :: f(:,:)
	real*16, intent(out) ::   potential(:)
	real*16 ::		  dx(size(x,1),size(x,1),3), distance(size(x,1),size(x,1)), distancesq(size(x,1),size(x,1))
	integer ::		  i, j, N
	
	N = size(x,1)

	do i=1,N
		do j=1,3
			dx(i,i,j) = 0
			distancesq(i,i) = 1                              ! distancesq is the distance !squared!.
		end do
		do j=i+1,N
			dx(i,j,:) = x(i,:) - x(j,:)
			dx(j,i,:) = -dx(i,j,:)
			distancesq(i,j) = dot_product(dx(i,j,:),dx(i,j,:))
			distancesq(j,i) = distancesq(j,i)
		end do
	end do
	f = 12*bondingEnergy*sum( dx*(rm**6/spread(distancesq**4,1,N)-rm**12/spread(distancesq**7,1,N)) , 2 )
	potential = bondingEnergy*sum( -rm**6/(distancesq**3) + rm**12/(distancesq**6) , 2 )

end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ForceCalculator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
