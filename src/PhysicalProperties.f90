!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module PhysicalProperties
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
private

public calculate_temperature                         ! This subroutine calculates the temperature.
public calculate_pressure                            ! A stub as of 12/02/2015
public correlation_function_scalar		     ! This one compiles, as opposed to the external array one;
						! this subroutine will be called by ForceCalculator.

contains


subroutine calculate_temperature (particleKineticEnergy,meanMomentumSq,mass,temperatureCalc)
	! This subroutine uses the equipartition theorem (in which (<p^2>-<p>^2)/2m = 3kB*T/2)
	! to calculate the temperature (T = (<p^2>-<p>^2)/(3kB*m)) = Ekin(2/3kB)-<p>^2/3kB*m.
	real*8, intent(in) ::   particleKineticEnergy(:),meanMomentumSq, mass
	real*8, intent(out) ::  temperatureCalc
	real*8 ::               kB

	kB = 1.3806488d-2/(1.65d0)
	
	temperatureCalc = 2*sum(particleKineticEnergy)/size(particleKineticEnergy) - meanMomentumSq/mass
	temperatureCalc = temperatureCalc/(3*kB)
end subroutine

subroutine calculate_pressure (x,f,L,N,currentTemp,pressure)
	! This subroutine uses the virial theorem to calculate the pressure of the multi-particle model.
	! PV = N*kB*T(1 - f(x,F)/3N) = NkBT + <r.F>/3
	real*8, intent(in)  :: x(:,:), f(:,:), L, currentTemp
	real*8, intent(out) :: pressure
	real*8              :: volume, kB, positionForce(N)
	integer, intent(in) :: N
	integer             :: i

	volume = L*L*L
	kB = 1.3806488d-2/(1.65d0)

	do i=1,N
		positionForce(i) = dot_product(x(i,:),f(i,:))
	end do

        pressure = N*kB*currentTemp + sum(positionForce)/(3*N)
        pressure = pressure/volume
        
end subroutine


subroutine correlation_function_scalar (particleDistanceSq,dr,L,g)
	! This subroutine divides one single relative distance between by its spherical shell depth,
	! adding it to histogram g(#,r;dr).

	real*8, intent(in)  :: dr, L, particleDistanceSq
	real*8, intent(inout) :: g(:)
	real*8              :: pi
	integer             :: i

	pi = atan2(0d0,-1d0)

	! The histogram is based on spherical shells centred -around- g(1,:).
	! Hence, their corresponding volumes are 4pi/3((r+dr/2)**3-(r-dr/2)**3) = 4pi/3(3r**2*dr+(dr**3)/4)
	! or dV = (4*r**2+dr**2/3)*dr*pi = (4pi*i**2*dr**3 +dr**2/3)*dr = (4*i**2 + 1/3)*pi*dr**3

	i = nint(sqrt(particleDistanceSq)/dr)
	g(i) = g(i) + (4.*i**2 + 1./3)*pi*dr**3

end subroutine




!subroutine correlation_function_array (particleDistanceSq,dr,L,g)
	! This subroutine divides all relative distances between by their spherical shell depths,
	! plotting them to histogram g(#,r;dr).
	! Note that the particle distance matrix includes N(N-1)/2 valid values; all others (N(N+1)/2) are zeroes
	! and should be subtracted from the correlation histogram.

!	real*8, intent(in)  :: particleDistanceSq(:,:), dr, L
!	real*8, intent(out) :: g(ceiling(sqrt(3.)*L/(2*dr))+1,2)
!	real*8              :: pi
!	integer             :: i, N, particleDistances(size(particleDistanceSq,1)*size(particleDistanceSq,1))

!	pi = atan2(0d0,-1d0)
!	N = size(particleDistanceSq,1)	
!	particleDistances = nint(sqrt(reshape(particleDistanceSq,N*N))))

!	do i=1,size(g,1)
!		g(i,1) = dr*i
!		g(i,2) = 0
!	end do
!
!	do i=1,N
!		g(particleDistances(i),2) = g(particleDistances(i),2)+1
!	end do
!	
!	g(1,2) = g(1,2) - N*(N+1)/2
!
!	! The histogram is based on spherical shells centred -around- g(1,:).
!	! Hence, their corresponding volumes are 4pi/3((r+dr/2)**3-(r-dr/2)**3) = 4pi/3(3r**2*dr+(dr**3)/4)
!	! or dV = (4*r**2+dr**2/3)*dr*pi = (4pi*i**2*dr**3 +dr**2/3)*dr = (4*i**2 + 1/3)*pi*dr**3
!
!	do i=1,N
!		g(particleDistances(i),2) = g(particleDistances(i),2)/(4.*particleDistances(i)**2+1/3.)
!		g(particleDistances(i),2) = g(particleDistances(i),2)/(pi*dr*dr*dr)
!	end do
!
!end subroutine
	
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module PhysicalProperties
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

