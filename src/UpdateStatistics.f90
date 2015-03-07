!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module UpdateStatistics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
private

public printValueStatistics
!public printGradientStatistics


contains



subroutine printValueStatistics(runTemp,runTempCalc,runKineticEnergy,runPotentialEnergy,runPressure,runVolume)

! A void function which only prints the values, stdevs and covariances with temperature
! of thermodynamical properties.
	real*8, intent(in):: runTempCalc(:),runKineticEnergy(:),runPotentialEnergy(:), &
 runPressure(:),runVolume(:),runTemp(:)
	real*8 :: meanTempCalc,stdevTempCalc,meanKineticEnergy,stdevKineticEnergy, &
 covKineticEnergy,meanPotentialEnergy,stdevPotentialEnergy,meanPressure,stdevPressure, &
 covPressure,meanVolume,stdevVolume,covVolume,stdevTotalEnergy,meanTotalEnergy,covTotalEnergy, &
 covPotentialEnergy
	integer :: N

	N = size(runVolume)
	
	meanTempCalc = sum(runTempCalc)/N
	meanKineticEnergy = sum(runKineticEnergy)/N
	meanPotentialEnergy = sum(runPotentialEnergy)/N
	meanPressure = sum(runPressure)/N
	meanVolume = sum(runVolume)/N
	meanTotalEnergy = meanKineticEnergy + meanPotentialEnergy

	stdevTempCalc = dsqrt(dot_product(runTempCalc,runTempCalc)/N - meanTempCalc*meanTempCalc)
	stdevKineticEnergy = dsqrt(dot_product(runKineticEnergy,runKineticEnergy)/N - &
 meanKineticEnergy*meanKineticEnergy)
	stdevPotentialEnergy = dsqrt(dot_product(runPotentialEnergy,runPotentialEnergy)/N - &
 meanPotentialEnergy*meanPotentialEnergy)
	stdevPressure = dsqrt(dot_product(runPressure,runPressure)/N - meanPressure*meanPressure)
	stdevVolume = dsqrt(dot_product(runVolume,runVolume)/N - meanVolume*meanVolume)
	stdevTotalEnergy = dsqrt(dot_product(runKineticEnergy+runPotentialEnergy, &
 runKineticEnergy+runPotentialEnergy)/N - (meanKineticEnergy+meanPotentialEnergy)**2)

	covKineticEnergy = dsqrt(dot_product(runTempCalc,runKineticEnergy)/N - &
 meanTempCalc*meanKineticEnergy)
	covPotentialEnergy = dsqrt(dot_product(runTempCalc,runPotentialEnergy)/N - &
 meanTempCalc*meanPotentialEnergy)
	covTotalEnergy = dsqrt(dot_product(runTempCalc,runPotentialEnergy)/N - &
 meanTempCalc*(meanPotentialEnergy+meanKineticEnergy))
	covPressure = dsqrt(dot_product(runTempCalc,runPressure)/N - meanTempCalc*meanPressure)
	covVolume = dsqrt(dot_product(runTempCalc,runPressure)/N - meanTempCalc*meanVolume)

	
	write (3,*) runTemp,";",meanTempCalc,";",stdevTempCalc,";",meanKineticEnergy,";", &
 stdevKineticEnergy,";",covKineticEnergy,";",meanPotentialEnergy,";",stdevPotentialEnergy,";", &
 covKineticEnergy,";",meanTotalEnergy,";",stdevTotalEnergy,";",covTotalEnergy,";",meanPressure, &
 ";",stdevPressure,";",covPressure,";",meanVolume,";",stdevVolume,";",covVolume

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module UpdateStatistics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
