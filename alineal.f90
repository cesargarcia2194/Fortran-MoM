MODULE alineal
use constantes
IMPLICIT NONE

contains
  	FUNCTION indexrot4(n) RESULT(res)
    	INTEGER, INTENT(IN) :: n
    	INTEGER :: res
	
    	res = n
	
    	IF(res>4) THEN
    	   res = res - ((res-1)/4)*4
    	END IF
END FUNCTION indexrot4
	FUNCTION indexrot3(n) RESULT(res)
    	INTEGER, INTENT(IN) :: n
    	INTEGER :: res
	
    	res = n
	
    	IF(res>3) THEN
    	   res = res - ((res-1)/3)*3
    	END IF
END FUNCTION indexrot3
	!Si los valores triplet 1 estan en triplet 2 devuel TRUE (no importa el orden)
	FUNCTION comp_trio(triplet1, triplet2) RESULT(res)
	    INTEGER, DIMENSION(3), INTENT(IN) :: triplet1, triplet2
	    LOGICAL :: res
	    ! Revisar esta declaracion y como es la matriz de verdad
	    INTEGER, PARAMETER, DIMENSION(3,6) :: indices = reshape((/1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1/),(/3,6/))
	    !1 1 2 2 3 3
	    !2 3 1 3 1 2
	    !3 2 3 1 2 1
	    INTEGER :: i
	
	    res = .FALSE.
	
	    DO i=1,6
	       IF(triplet1(1)==triplet2(indices(1,i)) .AND.&
	            triplet1(2)==triplet2(indices(2,i)) .AND.&
	            triplet1(3)==triplet2(indices(3,i))) THEN
	          res = .TRUE.
	          RETURN
	       END IF
	    END DO
END FUNCTION comp_trio

	FUNCTION comp_par(par1, par2)RESULT(res)

		integer,DIMENSION(2), INTENT(IN)::par1, par2
		logical :: res
		res= .FALSE.
		if( (par1(1)==par2(1) .AND. par1(2) == par2(2)) .OR.&
			(par1(1)==par2(2) .AND. par1(2) == par2(1)))then
			res= .TRUE.
		end if
END FUNCTION comp_par

	function prod_cruz(v1,v2) RESULT(v1xv2)
		real(KIND=dp),DIMENSION(3),INTENT(IN)::v1
		real(KIND=dp),DIMENSION(3),INTENT(IN)::v2
		real(KIND=dp),DIMENSION(3)::v1xv2

		v1xv2(1)=v1(2)*v2(3)-v2(2)*v1(3)
		v1xv2(2)=v2(1)*v1(3)-v1(1)*v2(3)
		v1xv2(3)=v1(1)*v2(2)-v2(1)-v1(2)
	end function prod_cruz
	function prod_escalar(v1,v2) RESULT(v1v2)
		real(KIND=dp),DIMENSION(3),INTENT(IN)::v1
		real(KIND=dp),DIMENSION(3),INTENT(IN)::v2
		real::v1v2

		v1v2 = SUM(v1*v2)
	end function prod_escalar

	function normav(v)RESULT(res)
		real(KIND=dp),DIMENSION(3),INTENT(IN)::v
		real(KIND=dp)::res
		res=SQRT(v(1)**2 +v(2)**2 +v(3)**2)
	end function normav

	function hace_vertor3(v1,v2)RESULT(v)
		real(KIND=dp),DIMENSION(3),INTENT(IN)::v1
		real(KIND=dp),DIMENSION(3),INTENT(IN)::v2
		real(KIND=dp),DIMENSION(3)::v
		v(1)=v2(1)-v1(1)
		v(2)=v2(2)-v1(2)
		v(3)=v2(3)-v1(3)
	end function hace_vertor3


END MODULE alineal