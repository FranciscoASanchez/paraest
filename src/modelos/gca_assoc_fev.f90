        SUBROUTINE ASSOCGCEOS(NCOD,NC,DC,RN,TEMPCRIT,V,TEMP,NG,NU,QAREA&
     &,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC&
     &,ENASS,RKASS,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,ierrassoc)
!
! HAY UN COMMON
!
! CALCULA LO MISMO QUE ASSOCGCEOS3, PERO HACE UN CHEQUEO, COMPARANDO
!   DERIVADAS CALCULADAS EN FORMA NUMERICA CON DERIVADAS CALCULADAS
!	EN FORMA ANALITICA en el caso en que NOPCION = 1
!
!    SI NOPCION = 0 ESTE PROGRAMA HACE EXACTAMENTE LO MISMO QUE ASSOCGCEOS3
!
!    IERRDER = 0  ===> NO PROBLEM
!    IERRDER = 1  ===> LA DERIVACION NUMERICA DIFIERE MUCHO DE LA ANALITICA
!				POSIBLE ERROR DE PROGRAMACION
!
!
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50, NGM = 20
!
        DIMENSION DC(NC), RN(NC)
        DIMENSION  TEMPCRIT(NC)
!
! DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL CARACTER GRUPAL
! DEL TERMINO DE DISPERSION
!
        DIMENSION NU(NGM,NC), QAREA(NGM), TSTR(NGM), GSTR(NGM)
        DIMENSION GDOT(NGM), GDDOT(NGM), RKIJ(NGM,NGM), RKTIJ(NGM,NGM)
        DIMENSION ALFAIJ(NGM,NGM)
!
! DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL TERMINO DE ASOCIACION
!
        DIMENSION ONUASOC(12,NC)
        DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
!
        DIMENSION FFI(NC)  !DIMENSIONAMIENTO COEFICIENTES DE FUGACIDAD
!
      DIMENSION RNCTRL(NCM), FFIANAL(NCM)
!
      COMMON/DERIVNUMER/NOPCION
!
      IERRDER = 0
!
      IF(NOPCION.EQ.1) NCODRES = NCOD
      IF(NOPCION.EQ.1) NCOD = 5
!
!   EN ESTE PRIMER LLAMADO SE CALCULA TODO: Egia de Helmholtz, factor de
!     compresibilidad y coeficiente de fugacidad SOLO si NOPCION = 1
!      Si NOPCION = 0 se respeta el valor de NCOD original
!
        CALL ASSOCGCEOS3(NCOD,NC,DC,RN,TEMPCRIT,V,TEMP,NG,NU,QAREA&
     &,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC&
     &,ENASS,RKASS,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,ierrassoc)
!
! FACTORES PARA CALCULO DE DERIVADAS EN FORMA NUMERICA
!
      FVOL = 1.D-8
      FMOL = 1.D-8
!
      IF(NOPCION.EQ.1) THEN
!		------------------------------------------
!
!		ALMACENAMIENTO DE INFORMACION ANALITICA
!
      	VCTRL = V  !variable de entrada
!
!		ZSISANAL = FACTOR DE COPRES. DEL SISTEMA CALCULADO EN
!		FORMA ANALITICA
!
      	ZSISANAL = ZSIS    !variable de salida
      	ARESCTRL = ARESTOT !variable de salida
!
      	SUMAR = 0.
!
      	DO I=1,NC
      		RNCTRL(I)=RN(I)  !variable de entrada
!
!			FFIANAL(I)= COEF. DE FUGAC. CALCULADOS EN FORMA
!					ANALITICA
!
      		FFIANAL(I)=FFI(I) !variable de salida
!
      		SUMAR = SUMAR + RN(I)
      	END DO
!
      	RMOLTOT = SUMAR !Nro. total de moles
!
!		Perturbacion del Volumen
!
      	DELTAV = FVOL * VCTRL
      	V = VCTRL + DELTAV
!
      	NCOD = 4 !calculo solo de energia residual
        CALL ASSOCGCEOS3(NCOD,NC,DC,RN,TEMPCRIT,V,TEMP,NG,NU,QAREA&
     &,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC&
     &,ENASS,RKASS,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,&
     &ierrassoc)
!
!     		DERIFV= DERIVADA DE ARESTOT CON RESPECTO AL VOLUMEN V
!
      	DERIFV = (ARESTOT - ARESCTRL)/ DELTAV
!
!		ZSISNUMER = Factor de compresib. del sistema calculado en forma
!				numerica
!
      	ZSISNUMER = (-1.) * (VCTRL/RMOLTOT) * DERIFV  + 1.
!
!
      	if(dabs(ZSISNUMER-ZSISANAL).eq.0) THEN
!
      	  OPEN(UNIT=1,STATUS='OLD',FILE='ASSOCPARA.RES',&
     & FORM='FORMATTED',&
     & POSITION='APPEND')
!
      	  write(1,*)&
     &' WARNING: El paso para la dif. numerica resp. al vol &
     &es demasiado peq.'
!
      	  CLOSE(UNIT=1)
      	 ELSE
      	  CONTINUE
      	END IF
!
      	IF(DABS((ZSISNUMER-ZSISANAL)/ZSISANAL).GT.1.D-4) THEN
!
      	  OPEN(UNIT=1,STATUS='OLD',FILE='ASSOCPARA.RES',&
     & FORM='FORMATTED',&
     & POSITION='APPEND')
!
      	  WRITE(1,*) ' LA DERIVACION NUMERICA CON RESPECTO AL VOLUMEN&
     & NO COINCIDE CON LA DERIVACION ANALITICA'
      	  WRITE(1,*) ' Zanal =',ZSISANAL, ' Znum = ',ZSISNUMER
!
      	  CLOSE(UNIT=1)
!
!		EN ESTE CASO LAS POSIBILIDADES SON:
!
!		 - POSIBLE ERROR DE PROGRAMACION
!     		 - VALOR DEMASIADO ALTO DEL PASO. Esto es FVOL  muy alto
!		 - Valor demasiado bajo del paso. Esto es FVOL  muy bajo
!							para la precision de la
!							maquina
!
      	   IERRDER=1
      	 ELSE
      	  CONTINUE
      	END IF
!
!		RETORNO el volumen a su valor inicial
!
      	V = VCTRL
!
!		CALCULO DE COEFICIENTES DE FUGACIDAD EN FORMA NUMERICA
!
      	DO I=1,NC
!			---------------------------------
!			Perturbacion del nro. de moles
!			Afecto solo al numero de moles del compuesto I
!
      		DELTAM = FMOL * RN(I)
      		RN(I) = RNCTRL(I) + DELTAM
!
      	        NCOD = 4 !calculo solo de energia residual
        CALL ASSOCGCEOS3(NCOD,NC,DC,RN,TEMPCRIT,V,TEMP,NG,NU,QAREA&
     &,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC&
     &,ENASS,RKASS,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP&
     &,ierrassoc)
!
      		DERIFN =  (ARESTOT - ARESCTRL)/ DELTAM
!
      		IF(ZSISANAL.GT.0) THEN
!                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      		FFINUMER = DEXP( DERIFN - DLOG(ZSISANAL) )
!
      		IF(DABS(FFINUMER-FFIANAL(I)).EQ.0) THEN
!
      		 OPEN(UNIT=1,STATUS='OLD',FILE='ASSOCPARA.RES',&
     & FORM='FORMATTED',&
     & POSITION='APPEND')
!
      		 WRITE(1,*)&
     &' WARNING: El paso para la diferenciacion numerica con respecto &
     &  al numero de moles es demasiado chico'
!
      		 CLOSE(UNIT=1)
!
      		 ELSE
      		  CONTINUE
      		END IF
!
      		IF(DABS((FFINUMER-FFIANAL(I))/FFIANAL(I)).GT.1.D-4) THEN
!
      		 OPEN(UNIT=1,STATUS='OLD',FILE='ASSOCPARA.RES',&
     & FORM='FORMATTED',&
     & POSITION='APPEND')
!
      	  	 WRITE(1,*) ' LA DERIVACION NUMERICA CON RESPECTO AL &
     & NUMERO DE MOLES DEL COMP ',I,' NO COINCIDE CON LA&
     & DERIVACION ANALITICA'
      	  	 WRITE(1,*) ' COEFanal = ',FFIANAL(I),&
     &' COEFnum = ', FFINUMER
!
      		 CLOSE(UNIT=1)
!
      		  IERRDER = 1
      		 ELSE
      		  CONTINUE
      		END IF
!
!                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      		 ELSE
      		END IF
!
!
!			Recupero el valor del numero de moles
!
      		RN(I) = RNCTRL(I)
!
!			---------------------------------
      	END DO
!
!		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!		Recuperacion de los valores adecuados en las variables
!		 de salida que fueron perturbadas

!
      	ARESTOT = ARESCTRL
      	ZSIS = ZSISANAL
!
!
      	DO I=1,NC
!
      		FFI(I) = FFIANAL(I)
!
      	END DO
!
!		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      	NCOD = NCODRES
!		------------------------------------------
!
       ELSE
        CONTINUE
      END IF
!
!
      END
!
!///////////////////////////////////////////////////////////////////////
