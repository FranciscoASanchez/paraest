        SUBROUTINE ASSOCGCEOS(NCOD,NC,DC,RN,TEMPCRIT,V,TEMP,NG,NU,QAREA
     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
     @,ENASS,RKASS,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,ierrassoc)
C
C HAY UN COMMON
C
C CALCULA LO MISMO QUE ASSOCGCEOS3, PERO HACE UN CHEQUEO, COMPARANDO
C   DERIVADAS CALCULADAS EN FORMA NUMERICA CON DERIVADAS CALCULADAS
C	EN FORMA ANALITICA en el caso en que NOPCION = 1
C
C    SI NOPCION = 0 ESTE PROGRAMA HACE EXACTAMENTE LO MISMO QUE ASSOCGCEOS3
C
C    IERRDER = 0  ===> NO PROBLEM
C    IERRDER = 1  ===> LA DERIVACION NUMERICA DIFIERE MUCHO DE LA ANALITICA
C				POSIBLE ERROR DE PROGRAMACION
C
C
        IMPLICIT REAL*8 (A-H,O-Z)
        integer, parameter :: NCM = 50, NGM = 20
C
        DIMENSION DC(NC), RN(NC)
        DIMENSION  TEMPCRIT(NC)
C
C DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL CARACTER GRUPAL
C DEL TERMINO DE DISPERSION
C
        DIMENSION NU(NGM,NC), QAREA(NGM), TSTR(NGM), GSTR(NGM)
        DIMENSION GDOT(NGM), GDDOT(NGM), RKIJ(NGM,NGM), RKTIJ(NGM,NGM)
        DIMENSION ALFAIJ(NGM,NGM)
C
C DIMENSIONAMIENTO DE ARREGLOS RELACIONADOS CON EL TERMINO DE ASOCIACION
C
        DIMENSION ONUASOC(12,NC)
        DIMENSION MASSOC(12),ENASS(6,12,6,12), RKASS(6,12,6,12)
C
        DIMENSION FFI(NC)  !DIMENSIONAMIENTO COEFICIENTES DE FUGACIDAD
C
	DIMENSION RNCTRL(NCM), FFIANAL(NCM)
C
	COMMON/DERIVNUMER/NOPCION
C
	IERRDER = 0
C
	IF(NOPCION.EQ.1) NCODRES = NCOD
	IF(NOPCION.EQ.1) NCOD = 5
C
C   EN ESTE PRIMER LLAMADO SE CALCULA TODO: Egia de Helmholtz, factor de
C     compresibilidad y coeficiente de fugacidad SOLO si NOPCION = 1
C      Si NOPCION = 0 se respeta el valor de NCOD original
C
        CALL ASSOCGCEOS3(NCOD,NC,DC,RN,TEMPCRIT,V,TEMP,NG,NU,QAREA
     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
     @,ENASS,RKASS,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,ierrassoc)
C
C FACTORES PARA CALCULO DE DERIVADAS EN FORMA NUMERICA
C
	FVOL = 1.D-8
	FMOL = 1.D-8
C
	IF(NOPCION.EQ.1) THEN
C		------------------------------------------
C
C		ALMACENAMIENTO DE INFORMACION ANALITICA
C
		VCTRL = V  !variable de entrada
C
C		ZSISANAL = FACTOR DE COPRES. DEL SISTEMA CALCULADO EN 
C		FORMA ANALITICA
C
		ZSISANAL = ZSIS    !variable de salida
		ARESCTRL = ARESTOT !variable de salida
C
		SUMAR = 0.
C
		DO I=1,NC
			RNCTRL(I)=RN(I)  !variable de entrada
C
C			FFIANAL(I)= COEF. DE FUGAC. CALCULADOS EN FORMA
C					ANALITICA
C
			FFIANAL(I)=FFI(I) !variable de salida
C
			SUMAR = SUMAR + RN(I)
		END DO
C
		RMOLTOT = SUMAR !Nro. total de moles
C
C		Perturbacion del Volumen
C
		DELTAV = FVOL * VCTRL 
		V = VCTRL + DELTAV 
C
		NCOD = 4 !calculo solo de energia residual
        CALL ASSOCGCEOS3(NCOD,NC,DC,RN,TEMPCRIT,V,TEMP,NG,NU,QAREA
     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
     @,ENASS,RKASS,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP,
     ,ierrassoc)
C
C     		DERIFV= DERIVADA DE ARESTOT CON RESPECTO AL VOLUMEN V
C
		DERIFV = (ARESTOT - ARESCTRL)/ DELTAV 
C
C		ZSISNUMER = Factor de compresib. del sistema calculado en forma
C				numerica
C	
		ZSISNUMER = (-1.) * (VCTRL/RMOLTOT) * DERIFV  + 1. 
C
C
		if(dabs(ZSISNUMER-ZSISANAL).eq.0) THEN
C
		  OPEN(UNIT=1,TYPE='OLD',NAME='ASSOCPARA.RES',
     @ FORM='FORMATTED',
     @ ACCESS='APPEND')
C
		  write(1,*) 
     @' WARNING: El paso para la dif. numerica resp. al vol 
     @es demasiado peq.'
C
		  CLOSE(UNIT=1)
		 ELSE
		  CONTINUE
		END IF
C
		IF(DABS((ZSISNUMER-ZSISANAL)/ZSISANAL).GT.1.D-4) THEN
C
		  OPEN(UNIT=1,TYPE='OLD',NAME='ASSOCPARA.RES',
     @ FORM='FORMATTED',
     @ ACCESS='APPEND')
C
		  WRITE(1,*) ' LA DERIVACION NUMERICA CON RESPECTO AL VOLUMEN
     @ NO COINCIDE CON LA DERIVACION ANALITICA'
		  WRITE(1,*) ' Zanal =',ZSISANAL, ' Znum = ',ZSISNUMER
C
		  CLOSE(UNIT=1)
C
C		EN ESTE CASO LAS POSIBILIDADES SON:
C
C		 - POSIBLE ERROR DE PROGRAMACION
C     		 - VALOR DEMASIADO ALTO DEL PASO. Esto es FVOL  muy alto
C		 - Valor demasiado bajo del paso. Esto es FVOL  muy bajo
C							para la precision de la
C							maquina
C
		   IERRDER=1
		 ELSE
		  CONTINUE
		END IF
C
C		RETORNO el volumen a su valor inicial
C
		V = VCTRL 
C
C		CALCULO DE COEFICIENTES DE FUGACIDAD EN FORMA NUMERICA
C
		DO I=1,NC
C			---------------------------------
C			Perturbacion del nro. de moles
C			Afecto solo al numero de moles del compuesto I
C
			DELTAM = FMOL * RN(I)
			RN(I) = RNCTRL(I) + DELTAM 
C
		        NCOD = 4 !calculo solo de energia residual
        CALL ASSOCGCEOS3(NCOD,NC,DC,RN,TEMPCRIT,V,TEMP,NG,NU,QAREA
     @,TSTR,GSTR,GDOT,GDDOT,RKIJ,RKTIJ,ALFAIJ,NGA,ONUASOC,MASSOC
     @,ENASS,RKASS,ARESTOT,ZSIS,FFI,IERRFUG,IERRZSIS,IERRREP
     @,ierrassoc)
C
			DERIFN =  (ARESTOT - ARESCTRL)/ DELTAM 
C
			IF(ZSISANAL.GT.0) THEN
C                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
			FFINUMER = DEXP( DERIFN - DLOG(ZSISANAL) ) 
C
			IF(DABS(FFINUMER-FFIANAL(I)).EQ.0) THEN
C
			 OPEN(UNIT=1,TYPE='OLD',NAME='ASSOCPARA.RES',
     @ FORM='FORMATTED',
     @ ACCESS='APPEND')
C
 			 WRITE(1,*) 
     @' WARNING: El paso para la diferenciacion numerica con respecto 
     @  al numero de moles es demasiado chico'
C
			 CLOSE(UNIT=1)
C
			 ELSE
			  CONTINUE
			END IF
C
			IF(DABS((FFINUMER-FFIANAL(I))/FFIANAL(I)).GT.1.D-4) THEN
C
			 OPEN(UNIT=1,TYPE='OLD',NAME='ASSOCPARA.RES',
     @ FORM='FORMATTED',
     @ ACCESS='APPEND')
C
		  	 WRITE(1,*) ' LA DERIVACION NUMERICA CON RESPECTO AL 
     @ NUMERO DE MOLES DEL COMP ',I,' NO COINCIDE CON LA
     @ DERIVACION ANALITICA'
		  	 WRITE(1,*) ' COEFanal = ',FFIANAL(I), 
     @' COEFnum = ', FFINUMER
C
			 CLOSE(UNIT=1)
C
			  IERRDER = 1
			 ELSE
			  CONTINUE
			END IF
C
C                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			 ELSE
			END IF
C
C
C			Recupero el valor del numero de moles
C
			RN(I) = RNCTRL(I) 
C 
C			---------------------------------
		END DO
C
C		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C		Recuperacion de los valores adecuados en las variables
C		 de salida que fueron perturbadas

C
		ARESTOT = ARESCTRL 
		ZSIS = ZSISANAL 
C
C
		DO I=1,NC
C
			FFI(I) = FFIANAL(I)
C
		END DO
C
C		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		NCOD = NCODRES 
C		------------------------------------------
C
	 ELSE
	  CONTINUE
	END IF
C
C
	END
C
C///////////////////////////////////////////////////////////////////////
