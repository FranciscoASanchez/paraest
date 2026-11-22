!--------------------------------------------------------------------------------------------------------------------------
!
!  Module gca_parameters: structure and parameters of the GCA-EoS model, and
!  the mapping between the fitted parameters and their place in the model.
!
!  Transitional module: when paraest uses the GCA package (with its
!  parameters in modules), this is replaced by the package's module.
!
!--------------------------------------------------------------------------------------------------------------------------
module gca_parameters
  use dimensions, only: NCM, NGM, NPTIPOM, NIDM, NPARM
  implicit none

  ! --- Structure ----------------------------------------------------------
  integer :: NG, NGA                       ! groups and associating groups
  integer :: NU(NGM,NCM)                   ! groups of each component
  integer :: MASSOC(12)                    ! sites of each associating group
  double precision :: ONUASOC(12,NCM)      ! associating groups of each component

  ! --- Dispersive parameters (per group and per pair of groups) -------------
  double precision :: QAREA(NGM), TSTR(NGM), GSTR(NGM), GDOT(NGM), GDDOT(NGM), TSPLINE(NGM)
  double precision :: RKIJ(NGM,NGM), RKTIJ(NGM,NGM), ALFAIJ(NGM,NGM)

  ! --- Association parameters ---------------------------------------------
  double precision :: ENASS(6,12,6,12), RKASS(6,12,6,12)

  ! --- Critical diameter --------------------------------------------------
  double precision :: DC(NCM)
  ! Vapor pressure point to which dc is fitted
  double precision :: TSAT(NCM), PSAT(NCM)

  ! --- Fitted parameters --------------------------------------------------
  ! Number of each type
  integer :: NNP, NNQ, NALFA, NENASS, NRKASS, numTc, numDC
  ! g_ii / k_ij: indices, how many identical ones each has and which they are
  integer :: IDG(NPTIPOM,2), ITEMP(NPTIPOM), NIDG(NPTIPOM), IDIDG(NIDM,2)
  integer :: IDQ(NPTIPOM)                                                     ! areas
  integer :: IDALF(NPTIPOM,2), NIDALF(NPTIPOM), IDIDALF(NIDM,2)               ! alpha
  integer :: IDENASS(NPTIPOM,4), NIDENASS(NPTIPOM), IDIDENASS(NIDM,4)         ! assoc. energy
  integer :: IDRKASS(NPTIPOM,4), NIDRKASS(NPTIPOM), IDIDRKASS(NIDM,4)         ! assoc. volume
  integer :: idTc(NPTIPOM), IDDC(NPTIPOM)                                     ! Tc and dc

contains

!------------------------------------------------------------------------------------------------------------------------------------
!
  ! Sets the model parameters from the new values of the fitted parameters,
  ! which come unscaled in XXAUX. tempCrit (the critical temperatures of the
  ! components) changes only when Tc is fitted.
!
SUBROUTINE set_model_params(XXAUX, tempCrit)
!
  integer :: I, J, L, NPAR2, NSSA, NSSE, NSSRK, NSSS
!
  double precision, intent(in) :: XXAUX(NPARM)       ! values of the fitted parameters
  double precision, intent(inout) :: tempCrit(NCM)   ! critical temperatures
!
  NPAR2 = 0
!
  ! Parameters of type gii or kij
  IF(NNP.NE.0) then !A

    DO I=1,NNP

      IF(IDG(I,1).EQ.IDG(I,2)) THEN

        GSTR(IDG(I,1)) = XXAUX(NPAR2+1)
        NPAR2 = NPAR2 + 1
        IF(ITEMP(I).EQ.0)then

          GDOT(IDG(I,1)) = 0.D0
          GDDOT(IDG(I,1)) = 0.D0

        endif
        IF(ITEMP(I).GE.1) GDOT(IDG(I,1)) = XXAUX(NPAR2+1)
        IF(ITEMP(I).GE.1) NPAR2 = NPAR2 + 1
        IF(ITEMP(I).EQ.1) GDDOT(IDG(I,1)) = 0.D0
        IF(ITEMP(I).EQ.2) GDDOT(IDG(I,1)) = XXAUX(NPAR2+1)
        IF(ITEMP(I).EQ.2) NPAR2 = NPAR2 + 1
!
      ELSE

        RKIJ(IDG(I,1),IDG(I,2)) = XXAUX(NPAR2+1)
        NPAR2 = NPAR2 + 1
        IF(ITEMP(I).EQ.0) RKTIJ(IDG(I,1),IDG(I,2))= 0.D0
1001        IF(ITEMP(I).EQ.1) RKTIJ(IDG(I,1),IDG(I,2))=XXAUX(NPAR2+1)
        IF(ITEMP(I).EQ.1) NPAR2 = NPAR2 + 1
        IF(ITEMP(I).EQ.2) THEN

          ! Avoid human error
          ITEMP(I) = 1
          GO TO 1001

        END IF
!
      END IF

    END DO

  endif
!
  ! Identical ones of type gii or kij
  IF(NNP.NE.0) THEN

    NSSS = 0
    DO I=1,NNP

      if(NIDG(I).ne.0) then

        IF(IDG(I,1).EQ.IDG(I,2)) THEN

          DO J=(1+NSSS),(NIDG(I)+NSSS)

            GSTR(IDIDG(J,1))  = GSTR(IDG(I,1))
            GDOT(IDIDG(J,1))  = GDOT(IDG(I,1))
            GDDOT(IDIDG(J,1)) = GDDOT(IDG(I,1))

           END DO
          NSSS = NSSS + NIDG(I)

        ELSE

          DO J=(1+NSSS),(NIDG(I)+NSSS)

            RKIJ(IDIDG(J,1),IDIDG(J,2)) = RKIJ(IDG(I,1),IDG(I,2))
            RKTIJ(IDIDG(J,1),IDIDG(J,2))= RKTIJ(IDG(I,1),IDG(I,2))

          END DO
          NSSS = NSSS + NIDG(I)
!
        END IF

      end if

    END DO

  END IF
!
  ! Area parameters
  IF(NNQ.NE.0) THEN

    DO I=1,NNQ

      QAREA(IDQ(I)) = XXAUX(NPAR2+1)
      NPAR2 = NPAR2 + 1

    END DO

  END IF
!
  ! Alpha parameters
  IF(NALFA.NE.0) THEN

    DO I=1,NALFA

      ALFAIJ(IDALF(I,1),IDALF(I,2)) = XXAUX(NPAR2+1)
      NPAR2 = NPAR2 + 1

    END DO

  END IF
!
  ! Identical alphas
  IF(NALFA.NE.0) THEN

    NSSA = 0
    DO I=1,NALFA

      if(NIDALF(I).ne.0) then

        DO J=(1+NSSA),(NIDALF(I)+NSSA)

          ALFAIJ(IDIDALF(J,1),IDIDALF(J,2))= ALFAIJ(IDALF(I,1),IDALF(I,2))

        END DO
        NSSA = NSSA + NIDALF(I)

      end if

    END DO

  END IF
!
  ! Association energy parameters
  IF(NENASS.NE.0) THEN
!
    DO I=1,NENASS

      ENASS(IDENASS(I,1),IDENASS(I,2),IDENASS(I,3),IDENASS(I,4)) = XXAUX(NPAR2+1)
      NPAR2 = NPAR2 + 1

    END DO

  END IF
!
  ! Identical association energies
  IF(NENASS.NE.0) THEN

    NSSE = 0
    DO L=1,NENASS

      if(NIDENASS(L).ne.0) then

        DO J=(1+NSSE),(NIDENASS(L)+NSSE)

          ENASS(IDIDENASS(J,1),IDIDENASS(J,2),IDIDENASS(J,3),IDIDENASS(J,4)) = ENASS(IDENASS(L,1),IDENASS(L,2) ,IDENASS(L,3),IDENASS(L,4))

        END DO
        NSSE = NSSE + NIDENASS(L)

      end if

    END DO

  END IF
!
  ! Association volume parameters
  if(NRKASS.ne.0) then

    DO I=1,NRKASS

      RKASS(IDRKASS(I,1),IDRKASS(I,2),IDRKASS(I,3),IDRKASS(I,4)) = XXAUX(NPAR2+1)
      NPAR2 = NPAR2 + 1

    END DO

  end if
!
  ! Identical association volumes
  IF(NRKASS.ne.0) THEN

    NSSRK = 0
    DO L=1,NRKASS

      if(NIDRKASS(L).ne.0) then

        DO J=(1+NSSRK),(NIDRKASS(L)+NSSRK)

          RKASS(IDIDRKASS(J,1),IDIDRKASS(J,2),IDIDRKASS(J,3),IDIDRKASS(J,4)) = RKASS(IDRKASS(L,1),IDRKASS(L,2),IDRKASS(L,3),IDRKASS(L,4))

        END DO
        NSSRK = NSSRK + NIDRKASS(L)

      end if

    END DO

  END IF
!
  ! Critical temperatures and diameters, when they are fitted parameters
  IF(numtc.NE.0) THEN

    DO I=1,numtc

      tempCrit(idtc(I)) = XXAUX(NPAR2+1)
      NPAR2 = NPAR2 + 1

    END DO

  END IF

  IF(NUMDC.NE.0) THEN

    DO I=1,NUMDC

      DC(IDDC(I)) = XXAUX(NPAR2+1)
      NPAR2 = NPAR2 + 1

    END DO

  END IF
!
ENDsubroutine set_model_params

endmodule gca_parameters
