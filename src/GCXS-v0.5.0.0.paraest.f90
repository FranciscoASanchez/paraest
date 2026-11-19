!-----------------------------------------------------------------------------------------
!
!	Cálculo de entalpías en exceso de mezclas binarias. Usa directamente todo el archivo
!	GCA-vX.Y.Z.for aunque sólo llame a 2 subrutinas
!
!	Francisco S. 19/04/2011
!-----------------------------------------------------------------------------------------
!
subroutine GCXS(ifl,NG,NST,T,P,feed,propEx,CZ,iErrEx)
!
!	Declaraciones:
	implicit real*8 (A-H,o-Z)
	parameter(R = 82.05,MaxD=500,MaxNS=50)
	dimension xn(15),xn1(15),xn2(15),NP(MaxNS),Z(2),feed(15)
	character*10::version='0.4.0',AName
	character*20 versGCA
!	real*8 MatrizDatos(MaxD,4)
	character*50 DataName(MaxNS)
!
!	COMMONs:
	common/versSUB/versGCA
	common/NAME/AName(15)
	common/GCPROM/PMM(15),Pen(15),HHA(15),HHB(15),HHC(15),HHD(15),HHE(15),HHF(15)
!
!
!
!
!	open(unit=5,file='entrada.txt',status='OLD')
!	open(unit=2,file='salida.txt')
!
!	Impresión de información sobre versiones
1	format(/,74("-"),/)
!	write(2,1)
!	write(2,'("GCXS version ",A10)')version
!	write(2,'("Calculo de entalpia y volumen de exceso de una mezcla binaria con la GCA-EoS")')
!	write(2,*)'Version del paquete termodinamico: ',versGCA
!	write(2,1)
!
!	La variable it hace referencia al tipo de entrada que se va a ingresar.
!	Si it = 0, el programa lee T,P del sistema, y se limita a calcular NP pun-
!	tos separados de forma homogénea.
!	Si it = 1, el usuario luego de especificar NP, debe incluir ND datos de x(1,j)
!	HE(j) a T y P
!
!	La variable iunid hace referencia a las unidades. iUnid = 0 usa las unidades
!	por defecto en la GCA-EoS (R = 82.05 atm cm3/(mol K)), [Hr] = atm cm3/mol.
!	iunid = 1, usa unidades SI: [Hr] = J/mol.
!
!	NSet hace referencia al nro de sets de datos que se va a ingresar. La varia-
!	ble es útil sólo si it = 1. Sino, cualquier valor puede ser ingresado.
!
!	read(5,*)NG,NGA,it,iUnid,NSet
	iUnid = 1
	NSet = 1
	it = 1
	NC = 2

!	Lectura de parámetros grupales y de moleculares
!	CALL PARGR(NC,NG,NGA,NST,2)
!	CALL PARMOL(NC,NG,NST,2)
!	write(2,100)
!
!	Inicio para el cálculo del volumen (subrutina call ZMAX)
	IGZ = 0
!	Se requieren propiedades derivadas respecto de T para calcular entalpías:
	if(ifl.EQ.2)NTEMP = 1
!
!
	if(it.EQ.0)then
!
!       ----------------------------------------------------------------------------
!	  it = 0: Cálculo de NSet curvas de HE y VE, de NPt puntos cada una
	  do io = 1,NSet
!
	    read(5,*)NPt
	    read(5,*)T,P
!
!	    Parámetros dependientes de T
	    call PARAGC(T,NC,NG,NST,NTEMP)
!	    Compuesto 1 puro a T y P
	    xn1(:NC) = (/1D0, 1D-12/)
!
!	    Entalpía residual del mismo
	    call NTLPY(NC,NG,NST,T,P,xn1(:NC),HRES1,CZ,IGZ,1,IC)
!	    Cambio de unidades acorde:
	    if(iUnid.EQ.1)then

	      write(2,'(2X,"Entalpías residuales de los puros (J/mol):")')
	      HRes1 = HRes1*.101325D0

	    elseif(iunid.EQ.0)then

	      write(2,'(2X,"Entalpías residuales de los puros (atm cm3/mol):")')

	    else

	      call error(1)
	      stop

	    endif
	    write(2,2)
2	    format(59("-"))
	    write(2,100)AName(1),HRes1
100	    format(A10,X,F10.2)
!
!	    Volumen molar de 1. La traslación de Peneloux se aplica en la impresión de valores
	    V1 = CZ*R*T/P
!
!	    Compuesto 2 puro a T y P
	    xn2(:NC) = (/1D-12, 1D0/)
!	    Entalpía residual del mismo:
	    call NTLPY(NC,NG,NST,T,P,xn2(:NC),HRES2,CZ,IGZ,1,IC)
	    if(iUnid.EQ.1)then

	      HRes2 = HRes2*.101325D0

	    endif
	    write(2,100)AName(2),HRes2
	    write(2,*)
!
!	    Volumen molar de 2:
	    V2 = CZ*R*T/P
	    xn(1) = 0D0
!
!	    La base es continuamente n = 1, por lo tanto, la E^R devuelta por la subrutina
!	    NTLPY es total (como todo el paquete termodinamico), pero = a la molar en magnitud.
	    if(iUnid.EQ.0)then

	      write(2,'(2X,"x(1)   V (cm3/mol)    V^E (cm3/mol)       H^E (atm cm3/mol)")')

	    elseif(iUnid.EQ.1)then

	      write(2,'(2X,"x(1)   V (cm3/mol)    V^E (cm3/mol)      H^E (J/mol)")')

	    endif
	    write(2,2)
!
!	    Escribe el primer punto comenzando por x1 = 0
	    write(2,10)xn(1),V2+Pen(2),0,0
10	    format(3X,F5.3,2X,F7.2,7X,F10.4,6X,F13.2)
!
!
!	    Comienzo lazo de entalpías residuales de la mezcla.
!	    NP = 50
	    do i = 1, NPt - 1

	      xn(1) = xn(1) + 1D0/NPt
	      xn(2) = 1 - xn(1)

!	      Posible chequeo de existencia de vapor:
!	      call PSAT(T,PS,Z,NC,NG,nga,IC,IGP,IGZ,IER)
!
!	      Llamado a HR, siempre como fase líquida
	      call NTLPY(NC,NG,NST,T,P,xn(:NC),HRES,CZ,IGZ,0,IC)
!
	      V = CZ*R*T/P
	      if(iUnid.EQ.1)HRes = HRes*.101325D0
!
!	      Entalpía y volumen en exceso:
	      HE = HRES - xn(1)*HRES1 - xn(2)*HRES2
	      VE = V - xn(1)*V1 - xn(2)*V2
!
!	      El volumen molar se imprime corregido con la traslación de Peneloux de
!	      la mezcla. El VE es invariante.
	      write(2,10)xn(1),V+xn(1)*Pen(1)+xn(2)*Pen(2),VE,HE

	    enddo
	    xn(1) = 1D0
!	    Impresión del último punto
	    write(2,10)Xn(1),V1+Pen(1),0,0
	    write(2,'(40("-"),//)')
!
!	    Fin lazo de sets si IT = 0
!         -----------------------------------------------------------
	  enddo
!
!	Condicional: it = 1, implica cálculo de desviación de datos experimentales.
	elseif(it.EQ.1)then
!
!	  call lectura(NSet,DataName,NP,MatrizDatos)
!	  NPT = 0
	  NPT = 1
!
!	  Lazo para el set i
	  do i = 1,NSet
!
!	    read(5,*)NP(i)
	    NP(i) = 1
!	    NPT = NP(i) + NPT
!
!	    read(5,*)DataName(i)
!
!	    Identificación del set de datos alimentado:
!	    write(2,'(2X,"Set de datos nro",X,I2,2X,A50,/)')i,DataName(i)
!
!	    if(iunid.EQ.0)then

!	      write(2,'(3X,"T(K)    P(atm)     x(1)   V^E(cm3/mol)  Hexp^E(atm cm3/mol)  H^E (atm cm3/mol)  Error rel. %")')

!	    elseif(iunid.EQ.1)then

!	      write(2,'(3X,"T(K)    P(atm)     x(1)   V^E(cm3/mol)    Hexp^E(J/mol)  H^E(J/mol)  Error rel. %")')

!	    else

!	      call error(2)

!	    endif
!	    write(2,1)
!
!	    aerrorH = ARD% en entalpía en exceso
!	    aErrorH = 0D0
!	    nError = 0
!
!	    Lazo interno para el punto NP(i) del set i
	    do j = 1,NP(i)
!
!	      read(5,*)T,P,x1exp,HEexp
	      call PARAGC(T,NC,NG,NST,NTEMP)
!
!	      Compuesto 1 puro a T y P
	      xn(:NC) = (/1D0, 1D-12/)
	      IF(ifl.EQ.1)then

!	        Cálculo sólo de volumen en exceso: Z1:
	        call ZMAX(Z,0,T,P,XN,iGz,NC,NG,NST,iC)
	        if(iC.EQ.-1)then

!	          Se sale de la subrutina si hay error.
	          iErrEx = 1
	          return

	        endif
	        Z1 = Z(2)
		  puro1 = Z1

	      else

!	        Cálculo de entalpía en exceso
	        call NTLPY(NC,NG,NST,T,P,xn(:NC),puro1,CZ,IGZ,0,IC)
	        if(iC.EQ.-1)then

!	          Se sale de la subrutina si hay error.
	          iErrEx = 1
	          return

	        endif
	        if(iunid.EQ.1)then

	          puro1 = puro1*.101325D0

	        elseif(iunid.NE.0)then

	          call error(1)
	          stop

	        endif
!	        Volumen molar de "1"
!	        V1 = CZ*R*T/P
	        Z1 = cZ

	      endif
!
!
!	      Compuesto 2 puro a T y P
	      xn(:NC) = (/1D-12, 1D0/)
	      IF(ifl.EQ.1)then

!	        Cálculo sólo de volumen en exceso: Z2:
	        call ZMAX(Z,0,T,P,XN,iGz,NC,NG,NST,iC)
	        if(iC.EQ.-1)then

!	          Se sale de la subrutina si hay error.
	          iErrEx = 1
	          return

	        endif
	        puro2 = Z(2)
		  Z2 = puro2

	      else

!	        Cálculo de entalpía en exceso
	        call NTLPY(NC,NG,NST,T,P,xn(:NC),puro2,CZ,IGZ,0,IC)
	        if(iC.EQ.-1)then

!	          Se sale de la subrutina si hay error.
	          iErrEx = 1
	          return

	        endif
	        if(iunid.EQ.1)then

	          puro2 = puro2*.101325D0

	        elseif(iunid.NE.0)then

	          call error(1)
	          stop

	        endif
!	        Factor de compresibilidad de "1" para inicio de cálculo de Z_mezcla
	        Z2 = CZ

	      endif
!	      Volumen molar de "2":
!	      V2 = CZ*R*T/P
!
!	      Cálculo para el punto NP del set i
!	      No se permiten x = 0:
!           modError aparece para cuando se ingresa una fracción molar = 0 ó 1.
!	      nError cuenta las veces uque aparece este error en un set: nError = 0, 1 ó 2.
!	      modError = 0
!	      if(x1exp.LT.1D-20)then
!
!	        Están ingresándose compuestos puros
!	        x1exp = 1D-12
!	        modError = 1
!	        nError = nError + 1

!	      elseif(x1exp.GE.1D0)then

!	        x1exp = 0.999999999990D0
!	        modError = 1
!	        nError = nError + 1
!
!	      endif
!
!	      Alimentación:
	      xn(:NC) = feed(:NC)
	      Z(2) = xn(1)*Z1 + xn(2)*Z2
	      Z(1) = 1
	      if(ifl.EQ.1)then

	        call ZMAX(Z,0,T,P,XN,0,NC,NG,NST,iC)
	        if(iC.EQ.-1)then

!	          Se sale de la subrutina si hay error.
	          iErrEx = 1
	          return

	        endif
	        propMezcla = Z(2)
	        CZ = Z(2)

	      elseif(ifl.EQ.2)then

	        call NTLPY(NC,NG,NST,T,P,xn(:NC),propMezcla,Z(2),0,0,IC)
	        if(iC.EQ.-1)then

!	          Se sale de la subrutina si hay error.
	          iErrEx = 1
	          return

	        endif
!	        V = CZ*R*T/P
	        CZ = Z(2)
	        if(iUnid.EQ.1)propMezcla = propMezcla*.101325D0

	      endif
!
!
!
!	      Propiedad en exceso (V o H)
	      propEx = propMezcla - xn(1)*puro1 - xn(2)*puro2
!
!	      Cálculo de volumen si ifl = 1:
	      if(ifl.EQ.1)propEx = propEx*R*T/P
!	      VE = V - xn(1)*V1 - xn(2)*V2
!
!	      if(modError.EQ.0)then

!	        Cálculo del error:
!	        errorH = 1D0 - HE/HEexp
!	        aerrorH = abs(errorH) + aerrorH
!	        if(IUnid.EQ.0)then

!	          write(2,211)T,P,xn(1),VE,HEexp,HE,errorH*1D2
211	          format(X,F7.2 ,X, F9.4, 3X, F5.3, 4X, F8.3, 3X, 2(F13.2,10X), F8.3)

!	        else

!	          write(2,213)T,P,xn(1),VE,HEexp,HE,errorH*1D2
213	          format(X, F7.2, X, F9.4, 3X, F5.3, 4X, F8.3, 3X, 2(F13.2,X), 5X, F8.3)

!	        endif

!	      else
!	        if(IUnid.EQ.0)then

!	          write(2,212)T,P,xn(1),VE,HEexp,HE
212	          format(X, F7.2, X, F9.4, 3X, F5.3, 4X, F8.3, 3X, 2(F13.2,10X), "------")

!	        else

!	          write(2,214)T,P,xn(1),VE,HEexp,HE
214	          format(X, F7.2, X, F9.4, 3X, F5.3, 4X, F8.3, 3X, 2(F13.2,X), 5X, "------")

!	        endif
!	      endif
!
!	    Fin del lazo del punto NP(i)
	    enddo
!	    write(2,210)aerrorH*1D2/(NP(i) - nError)
210	    format(/,2X,23HError relativo en H^E =,F8.3,X,"%",//)
!	    if(modError.EQ.1)then

!	      write(2,'("Atención: tenga en cuenta que no se consideraron los puntos correspondientes a",/, &
!&	                "compuestos puros, es decir, x = 0 ó 1",/)')

!	    endif
!
!	  Fin del lazo del set i
	  enddo
!
	else

	  call error(2)

	endif
!
	return
!	stop
!
101	format(10H H^R(2) = ,F10.2,9H J mol^-1)
102	format(10H H^R(1) = ,F10.2,16H atm cm^3 mol^-1)
103	format(10H H^R(2) = ,F10.2,16H atm cm^3 mol^-1)
!
end
!	----------------------------------------------------------------------------
!
!
!
subroutine error(tipo)
!
!	Esta subrutina informa el tipo de error en el archivo de entrada, cuando es-
!	tá disponible.
!
!	Planeo ir agregándolos a medida que me vaya dando cuenta cómo y cuales con-
!	viene agregar. En general se trata de errores en el ingreso de parámetros
!	pero pueden haber también errores en la convergencia de la subrutina termo-
!	dinámica.
!				Francisco, Junio de 2011. v0.2.0
!
	integer tipo
!
	write(*,'(//)')
	if(tipo.EQ.1)then
!
	  write(*,'("Error de lectura en la primera línea: usted ha in-",/, &
&	            "gresado un tipo de unidad incorrecto",//,              &
&	            "La variable mal especificada es IUnid. Esta varia-",/  &
&	            "ble se ingresa 4º en la 1ª línea:",/)')
	  write(*,'("NG   NGA   IT   IUnid   NSet",//)')
	  write(*,'("Los tipos de unidades soportados son",/,               &
&	            "IUnid = 0: por defecto en la GCA:",/)')
	  write(*,'(25X,"[energía] = atm cm3/mol",//)')
	  write(*,'("IUnid = 1: común a la mayoría de las publicaciones ",/,&
&	            "de datos experimentales de entalpías en exceso:",//)')
	  write(*,'(25X,"[energía] = J/mol")')
!
	elseif(tipo.EQ.2)then
!
	  write(*,'("Error de lectura en la primera línea: usted ha in-",/, &
&	            "gresado un tipo de cálculo incorrecto.",/)')
	  write(*,'("La variable mal especificada es IT. Esta variable se ingresa 3º en la 1ª línea:",/)')
	  write(*,'("NG   NGA   IT   IUnid   NSet",//)')
	  write(*,'("Los tipos de cálculo se refieren a:",///)')
!
	  write(*,'("IT = 0: calcular NP puntos de VE y HE de la mezcla, ",/&
&	            "separados una distancia x(1) = 1/NP",///)')
	  write(*,'("IT = 1: dados NSet set de dates esperimentales, cal-",/&
&	            "cular la HE y VE predichos dados T, P y x emperimentales",/)')

	endif
!
	return
!
endsubroutine
!	----------------------------------------------------------------------------
!
subroutine lectura(MaxD,NSet,DataName,NP,MatrizDatos)
!
!	Esta subrutina tiene el propósito de leer NSet de datos experimentales cada
!	uno con NP datos.
!
!	Actualmente no está en uso. Francisco, Junio de 2011. Durante versión 0.2.0
!
	real*8 MatrizDatos(MaxD,4)
	character*50 DataName(NSet)
	dimension NP(NSet)
!
!
	NPT = 0
	do i = 1,NSet
!
	  read(1,*)NP(i)
	  NPT = NP(i) + NPT
	  read(1,*)DataName(i)
	  do j = 1,NP(i)
!
	    read(1,*)MatrizDatos(i,1),MatrizDatos(i,2),MatrizDatos(i,3),MatrizDatos(i,4)
!
	  enddo
!
	enddo
!
	return

endsubroutine