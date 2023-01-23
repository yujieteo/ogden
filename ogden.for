      SUBROUTINE UHYPER(BI1,BI2,AJ,U,UI1,UI2,UI3,TEMP,NOEL,
     1 CMNAME,INCMPFLAG,NUMSTATEV,STATEV,NUMFIELDV,FIELDV,
     2 FIELDVINC,NUMPROPS,PROPS)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION U(2),UI1(3),UI2(6),UI3(6),STATEV(*),FIELDV(*),
     1 FIELDVINC(*),PROPS(*)
      
C     Ogden parameters
      REAL(8) mu, alpha
      
C     Fitting parameters.
      REAL(8) r, c0, cr, cu, beta
      REAL(8) nu1, nu2, nu3, lambda0
      
C     Increment.
      REAL(8) N, Nmax
      
C     Calculated properties:
      REAL(8) kappa, gamma
      
C     Ogden properties
C     This is arbitrary.
C     mu      = 3.97  (kPa)
C     alpha   = 2.158 (-)
      mu      = PROPS(1)
      alpha   = PROPS(2)

C     Fitting parameters
C     How to obtain r, c0, beta?
C     First, use 1/r * ERF(Wmax - W(lambda1, lambda2)/c0 + beta * Wmax))
C     This is done as a separate curve fitting code, possible in MATLAB.
C     Then this is keyed into properties.
C     r       = 1.75
C     c0      = 10
C     beta    = 0.1
      r       = PROPS(3)
      c0      = PROPS(4)
      beta    = PROPS(5)
      
C     List of auxiliary and derived properties.
C     How to obtain nu1 and nu3?
C     
      nu1     = PROPS(6)
      nu3     = PROPS(7)

C     Do not understand how they calculate nu2? 
C     Is there an anisotropic formula.
      nu2     = nu1
      
C     Lambda-max for experiment
C     Lambda-max determines directions
      
      lambda0 = PROPS(8)
      lambda1 = lambda0 
      lambda2 = 1
      lambda3 = 1/lambda0
      
C     Kappa is calculated from the first cycle.
C     Equation (10) is equal to kappa from the first cycle.
      kappa   = mu * (lambda0 ** (alpha - 1) - lambda0 ** (- alpha - 1)
      
	U(1)    = 
	U(2)    =   

C     Initialise partial derivatives as zeroes first.

      UI1     = 0
      UI2     = 0
      UI3     = 0

	UI1(1)  = 
	UI1(2)  = 
	UI1(3)  = 
	
	UI2(3)  =
      UI2(5)  = 
	UI3(4)  = 
	UI3(6)  = 

      RETURN
      END