MODULE lim_Sdom

USE mod_flux_x
USE mod_flux_y

IMPLICIT NONE
PUBLIC

CONTAINS



SUBROUTINE RK4_k1_L_ux(rho,p,alpha,gamma1,gamma2,dpl,u_k1,rho_k1,alpha_k1)
!------------------------------------------------!
!                                                !
! RK4 k1 step for left wave                      !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,gamma1,gamma2,dpl
DOUBLE PRECISION, INTENT(OUT) :: u_k1,rho_k1,alpha_k1
DOUBLE PRECISION :: gamma,a,b,p_k1

p_k1     = p
rho_k1   = rho
alpha_k1 = alpha
gamma    = 1D0/(alpha_k1/(gamma1-1D0)+(1D0-alpha_k1)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k1/rho_k1)
b        = alpha_k1*(1D0-alpha_k1)*(gamma1-gamma2)/((1D0-alpha_k1)*gamma1+alpha_k1*gamma2)
u_k1     = -dpl/(rho_k1*a)
alpha_k1 =  dpl*b/(rho_k1*a**(2D0))
rho_k1   =  dpl/(a**(2D0))

END SUBROUTINE RK4_k1_L_ux


SUBROUTINE RK4_k1_R_ux(rho,p,alpha,gamma1,gamma2,dpr,u_k1,rho_k1,alpha_k1)
!------------------------------------------------!
!                                                !
! RK4 k1 step for right wave                     !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,gamma1,gamma2,dpr
DOUBLE PRECISION, INTENT(OUT) :: u_k1,rho_k1,alpha_k1
DOUBLE PRECISION :: gamma,a,b,p_k1

p_k1     = p
alpha_k1 = alpha
rho_k1   = rho
gamma    = 1D0/(alpha_k1/(gamma1-1D0)+(1D0-alpha_k1)/(gamma2-1))+1D0
a        = sqrt(gamma*p_k1/rho_k1)
b        = alpha_k1*(1D0-alpha_k1)*(gamma1-gamma2)/((1D0-alpha_k1)*gamma1+alpha_k1*gamma2)
u_k1     =  dpr/(rho_k1*a)
alpha_k1 =  dpr*b/(rho_k1*a**(2D0))
rho_k1   =  dpr/(a**(2D0))

END SUBROUTINE RK4_k1_R_ux


SUBROUTINE RK4_k2_L_ux(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpl,u_k2,rho_k2,alpha_k2)
!------------------------------------------------!
!                                                !
! RK4 k2 step for left wave                      !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpl
DOUBLE PRECISION, INTENT(OUT) :: u_k2,rho_k2,alpha_k2
DOUBLE PRECISION :: gamma,a,b,p_k2

rho_k2   = rho+rho_k1/2D0
p_k2     = p+1D0/2D0*dpl
alpha_k2 = alpha+alpha_k1/2D0
gamma    = 1D0/(alpha_k2/(gamma1-1D0)+(1D0-alpha_k2)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k2/rho_k2)
b        = alpha_k2*(1-alpha_k2)*(gamma1-gamma2)/((1D0-alpha_k2)*gamma1+alpha_k2*gamma2)
u_k2     = -dpl/(rho_k2*a)
alpha_k2 =  dpl*b/(rho_k2*a**(2D0))
rho_k2   =  dpl/(a**(2D0))

END SUBROUTINE RK4_k2_L_ux


SUBROUTINE RK4_k2_R_ux(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpr,u_k2,rho_k2,alpha_k2)
!------------------------------------------------!
!                                                !
! RK4 k2 step for right wave                     !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpr
DOUBLE PRECISION, INTENT(OUT) :: u_k2,rho_k2,alpha_k2
DOUBLE PRECISION :: gamma,a,b,p_k2

p_k2     = p+1D0/2D0*dpr
alpha_k2 = alpha+alpha_k1/2D0
rho_k2   = rho+rho_k1/2D0
gamma    = 1D0/(alpha_k2/(gamma1-1D0)+(1D0-alpha_k2)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k2/rho_k2)
b        = alpha_k2*(1D0-alpha_k2)*(gamma1-gamma2)/((1D0-alpha_k2)*gamma1+alpha_k2*gamma2)
u_k2     =  dpr/(rho_k2*a)
alpha_k2 =  dpr*b/(rho_k2*a**(2D0))
rho_k2   =  dpr/(a**(2D0))

END SUBROUTINE RK4_k2_R_ux


SUBROUTINE RK4_k3_L_ux(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpl,u_k3,rho_k3,alpha_k3)
!------------------------------------------------!
!                                                !
! RK4 k3 step for left wave                      !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpl
DOUBLE PRECISION, INTENT(OUT) :: u_k3,rho_k3,alpha_k3
DOUBLE PRECISION :: gamma,a,b,p_k3

p_k3     = p+1D0/2D0*dpl
alpha_k3 = alpha+alpha_k2/2D0
rho_k3   = rho+rho_k2/2D0
gamma    = 1D0/(alpha_k3/(gamma1-1D0)+(1D0-alpha_k3)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k3/rho_k3)
b        = alpha_k3*(1D0-alpha_k3)*(gamma1-gamma2)/((1D0-alpha_k3)*gamma1+alpha_k3*gamma2)
u_k3     = -dpl/(rho_k3*a)
alpha_k3 =  dpl*b/(rho_k3*a**(2D0))
rho_k3   =  dpl/(a**(2D0))

END SUBROUTINE RK4_k3_L_ux


SUBROUTINE RK4_k3_R_ux(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpr,u_k3,rho_k3,alpha_k3)
!------------------------------------------------!
!                                                !
! RK4 k3 step for right wave                     !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpr
DOUBLE PRECISION, INTENT(OUT) :: u_k3,rho_k3,alpha_k3
DOUBLE PRECISION :: gamma,a,b,p_k3

p_k3     = p+1D0/2D0*dpr
alpha_k3 = alpha+alpha_k2/2D0
rho_k3   = rho+rho_k2/2D0
gamma    = 1D0/(alpha_k3/(gamma1-1D0)+(1D0-alpha_k3)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k3/rho_k3)
b        = alpha_k3*(1D0-alpha_k3)*(gamma1-gamma2)/((1D0-alpha_k3)*gamma1+alpha_k3*gamma2)
u_k3     =  dpr/(rho_k3*a)
alpha_k3 =  dpr*b/(rho_k3*a**(2D0))
rho_k3   =  dpr/(a**(2D0))

END SUBROUTINE RK4_k3_R_ux




SUBROUTINE RK4_k4_L_ux(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpl,u_k4,rho_k4,alpha_k4)
!------------------------------------------------!
!                                                !
! RK4 k4 step for left wave                      !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpl
DOUBLE PRECISION, INTENT(OUT) :: u_k4,rho_k4,alpha_k4
DOUBLE PRECISION :: gamma,a,b,p_k4

p_k4     = p+dpl
alpha_k4 = alpha+alpha_k3
rho_k4   = rho+rho_k3
gamma    = 1D0/(alpha_k4/(gamma1-1D0)+(1D0-alpha_k4)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k4/rho_k4)
b        = alpha_k4*(1D0-alpha_k4)*(gamma1-gamma2)/((1D0-alpha_k4)*gamma1+alpha_k4*gamma2)
u_k4     = -dpl/(rho_k4*a)
alpha_k4 =  dpl*b/(rho_k4*a**(2D0))
rho_k4   =  dpl/(a**(2D0))

END SUBROUTINE RK4_k4_L_ux



SUBROUTINE RK4_k4_R_ux(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpr,u_k4,rho_k4,alpha_k4)
!------------------------------------------------!
!                                                !
! RK4 k4 step for right wave                     !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpr
DOUBLE PRECISION, INTENT(OUT) :: u_k4,rho_k4,alpha_k4
DOUBLE PRECISION :: gamma,a,b,p_k4

p_k4     = p+dpr
alpha_k4 = alpha+alpha_k3
rho_k4   = rho+rho_k3
gamma    = 1D0/(alpha_k4/(gamma1-1D0)+(1D0-alpha_k4)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k4/rho_k4)
b        = alpha_k4*(1D0-alpha_k4)*(gamma1-gamma2)/((1D0-alpha_k4)*gamma1+alpha_k4*gamma2)
u_k4     =  dpr/(rho_k4*a)
alpha_k4 =  dpr*b/(rho_k4*a**(2D0))
rho_k4   =  dpr/(a**(2D0))

END SUBROUTINE RK4_k4_R_ux



SUBROUTINE RK4_k1_L_vy(rho,p,alpha,gamma1,gamma2,dpl,v_k1,rho_k1,alpha_k1)
!------------------------------------------------!
!                                                !
! RK4 k1 step for left wave                      !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,gamma1,gamma2,dpl
DOUBLE PRECISION, INTENT(OUT) :: v_k1,rho_k1,alpha_k1
DOUBLE PRECISION :: gamma,a,b,p_k1

p_k1     = p
rho_k1   = rho
alpha_k1 = alpha
gamma    = 1D0/(alpha_k1/(gamma1-1D0)+(1D0-alpha_k1)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k1/rho_k1)
b        = alpha_k1*(1D0-alpha_k1)*(gamma1-gamma2)/((1D0-alpha_k1)*gamma1+alpha_k1*gamma2)
v_k1     = -dpl/(rho_k1*a)
alpha_k1 =  dpl*b/(rho_k1*a**(2D0))
rho_k1   =  dpl/(a**(2D0))

END SUBROUTINE RK4_k1_L_vy


SUBROUTINE RK4_k1_R_vy(rho,p,alpha,gamma1,gamma2,dpr,v_k1,rho_k1,alpha_k1)
!------------------------------------------------!
!                                                !
! RK4 k1 step for right wave                     !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,gamma1,gamma2,dpr
DOUBLE PRECISION, INTENT(OUT) :: v_k1,rho_k1,alpha_k1
DOUBLE PRECISION :: gamma,a,b,p_k1

p_k1     = p
alpha_k1 = alpha
rho_k1   = rho
gamma    = 1D0/(alpha_k1/(gamma1-1D0)+(1D0-alpha_k1)/(gamma2-1))+1D0
a        = sqrt(gamma*p_k1/rho_k1)
b        = alpha_k1*(1D0-alpha_k1)*(gamma1-gamma2)/((1D0-alpha_k1)*gamma1+alpha_k1*gamma2)
v_k1     =  dpr/(rho_k1*a)
alpha_k1 =  dpr*b/(rho_k1*a**(2D0))
rho_k1   =  dpr/(a**(2D0))

END SUBROUTINE RK4_k1_R_vy


SUBROUTINE RK4_k2_L_vy(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpl,v_k2,rho_k2,alpha_k2)
!------------------------------------------------!
!                                                !
! RK4 k2 step for left wave                      !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpl
DOUBLE PRECISION, INTENT(OUT) :: v_k2,rho_k2,alpha_k2
DOUBLE PRECISION :: gamma,a,b,p_k2

rho_k2   = rho+rho_k1/2D0
p_k2     = p+1D0/2D0*dpl
alpha_k2 = alpha+alpha_k1/2D0
gamma    = 1D0/(alpha_k2/(gamma1-1D0)+(1D0-alpha_k2)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k2/rho_k2)
b        = alpha_k2*(1-alpha_k2)*(gamma1-gamma2)/((1D0-alpha_k2)*gamma1+alpha_k2*gamma2)
v_k2     = -dpl/(rho_k2*a)
alpha_k2 =  dpl*b/(rho_k2*a**(2D0))
rho_k2   =  dpl/(a**(2D0))

END SUBROUTINE RK4_k2_L_vy


SUBROUTINE RK4_k2_R_vy(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpr,v_k2,rho_k2,alpha_k2)
!------------------------------------------------!
!                                                !
! RK4 k2 step for right wave                     !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpr
DOUBLE PRECISION, INTENT(OUT) :: v_k2,rho_k2,alpha_k2
DOUBLE PRECISION :: gamma,a,b,p_k2

p_k2     = p+1D0/2D0*dpr
alpha_k2 = alpha+alpha_k1/2D0
rho_k2   = rho+rho_k1/2D0
gamma    = 1D0/(alpha_k2/(gamma1-1D0)+(1D0-alpha_k2)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k2/rho_k2)
b        = alpha_k2*(1D0-alpha_k2)*(gamma1-gamma2)/((1D0-alpha_k2)*gamma1+alpha_k2*gamma2)
v_k2     =  dpr/(rho_k2*a)
alpha_k2 =  dpr*b/(rho_k2*a**(2D0))
rho_k2   =  dpr/(a**(2D0))

END SUBROUTINE RK4_k2_R_vy


SUBROUTINE RK4_k3_L_vy(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpl,v_k3,rho_k3,alpha_k3)
!------------------------------------------------!
!                                                !
! RK4 k3 step for left wave                      !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpl
DOUBLE PRECISION, INTENT(OUT) :: v_k3,rho_k3,alpha_k3
DOUBLE PRECISION :: gamma,a,b,p_k3

p_k3     = p+1D0/2D0*dpl
alpha_k3 = alpha+alpha_k2/2D0
rho_k3   = rho+rho_k2/2D0
gamma    = 1D0/(alpha_k3/(gamma1-1D0)+(1D0-alpha_k3)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k3/rho_k3)
b        = alpha_k3*(1D0-alpha_k3)*(gamma1-gamma2)/((1D0-alpha_k3)*gamma1+alpha_k3*gamma2)
v_k3     = -dpl/(rho_k3*a)
alpha_k3 =  dpl*b/(rho_k3*a**(2D0))
rho_k3   =  dpl/(a**(2D0))

END SUBROUTINE RK4_k3_L_vy


SUBROUTINE RK4_k3_R_vy(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpr,v_k3,rho_k3,alpha_k3)
!------------------------------------------------!
!                                                !
! RK4 k3 step for right wave                     !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpr
DOUBLE PRECISION, INTENT(OUT) :: v_k3,rho_k3,alpha_k3
DOUBLE PRECISION :: gamma,a,b,p_k3

p_k3     = p+1D0/2D0*dpr
alpha_k3 = alpha+alpha_k2/2D0
rho_k3   = rho+rho_k2/2D0
gamma    = 1D0/(alpha_k3/(gamma1-1D0)+(1D0-alpha_k3)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k3/rho_k3)
b        = alpha_k3*(1D0-alpha_k3)*(gamma1-gamma2)/((1D0-alpha_k3)*gamma1+alpha_k3*gamma2)
v_k3     =  dpr/(rho_k3*a)
alpha_k3 =  dpr*b/(rho_k3*a**(2D0))
rho_k3   =  dpr/(a**(2D0))

END SUBROUTINE RK4_k3_R_vy




SUBROUTINE RK4_k4_L_vy(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpl,v_k4,rho_k4,alpha_k4)
!------------------------------------------------!
!                                                !
! RK4 k4 step for left wave                      !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpl
DOUBLE PRECISION, INTENT(OUT) :: v_k4,rho_k4,alpha_k4
DOUBLE PRECISION :: gamma,a,b,p_k4

p_k4     = p+dpl
alpha_k4 = alpha+alpha_k3
rho_k4   = rho+rho_k3
gamma    = 1D0/(alpha_k4/(gamma1-1D0)+(1D0-alpha_k4)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k4/rho_k4)
b        = alpha_k4*(1D0-alpha_k4)*(gamma1-gamma2)/((1D0-alpha_k4)*gamma1+alpha_k4*gamma2)
v_k4     = -dpl/(rho_k4*a)
alpha_k4 =  dpl*b/(rho_k4*a**(2D0))
rho_k4   =  dpl/(a**(2D0))

END SUBROUTINE RK4_k4_L_vy



SUBROUTINE RK4_k4_R_vy(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpr,v_k4,rho_k4,alpha_k4)
!------------------------------------------------!
!                                                !
! RK4 k4 step for right wave                     !
!                                                !
!------------------------------------------------!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpr
DOUBLE PRECISION, INTENT(OUT) :: v_k4,rho_k4,alpha_k4
DOUBLE PRECISION :: gamma,a,b,p_k4

p_k4     = p+dpr
alpha_k4 = alpha+alpha_k3
rho_k4   = rho+rho_k3
gamma    = 1D0/(alpha_k4/(gamma1-1D0)+(1D0-alpha_k4)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k4/rho_k4)
b        = alpha_k4*(1D0-alpha_k4)*(gamma1-gamma2)/((1D0-alpha_k4)*gamma1+alpha_k4*gamma2)
v_k4     =  dpr/(rho_k4*a)
alpha_k4 =  dpr*b/(rho_k4*a**(2D0))
rho_k4   =  dpr/(a**(2D0))

END SUBROUTINE RK4_k4_R_vy

END MODULE lim_Sdom
