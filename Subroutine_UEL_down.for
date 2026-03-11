C ======================================================================
C
C                     Flexo Electric UEL Subroutine
C     
C                                                         2023.10.30
C                        
C
C ---------------------------  PROPS  -----------------------------------
C
C   For UEL Element (Flexoelectric)
C     PROPS(1) = Elastic modulus (E)
C     PROPS(2) = Poison's ratio (nu)
C     PROPS(3) = Length scale of strain gradient (l)
C     PROPS(4) = Flexoelectric coefficient (f1)
C     PROPS(5) = Flexoelectric coefficient (f2)
C     PROPS(6) = Electric permisivity (k)
C
C   For UMAT Element (eummy element)
C     PROPS(1) = Young's modulus (E)
C     PROPS(2) = Poison's ratio (nu)
C
C ------------------  Solution Dependent Variables  ---------------------
C
C - UEL Element (Phase Field)     
C     SVARS(1) = Empty
C     SVARS(2) = Empty
C             ��
C
C - UMAT Element (dummy element)
C     STATEV(1)  = Empty
C     STATEV(2)  = Empty
C             ��
C    
C ---------------------- Module Global SDV -------------------------
C
C     N_ELEM - # of element 
C     N_INPT - # of integration point per element
C     N_SDV  - # of global variables
C
C  - Mechanical Element
C     globalSDV(i,j, ) >> i-th element j-th integration point variables
C
C     globalSDV(i,j,1)  = STRESS ��11
C     globalSDV(i,j,2)  = STRESS ��22
C     globalSDV(i,j,3)  = STRESS ��12
C     globalSDV(i,j,4)  = High order STRESS ��111
C     globalSDV(i,j,5)  = High order STRESS ��122
C     globalSDV(i,j,6)  = High order STRESS ��112
C     globalSDV(i,j,7)  = High order STRESS ��211
C     globalSDV(i,j,8)  = High order STRESS ��222
C     globalSDV(i,j,9)  = High order STRESS ��212
C     globalSDV(i,j,10) = Electric Displacement d1
C     globalSDV(i,j,11) = Electric Displacement d2
C                 ��
C 
C ======================================================================
      MODULE GLOBAL
      
C     Paramter ,  PARAMETER�� ������ �������� ���� �ȵ�. �ؿ��� �ߺ��Ǿ �ȵ�.
      INTEGER N_ELEM, N_INPT, N_SDV, N_offset
      PARAMETER( N_ELEM   = 3000 )
      PARAMETER( N_offset = 3000 )
      PARAMETER( N_INPT   = 9    )
      PARAMETER( N_SDV    = 50   )


C     global variable
      !DOUBLE PRECISION globalSDV(N_ELEM,N_INPT,N_SDV)
      ! Common Block�� ���� �ʾƵ� Module�� ���� ������ ������
      ! Increment�� �ٲ� ���� ���� ��� ����Ǿ� ����. �ʱ�ȭx
      
      END
      
C
C ======================================================================
C     UEL Subroutine for Phase Field Element
C ======================================================================
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
C
C     ==================================================================
C     UEL Element
C     ==================================================================

      CALL UEL_2D_FLEXO_Q59_Axis_sym(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     1     NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,
     2     DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)

      END


C
C ======================================================================
C     FlEXO UEL 2D Q59  (Axis-Symmetric)
C     Requirement Subroutine 
C      - MODULE GLOBAL
C      - SUBROUITNE GAUSS_QUADRATURE_3
C      - SUBROUTINE Q9_SHAPE_FUNCTIONS
C      - SUBROUTINE Q9_SHAPE_FUNCTIONS
C      - SUBROUTINE MatInversion
C ======================================================================  
      SUBROUTINE UEL_2D_FLEXO_Q59_Axis_sym(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     1     NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,
     2     DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C
      USE GLOBAL     
C
      INCLUDE 'ABA_PARAM.INC'

      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
C
C     ==================================================================
C     Declare variables and constants
C     ==================================================================
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 FIVE=5.D0,SIX=6.D0,SEVEN=7.D0,EIGHT=8.D0,NINE=9.D0,HALF=0.5D0)

      PARAMETER(PI=4.0*atan(1.0))
      DOUBLE PRECISION Axis_radius

      DOUBLE PRECISION STATE(59)
      DOUBLE PRECISION STATE_u(18), STATE_e(9)
      DOUBLE PRECISION STATE_n(16), STATE_a(16)
      DOUBLE PRECISION inc_u(18), inc_e(9) 
      DOUBLE PRECISION inc_n(16), inc_a(16) 

      DOUBLE PRECISION node9_xy(9,2), node4_xy(4,2)
      
      DOUBLE PRECISION GP4(4,3), GP9(9,3), GP16(16,3), GP25(25,3)
      DOUBLE PRECISION xi, eta, wght
      DOUBLE PRECISION N9(9), dN9dxi(2,9), dN9dx(2,9)
      DOUBLE PRECISION N4(4), dN4dxi(2,4), dN4dx(2,4)
      DOUBLE PRECISION dxidx(2,2), dxdxi(2,2), detJ4, detJ9
      
      DOUBLE PRECISION N_u(2,18), N_e(1,9) , N_n(4,16) , N_a(4,16)
      DOUBLE PRECISION B_u(4,18), B_n(10,16), B_e(2,9)
      
      DOUBLE PRECISION E, v, sl, f1, f2, permi
      DOUBLE PRECISION xlam,mu
      DOUBLE PRECISION D_uu(4,4), D_en(2,10), D_ee(2,2), D_nn(10,10), D_eu(2,4)
      
      DOUBLE PRECISION K_el(59,59)
      DOUBLE PRECISION K_uu(18,18), K_ee(9,9)  , K_nn(16,16)
      DOUBLE PRECISION K_en(9,16) , K_ua(18,16), K_na(16,16), K_eu(9,18)
      DOUBLE PRECISION f_el(59)
      DOUBLE PRECISION f_u(18) , f_e(9), f_n(16), f_a(16)
      
      DOUBLE PRECISION STRESS(4), STRESS_h(10), ELEC_disp(2)
C
C     ==================================================================
C     Initialize
C     ==================================================================
      RHS = ZERO
      AMATRX = ZERO
      
      K_el = ZERO
      K_uu = ZERO
      K_ee = ZERO
      K_nn = ZERO
      K_en = ZERO
      K_ua = ZERO
      K_na = ZERO
      K_eu = ZERO
      
      f_el = ZERO
      f_u  = ZERO
      f_e  = ZERO
      f_n  = ZERO
      f_a  = ZERO
      
C
C     ==================================================================
C     Current State & Increment
C     ==================================================================
      STATE   = U(1:59)
      STATE_u = STATE(1:18)
      STATE_e = STATE(19:27)
      STATE_n = STATE(28:43)
      STATE_a = STATE(44:59)
      
      inc_u = DU(1:18 ,1)
      inc_e = DU(19:27,1)
      inc_n = DU(28:43,1)
      inc_a = DU(44:59,1)
      
      node9_xy = transpose(COORDS(1:2,1:9))
      node4_xy = transpose(COORDS(1:2,1:4))
            
C
C     ==================================================================
C     Material Property & Matrix
C     ==================================================================
      !E       = PROPS(1)
      !xnu     = PROPS(2)
      !sl      = PROPS(3)
      !f1      = PROPS(4)
      !f2      = PROPS(5)
      !permi   = PROPS(6)
      
      E = 100 * 1e-3
      xnu = 0.35
      sl = 2
      f1 = 30e-9
      f2 = 30e-9
      permi = 0.1*1e-9*1e-6
       
      xlam = E*xnu/((ONE+xnu)*(ONE-TWO*xnu))
      xmu = E/TWO/(ONE+xnu)
      D_uu(1,:) = [ xlam+TWO*xmu ,         xlam,         xlam, ZERO    ]
      D_uu(2,:) = [         xlam , xlam+TWO*xmu,         xlam, ZERO    ]
      D_uu(3,:) = [         xlam ,         xlam, xlam+TWO*xmu, ZERO    ]
      D_uu(4,:) = [         ZERO ,         ZERO,         ZERO,  xmu    ]
      
      D_nn = ZERO 
      D_nn(1:4,1:4) = D_uu
      D_nn(5:8,5:8) = D_uu
      D_nn(9,9)     = xmu
      D_nn(10,10)   = xmu
      D_nn = D_nn*sl*sl
      
      D_ee(1,:) = [  permi,  ZERO ]
      D_ee(2,:) = [   ZERO, permi ]
      
      D_en(1,:) = [ f1+TWO*f2,   f1,   f1,   ZERO,  ZERO,      ZERO,  ZERO,     f2,   f2, ZERO ]
      D_en(2,:) = [      ZERO, ZERO, ZERO,     f2,    f1, f1+TWO*f2,    f1,   ZERO, ZERO,   f2 ]

      D_eu(1,:) = [     0.,     0.,     0.,  3.26 ]
      D_eu(2,:) = [  1.644,  2.792,  1.644,    0. ]
      D_eu      = D_eu * 1e-9
      D_eu      = D_eu * (-ONE)
      
C
C     ==================================================================
C     Integration Points & Integration weights
C     ==================================================================
      CALL GAUSS_QUADRATURE_3(GP9)      
C
C     ==================================================================
C     Formulate Stiffness & internal force
C     ==================================================================
      DO i=1,size(GP9,1)
      
      wght = GP9(i,3)
      xi   = GP9(i,1)
      eta  = GP9(i,2)
      
      ! Quardratic shape function for displacement & electric potential
      CALL Q9_SHAPE_FUNCTIONS(xi, eta, N9, dN9dxi)      
      dxdxi = matmul(dN9dxi,node9_xy)                 !dxdxi = [ dx/dxi, dy/dxi ; dx/deta, dy/deta ]
      CALL MatInversion(dxdxi, dxidx, detJ9,2)        !dxidx = [ dxi/dx, deta/dx ; dxi/dy, deta/dy ]
      dN9dx = matmul(dxidx,dN9dxi)
      
      CALL Q4_SHAPE_FUNCTIONS(xi, eta, N4, dN4dxi)
      dxdxi = matmul(dN4dxi,node4_xy)                 !dxdxi = [ dx/dxi, dy/dxi ; dx/deta, dy/deta ]
      CALL MatInversion(dxdxi, dxidx, detJ4,2)        !dxidx = [ dxi/dx, deta/dx ; dxi/dy, deta/dy ]
      dN4dx = matmul(dxidx,dN4dxi)
      
      Axis_radius = dot_product(N9,node9_xy(:,1))
      
      N_u = ZERO
      N_u(1,1:9)   = N9
      N_u(2,10:18) = N9
      
      N_e = ZERO
      N_e(1,:) = N9
      
      N_n = ZERO
      N_n(1,1:4)    = N4
      N_n(2,5:8)    = N4
      N_n(3,9:12)   = N4
      N_n(4,13:16)  = N4
      
      N_a = ZERO
      N_a(1,1:4)    = N4
      N_a(2,5:8)    = N4
      N_a(3,9:12)   = N4
      N_a(4,13:16)  = N4

      B_e = ZERO
      B_e(1,:) = - dN9dx(1,:)
      B_e(2,:) = - dN9dx(2,:)
      
      B_u = ZERO
      B_u(1,1:9)   = dN9dx(1,:)
      B_u(2,10:18) = dN9dx(2,:)
      B_u(3,1:9)   = N9 / Axis_radius
      B_u(4,1:9)   = dN9dx(2,:)
      B_u(4,10:18) = dN9dx(1,:)

      B_n = ZERO
      B_n(1,1:4)    = dN4dx(1,:)
      B_n(2,5:8)    = dN4dx(1,:)
      B_n(3,9:12)   = dN4dx(1,:)
      B_n(4,13:16)  = dN4dx(1,:)
      B_n(5,1:4)    = dN4dx(2,:)
      B_n(6,5:8)    = dN4dx(2,:)
      B_n(7,9:12)   = dN4dx(2,:)  
      B_n(8,13:16)  = dN4dx(2,:)  
      B_n(9,1:4)    =   N4 / Axis_radius
      B_n(9,9:12)   = - N4 / Axis_radius
      B_n(10,13:16) =   N4 / Axis_radius
      
      K_uu = K_uu + matmul(matmul(transpose(B_u),D_uu),B_u) *detJ9*wght *Axis_radius*TWO*PI
      K_ee = K_ee - matmul(matmul(transpose(B_e),D_ee),B_e) *detJ9*wght *Axis_radius*TWO*PI
      K_nn = K_nn + matmul(matmul(transpose(B_n),D_nn),B_n) *detJ9*wght *Axis_radius*TWO*PI
      K_en = K_en - matmul(matmul(transpose(B_e),D_en),B_n) *detJ9*wght *Axis_radius*TWO*PI
      K_ua = K_ua - matmul(transpose(B_u),N_a)              *detJ9*wght *Axis_radius*TWO*PI
      K_na = K_na + matmul(transpose(N_n) ,N_a)             *detJ9*wght *Axis_radius*TWO*PI
      K_eu = K_eu - matmul(matmul(transpose(B_e),D_eu),B_u) *detJ9*wght *Axis_radius*TWO*PI

      
      STRESS    = matmul(D_uu,matmul(B_u,STATE_u)) - matmul(transpose(D_eu),matmul(B_e,STATE_e))
      STRESS_h  = matmul(D_nn,matmul(B_n,STATE_n)) - matmul(transpose(D_en),matmul(B_e,STATE_e))
      ELEC_disp = matmul(D_ee,matmul(B_e,STATE_e)) + matmul(D_en,matmul(B_n,STATE_n)) + matmul(D_eu,matmul(B_u,STATE_u))
      
      
      f_u = f_u + matmul(transpose(B_u),STRESS)              *detJ9*wght *Axis_radius*TWO*PI
     1          - matmul(transpose(B_u),matmul(N_a,STATE_a)) *detJ9*wght *Axis_radius*TWO*PI
      f_e = f_e - matmul(transpose(B_e),ELEC_disp)           *detJ9*wght *Axis_radius*TWO*PI
      f_n = f_n + matmul(transpose(B_n),STRESS_h)            *detJ9*wght *Axis_radius*TWO*PI
     1          + matmul(transpose(N_n),matmul(N_a,STATE_a)) *detJ9*wght *Axis_radius*TWO*PI
      f_a = f_a + matmul(transpose(N_a),matmul(N_n,STATE_n)) *detJ9*wght *Axis_radius*TWO*PI
     1          - matmul(transpose(N_a),matmul(B_u,STATE_u)) *detJ9*wght *Axis_radius*TWO*PI
      
      
      !Store Updated Variables at integration Point
      !globalSDV(JELEM-N_offset,i,1:3)   = STRESS
      !globalSDV(JELEM-N_offset,i,4:9)   = STRESS_h
      !globalSDV(JELEM-N_offset,i,10:11) = ELEC_disp
      !globalSDV(JELEM-N_offset,i,21:23) = matmul(B_u,STATE_u)
      !globalSDV(JELEM-N_offset,i,24:26) = matmul(N_n,STATE_n)
      
      END DO

      K_el = ZERO
      K_el( 1:18, 1:18) = K_uu
      K_el( 1:18,19:27) = transpose(K_eu)
      K_el( 1:18,44:59) = K_ua
      K_el(19:27, 1:18) = K_eu
      K_el(19:27,19:27) = K_ee
      K_el(19:27,28:43) = K_en
      K_el(28:43,19:27) = transpose(K_en)
      K_el(28:43,28:43) = K_nn
      K_el(28:43,44:59) = K_na
      K_el(44:59, 1:18) = transpose(K_ua)
      K_el(44:59,28:43) = transpose(K_na)
      
      f_el = ZERO
      f_el( 1:18) = f_u
      f_el(19:27) = f_e
      f_el(28:43) = f_n
      f_el(44:59) = f_a
      
      RHS(1:59,1) = - f_el
      AMATRX = K_el

      
      END


C
C ======================================================================
C     UMAT Subroutine for dummy element
C ======================================================================
       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)      
C
C     ==================================================================
C     Select  UMAT ELEMENT
C     ==================================================================

      CALL UMAT_2D_CPS(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)

       END      

C
C ======================================================================
C     UMAT 2D Plane Stress Element
C ======================================================================  
      SUBROUTINE UMAT_2D_CPS(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      USE GLOBAL     
C
      INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

C     ==================================================================
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 FIVE=5.D0,SIX=6.D0,SEVEN=7.D0,EIGHT=8.D0,NINE=9.D0,HALF=0.5D0)
      
      E=PROPS(1)
      xnu=PROPS(2)
      
      DDSDDE = ZERO
      STRESS = ZERO
      
      !Store variables in STATEV (STRESS, etc.)
      !STATEV(1:30) = globalSDV(NOEL,NPT,1:30)

      END

C
C     ==================================================================  
C     Print Array
C     ==================================================================      
      SUBROUTINE Print_Array(Description,A,n)
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*(*) Description
      DIMENSION A(n)

      write(6,*) Description
      
      write(6,999) (A(i) , i=1,n)

 999  format(100(E16.4,' '))

      RETURN
      END


C
C     ==================================================================  
C     Print Matrix  
C     ==================================================================      
      SUBROUTINE Print_Matrix(Description,A,m,n)
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*(*) Description
      DIMENSION A(m,n)

      write(6,*) Description
      
      DO i=1,m
          write(6,999) (A(i,j) , j=1,n)
      END DO

 999  format(100(E16.4,' '))

      RETURN
      END

      
      
C     =====================================================================
C      Q9 Element Shape Function & Gaussian Quadrature
C     =====================================================================
C
C       4            7            3
C       ��------------��------------��
C       ��                         �� 
C       ��    x7     x8      x9    ��
C       ��                         ��
C       ��                         ��
C      8��    x4      ��5      x6  ��6   
C       ��            9            ��
C       ��                         ��
C       ��    x1     x2      x3    ��
C       ��                         ��
C       ��------------��------------��
C       1            5            2
C
C     =====================================================================
      
      SUBROUTINE Q9_SHAPE_FUNCTIONS(xi, eta, N, dN)

      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 FIVE=5.D0,SIX=6.D0,SEVEN=7.D0,EIGHT=8.D0,NINE=9.D0,HALF=0.5D0)
      
      DOUBLE PRECISION  xi, eta
      DOUBLE PRECISION  N(9), dN(2, 9)
  
      ! Shape functions for a Q8 element
      N(1) =   ONE/FOUR * (ONE - xi) * xi * (ONE - eta) * eta
      N(2) = - ONE/FOUR * (ONE + xi) * xi * (ONE - eta) * eta
      N(3) =   ONE/FOUR * (ONE + xi) * xi * (ONE + eta) * eta
      N(4) = - ONE/FOUR * (ONE - xi) * xi * (ONE + eta) * eta
      N(5) = - HALF * (ONE + xi) * (ONE - xi)  * (ONE - eta) * eta
      N(6) =   HALF * (ONE + xi) * xi * (ONE - eta) * (ONE + eta)
      N(7) =   HALF * (ONE + xi) * (ONE - xi)  * (ONE + eta) * eta
      N(8) = - HALF * (ONE - xi) * xi * (ONE - eta) * (ONE + eta)
      N(9) =   (ONE + xi) * (ONE - xi) * (ONE - eta) * (ONE + eta)
  
      ! Derivatives of shape functions with respect to xi and eta
      dN(1, 1) =   ONE/FOUR * (ONE - TWO*xi) * (ONE - eta) * eta
      dN(1, 2) = - ONE/FOUR * (ONE + TWO*xi) * (ONE - eta) * eta
      dN(1, 3) =   ONE/FOUR * (ONE + TWO*xi) * (ONE + eta) * eta
      dN(1, 4) = - ONE/FOUR * (ONE - TWO*xi) * (ONE + eta) * eta
      dN(1, 5) = - HALF * (-TWO*xi) * (ONE - eta) * eta
      dN(1, 6) =   HALF * (ONE + TWO*xi) * (ONE - eta) * (ONE + eta)
      dN(1, 7) =   HALF * (-TWO*xi) * (ONE + eta) * eta
      dN(1, 8) = - HALF * (ONE - TWO*xi) * (ONE - eta) * (ONE + eta)
      dN(1, 9) =   (-TWO*xi) * (ONE - eta) * (ONE + eta)
  
      dN(2, 1) =   ONE/FOUR * (ONE - xi) * xi * (ONE - TWO*eta)
      dN(2, 2) = - ONE/FOUR * (ONE + xi) * xi * (ONE - TWO*eta)
      dN(2, 3) =   ONE/FOUR * (ONE + xi) * xi * (ONE + TWO*eta)
      dN(2, 4) = - ONE/FOUR * (ONE - xi) * xi * (ONE + TWO*eta)
      dN(2, 5) = - HALF * (ONE + xi) * (ONE - xi)  * (ONE - TWO*eta)
      dN(2, 6) =   HALF * (ONE + xi) * xi * (-TWO*eta)
      dN(2, 7) =   HALF * (ONE + xi) * (ONE - xi)  * (ONE + TWO*eta)
      dN(2, 8) = - HALF * (ONE - xi) * xi * (-TWO*eta)
      dN(2, 9) =   (ONE + xi) * (ONE - xi) * (-TWO*eta)
  
      END 
      
      
      SUBROUTINE Q9_GAUSS_QUADRATURE(Gauss_Point)
      
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 FIVE=5.D0,SIX=6.D0,SEVEN=7.D0,EIGHT=8.D0,NINE=9.D0,HALF=0.5D0)
      
      DOUBLE PRECISION  Gauss_Point(9,3)
      
      Gauss_Point(1,:) = [-sqrt(THREE/FIVE),-sqrt(THREE/FIVE), 25./81. ]
      Gauss_Point(2,:) = [            ZERO ,-sqrt(THREE/FIVE), 40./81. ]
      Gauss_Point(3,:) = [ sqrt(THREE/FIVE),-sqrt(THREE/FIVE), 25./81. ]
      Gauss_Point(4,:) = [-sqrt(THREE/FIVE),            ZERO , 40./81. ]
      Gauss_Point(5,:) = [            ZERO ,            ZERO , 64./81. ]
      Gauss_Point(6,:) = [ sqrt(THREE/FIVE),            ZERO , 40./81. ]
      Gauss_Point(7,:) = [-sqrt(THREE/FIVE), sqrt(THREE/FIVE), 25./81. ]
      Gauss_Point(8,:) = [            ZERO , sqrt(THREE/FIVE), 40./81. ]
      Gauss_Point(9,:) = [ sqrt(THREE/FIVE), sqrt(THREE/FIVE), 25./81. ]
      
      END 
      
      
      

C     =====================================================================
C      Q4 Element Shape Function & Gaussian Quadrature
C     =====================================================================
C
C       4                         3
C       ��-------------------------��
C       ��                         ��
C       ��    x3             x4    ��
C       ��                         ��
C       ��                         ��
C       ��                         ��
C       ��                         ��
C       ��                         ��
C       ��    x1             x2    ��
C       ��                         ��
C       ��-------------------------��
C       1                         2
C
C     =====================================================================
      
      SUBROUTINE Q4_SHAPE_FUNCTIONS(xi, eta, N, dN)

      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 FIVE=5.D0,SIX=6.D0,SEVEN=7.D0,EIGHT=8.D0,NINE=9.D0,HALF=0.5D0)
      
      DOUBLE PRECISION  xi, eta
      DOUBLE PRECISION  N(4), dN(2, 4)
  
      ! Shape functions for a Q8 element
      N(1) = ONE/FOUR * (ONE - xi) * (ONE - eta)
      N(2) = ONE/FOUR * (ONE + xi) * (ONE - eta)
      N(3) = ONE/FOUR * (ONE + xi) * (ONE + eta)
      N(4) = ONE/FOUR * (ONE - xi) * (ONE + eta)
  
      ! Derivatives of shape functions with respect to xi and eta
      dN(1, 1) = - ONE/FOUR * (ONE - eta)
      dN(1, 2) =   ONE/FOUR * (ONE - eta)
      dN(1, 3) =   ONE/FOUR * (ONE + eta)
      dN(1, 4) = - ONE/FOUR * (ONE + eta)
  
      dN(2, 1) = - ONE/FOUR * (ONE - xi)
      dN(2, 2) = - ONE/FOUR * (ONE + xi)
      dN(2, 3) =   ONE/FOUR * (ONE + xi)
      dN(2, 4) =   ONE/FOUR * (ONE - xi)
  
      END 
      
      SUBROUTINE Q4_GAUSS_QUADRATURE(Gauss_Point)
      
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 FIVE=5.D0,SIX=6.D0,SEVEN=7.D0,EIGHT=8.D0,NINE=9.D0,HALF=0.5D0)
      
      DOUBLE PRECISION  Gauss_Point(4,3)
      
      Gauss_Point(1,:) = [-sqrt(ONE/THREE), -sqrt(ONE/THREE), 1. ]
      Gauss_Point(2,:) = [ sqrt(ONE/THREE), -sqrt(ONE/THREE), 1. ]
      Gauss_Point(3,:) = [-sqrt(ONE/THREE),  sqrt(ONE/THREE), 1. ]
      Gauss_Point(4,:) = [ sqrt(ONE/THREE),  sqrt(ONE/THREE), 1. ]
      
      END 
      
      
      SUBROUTINE MatInversion(A,A_inv,det_A,size)
      INTEGER size
      DOUBLE PRECISION A(size,size),A_inv(size,size),det_A
C
      IF (size.EQ.2) THEN
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      A_inv(1,1) =  A(2,2) / det_A
      A_inv(1,2) = -A(1,2) / det_A
      A_inv(2,1) = -A(2,1) / det_A
      A_inv(2,2) =  A(1,1) / det_A
C
      ELSEIF (size.EQ.3) THEN
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     1        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     1        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      A_inv(1,1) = (A(2,2)*A(3,3)-A(3,2)*A(2,3)) / det_A
      A_inv(1,2) = (A(3,2)*A(1,3)-A(1,2)*A(3,3)) / det_A
      A_inv(1,3) = (A(1,2)*A(2,3)-A(2,2)*A(1,3)) / det_A
      A_inv(2,1) = (A(3,1)*A(2,3)-A(2,1)*A(3,3)) / det_A
      A_inv(2,2) = (A(1,1)*A(3,3)-A(3,1)*A(1,3)) / det_A
      A_inv(2,3) = (A(2,1)*A(1,3)-A(1,1)*A(2,3)) / det_A
      A_inv(3,1) = (A(2,1)*A(3,2)-A(3,1)*A(2,2)) / det_A
      A_inv(3,2) = (A(3,1)*A(1,2)-A(1,1)*A(3,2)) / det_A
      A_inv(3,3) = (A(1,1)*A(2,2)-A(2,1)*A(1,2)) / det_A
C
      ELSEIF (size.EQ.4) THEN
      det_A = 
     1    A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*
     1     A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
     1  - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*
     1     A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
     1  + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*
     1     A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
     1  - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*
     1     A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
      
      A_inv(1,1) = A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*
     1        A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
      A_inv(2,1) = A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*
     1        A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
      A_inv(3,1) = A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*
     1        A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
      A_inv(4,1) = A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*
     1        A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
      A_inv(1,2) = A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*
     1        A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
      A_inv(2,2) = A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*
     1        A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
      A_inv(3,2) = A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*
     1        A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
      A_inv(4,2) = A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*
     1        A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
      A_inv(1,3) = A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*
     1        A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
      A_inv(2,3) = A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*
     1        A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
      A_inv(3,3) = A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*
     1        A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
      A_inv(4,3) = A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*
     1        A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
      A_inv(1,4) = A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*
     1        A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
      A_inv(2,4) = A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*
     1        A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      A_inv(3,4) = A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*
     1        A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
      A_inv(4,4) = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*
     1        A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      A_inv = A_inv / det_A
      ENDIF

      RETURN
      END
C
C     =====================================================================
C      Q4 Element Shape Function & Gaussian Quadrature
C     =====================================================================
      
      SUBROUTINE GAUSS_QUADRATURE_2(GP)
      INTEGER I,J, Count
      DOUBLE PRECISION GP(4,3) , P(2,2)
      
      P(1,:) = [ -0.577350269189626 , 1.0 ]
      P(2,:) = [  0.577350269189626 , 1.0 ]
      
      GP = 0
      Count  = 0
      DO I=1,2
          DO J=1,2
              Count = Count + 1
              GP(Count,:) = [ P(J,1) , P(I,1), P(I,2)*P(J,2) ]
          END DO
      END DO
      
      END

      SUBROUTINE GAUSS_QUADRATURE_3(GP)
      INTEGER I,J, Count
      DOUBLE PRECISION GP(9,3) , P(3,2)
      
      P(1,:) = [ -0.774596669241483 , 0.555555555555556 ]
      P(2,:) = [                 0. , 0.888888888888889 ]
      P(3,:) = [  0.774596669241483 , 0.555555555555556 ]
      
      GP = 0
      Count  = 0
      DO I=1,3
          DO J=1,3
              Count = Count + 1
              GP(Count,:) = [ P(J,1) , P(I,1), P(I,2)*P(J,2) ]
          END DO
      END DO
      
      END
      
      SUBROUTINE GAUSS_QUADRATURE_4(GP)
      INTEGER I,J, Count
      DOUBLE PRECISION GP(16,3) , P(4,2)
      
      P(1,:) = [ -0.861136311594053 , 0.347854845137454 ]
      P(2,:) = [ -0.339981043584856 , 0.652145154862546 ]
      P(3,:) = [  0.339981043584856 , 0.652145154862546 ]
      P(4,:) = [  0.861136311594053 , 0.347854845137454 ]
      
      GP = 0
      Count  = 0
      DO I=1,4
          DO J=1,4
              Count = Count + 1
              GP(Count,:) = [ P(J,1) , P(I,1), P(I,2)*P(J,2) ]
          END DO
      END DO
      
      END
      
      
      SUBROUTINE GAUSS_QUADRATURE_5(GP)
      INTEGER I,J, Count
      DOUBLE PRECISION GP(25,3) , P(5,2)
      
      P(1,:) = [ -0.906179845938664 , 0.236926885056189 ]
      P(2,:) = [ -0.538469310105683 , 0.478628670499366 ]
      P(3,:) = [                0.0 , 0.568888888888889 ]
      P(4,:) = [  0.538469310105683 , 0.478628670499366 ]
      P(5,:) = [  0.906179845938664 , 0.236926885056189 ]

      GP = 0
      Count  = 0
      DO I=1,5
          DO J=1,5
              Count = Count + 1
              GP(Count,:) = [ P(J,1) , P(I,1), P(I,2)*P(J,2) ]
          END DO
      END DO
      
      END
      
      SUBROUTINE GAUSS_QUADRATURE_6(GP)
      INTEGER I,J, Count
      DOUBLE PRECISION GP(36,3) , P(6,2)
      
      P(1,:) = [ -0.932469514203152 , 0.171324492379170 ]
      P(2,:) = [ -0.661209386466265 , 0.360761573048139 ]
      P(3,:) = [ -0.238619186083197 , 0.467913934572691 ]
      P(4,:) = [  0.238619186083197 , 0.467913934572691 ]
      P(5,:) = [  0.661209386466265 , 0.360761573048139 ]
      P(6,:) = [  0.932469514203152 , 0.171324492379170 ]
      
      GP = 0
      Count  = 0
      DO I=1,6
          DO J=1,6
              Count = Count + 1
              GP(Count,:) = [ P(J,1) , P(I,1), P(I,2)*P(J,2) ]
          END DO
      END DO
      
      END
      
      SUBROUTINE GAUSS_QUADRATURE_7(GP)
      INTEGER I,J, Count
      DOUBLE PRECISION GP(49,3) , P(7,2)
      
      P(1,:) = [ -0.949107912342759 , 0.129484966168870 ]
      P(2,:) = [ -0.741531185599394 , 0.279705391489277 ]
      P(3,:) = [ -0.405845151377397 , 0.381830050505119 ]
      P(4,:) = [  0.000000000000000 , 0.417959183673469 ]
      P(5,:) = [  0.405845151377397 , 0.381830050505119 ]
      P(6,:) = [  0.741531185599394 , 0.279705391489277 ]
      P(7,:) = [  0.949107912342759 , 0.129484966168870 ]
      
      GP = 0
      Count  = 0
      DO I=1,7
          DO J=1,7
              Count = Count + 1
              GP(Count,:) = [ P(J,1) , P(I,1), P(I,2)*P(J,2) ]
          END DO
      END DO
      
      END
      
      SUBROUTINE GAUSS_QUADRATURE_8(GP)
      INTEGER I,J, Count
      DOUBLE PRECISION GP(64,3) , P(8,2)
      
      P(1,:) = [ -0.960289856497536 , 0.101228536290376 ]
      P(2,:) = [ -0.796666477413627 , 0.222381034453374 ]
      P(3,:) = [ -0.525532409916329 , 0.313706645877887 ]
      P(4,:) = [ -0.183434642495650 , 0.362683783378362 ]
      P(5,:) = [  0.183434642495650 , 0.362683783378362 ]
      P(6,:) = [  0.525532409916329 , 0.313706645877887 ]
      P(7,:) = [  0.796666477413627 , 0.222381034453374 ]
      P(8,:) = [  0.960289856497536 , 0.101228536290376 ]
                 
      GP = 0
      Count  = 0
      DO I=1,8
          DO J=1,8
              Count = Count + 1
              GP(Count,:) = [ P(J,1) , P(I,1), P(I,2)*P(J,2) ]
          END DO
      END DO
      
      END
      
      
      SUBROUTINE GAUSS_QUADRATURE_9(GP)
      INTEGER I,J, Count
      DOUBLE PRECISION GP(81,3) , P(9,2)
      
      P(1,:) = [ -0.968160239507626 , 0.081274388361574 ]
      P(2,:) = [ -0.836031107326636 , 0.180648160694857 ]
      P(3,:) = [ -0.613371432700590 , 0.260610696402935 ]
      P(4,:) = [ -0.324253423403809 , 0.312347077040003 ]
      P(5,:) = [                0.0 , 0.330239355001260 ]
      P(6,:) = [  0.324253423403809 , 0.312347077040003 ]
      P(7,:) = [  0.613371432700590 , 0.260610696402935 ]
      P(8,:) = [  0.836031107326636 , 0.180648160694857 ]
      P(9,:) = [  0.968160239507626 , 0.081274388361574 ]
      
      GP = 0
      Count  = 0
      DO I=1,9
          DO J=1,9
              Count = Count + 1
              GP(Count,:) = [ P(J,1) , P(I,1), P(I,2)*P(J,2) ]
          END DO
      END DO
      
      END
      
      SUBROUTINE GAUSS_QUADRATURE_10(GP)
      INTEGER I,J, Count
      DOUBLE PRECISION GP(100,3) , P(10,2)
      
      P(1,:)  = [ -0.973906528517172 , 0.066671344308688 ]
      P(2,:)  = [ -0.865063366688985 , 0.149451349150581 ]
      P(3,:)  = [ -0.679409568299024 , 0.219086362515982 ]
      P(4,:)  = [ -0.433395394129247 , 0.269266719309996 ]
      P(5,:)  = [ -0.148874338981631 , 0.295524224714753 ]
      P(6,:)  = [  0.148874338981631 , 0.295524224714753 ]
      P(7,:)  = [  0.433395394129247 , 0.269266719309996 ]
      P(8,:)  = [  0.679409568299024 , 0.219086362515982 ]
      P(9,:)  = [  0.865063366688985 , 0.149451349150581 ]
      P(10,:) = [  0.973906528517172 , 0.066671344308688 ]
      
      GP = 0
      Count  = 0
      DO I=1,10
          DO J=1,10
              Count = Count + 1
              GP(Count,:) = [ P(J,1) , P(I,1), P(I,2)*P(J,2) ]
          END DO
      END DO
      
      END
      
      