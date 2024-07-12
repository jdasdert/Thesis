module tensorLib
  implicit none
  real*8, parameter :: sqrt2 = dsqrt(2.d0), invsqrt2=0.5d0*dsqrt(2.d0)
  real*8, parameter :: PI = 2.d0*asin(1.d0)
  ! idexing and factors for Voigt notation
  integer, dimension(3,2), parameter :: idxVgt2D = reshape((/ 1,2,1,1,2,2 /), (/3,2/))
  integer, dimension(6,2), parameter :: idxVgt3D = reshape((/ 1,2,3,2,1,1,1,2,3,3,3,2 /), (/6,2/))
  real*8, dimension(3), parameter :: fctVgt2D=(/ 1.d0, 1.d0, sqrt2 /)
  real*8, dimension(6), parameter :: fctVgt3D=(/ 1.d0, 1.d0, 1.d0, sqrt2, sqrt2, sqrt2 /), & 
                                     invFctVgt3D = (/ 1.d0, 1.d0, 1.d0, invsqrt2, invsqrt2, invsqrt2 /)
  
  interface ddot
    module procedure ddot22, ddot23, ddot32, ddot24, ddot42, ddot44
  end interface

  interface dyadic
    module procedure dyadic11, dyadic21, dyadic22
  end interface

  interface getDistance
      module procedure getDistancePP, getDistancePPV, getDistancePPVV
  end interface

  interface one
    module procedure one3D, oneND
  end interface

  interface oneV
    module procedure one3DVgt
  end interface

  interface randZYZangle
    module procedure randZYZangle, randZYZangleMax
  end interface

  interface rotate
    module procedure rotate2, rotate4
  end interface

  interface rotateVgt
    module procedure rotate4Vgt
  end interface

  interface sbox
    module procedure ATsBoxA
  end interface

  interface sboxsI
    module procedure AsboxsI
  end interface

  interface writeMat
    module procedure writeMat, writeMatFormatted, writeMat2File
  end interface

  interface writeTens
    module procedure writeTens, writeTensFormatted
  end interface

  contains

!----------------------------------------------------------------------------------------------------------------------------------!
!--[SCALAR FUNCTIONS]--------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!

  function sg(a)
    implicit none
    real*8 :: a, sg
    if(a.gt.0.d0) then
      sg = 1.d0
    elseif(a.lt.0.d0) then
      sg = -1.d0
    else
      sg = 0.d0
    endif
  end function

!----------------------------------------------------------------------------------------------------------------------------------!
!--[IENTITIES]---------------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!

  ! returns identitiy for 3x3 matrix
  function one3D()
    implicit none
    real*8, dimension(3,3) :: one3D
    one3D = 0.d0
    one3D(1,1) = 1.d0; one3D(2,2) = 1.d0; one3D(3,3) = 1.d0
  end function

  ! returns identity for nxn matrix
  function oneND(N)
    implicit none
    integer :: N, i
    real*8, dimension(N,N) :: oneND
    oneND = 0.d0
    do i=1, N
      oneND(i,i) = 1.d0
    enddo
  end function

  ! returns identity for 3x3 matrix in voigt notation
  function one3DVgt()
    implicit none
    real*8, dimension(6) :: one3DVgt
    one3DVgt(1:3) = 1.d0
    one3DVgt(4:6) = 0.d0
  end function

!----------------------------------------------------------------------------------------------------------------------------------!
!--[ROTATION MATRICES]-------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!

  function rotMatZ(a) result(R)
    implicit none
    real*8 :: a, dmmy
    real*8, dimension(3,3) :: R
    dmmy = PI/180.d0*a
    R(1,:) = (/ cos(dmmy), sin(dmmy), 0.d0 /)
    R(2,:) = (/ -sin(dmmy), cos(dmmy), 0.d0 /)
    R(3,:) = (/ 0.d0, 0.d0, 1.d0 /)
  end function

  function rotMatY(a) result(R)
    implicit none
    real*8 :: a, dmmy
    real*8, dimension(3,3) :: R
    dmmy = PI/180.d0*a
    R(1,:) = (/ cos(dmmy), 0.d0, sin(dmmy) /)
    R(2,:) = (/ 0.d0, 1.d0, 0.d0 /)
    R(3,:) = (/ -sin(dmmy), 0.d0, cos(dmmy) /)
  end function

  function rotMatX(a) result(R)
    implicit none
    real*8 :: a, dmmy
    real*8, dimension(3,3) :: R
    dmmy = PI/180.d0*a
    R(1,:) = (/ 1.d0, 0.d0, 0.d0 /)
    R(2,:) = (/ 0.d0, cos(dmmy), sin(dmmy) /)
    R(3,:) = (/ 0.d0, -sin(dmmy), cos(dmmy) /)
  end function

  !------------------------------------------------------------------------------------------------------------------------------!

 
  function rotateZYZ(M,a1,a2,a3) result(Mrot)
    implicit none
    real*8 :: a1, a2, a3
    real*8, dimension(3,3) :: M
    real*8, dimension(size(M,1), size(M,1)) :: Mrot
    if(size(M,1).eq.size(M,2)) then
      Mrot = rotate(rotate(rotate(M,rotMatZ(a1)),rotMatY(a2)),rotMatZ(a3))
    else
      print*, 'ERROR STOP - rotateZYZ: Matrix dimensions not equal'
      ERROR STOP
    endif
  end function

  !------------------------------------------------------------------------------------------------------------------------------!

  function rotate2(M, R) result(Mrot)
    implicit none
    real*8, dimension(:,:) :: M
    real*8, dimension(size(M,1), size(M,2)) :: Mrot, R
    if (size(M,1).eq.size(M,2)) then
      Mrot = matmul(matmul(R, M), transpose(R))
    else
      print*, 'ERROR STOP - rotate22: Matrix not quadratic'
    endif
  end function

  function rotate4(T, R) result(Trot)
    implicit none
    real*8, dimension(:,:,:,:) :: T
    real*8, dimension(size(T,1), size(T,1), size(T,1), size(T,1)) :: Trot
    real*8, dimension(size(T,1), size(T,1)) :: R
    integer :: ndm, i,j,k,l
    ndm = size(T,1)
    if(ndm.eq.size(T,2).and.ndm.eq.size(T,3).and.ndm.eq.size(T,4)) then
      do i=1, ndm
        do j=1, ndm
          Trot(i,j,:,:) = rotate(T(i,j,:,:),R)
        enddo
      enddo
      do k=1, ndm
        do l=1, ndm
          !Trot(:,:,k,l) = rotate(Trot(:,:,k,l),transpose(R))
          Trot(:,:,k,l) = rotate(Trot(:,:,k,l),R)
        enddo
      enddo
    else
      print*, 'ERROR STOP - rotate44: Tensor with different dimensions'
      ERROR STOP
    endif
  end function

  function rotate4Vgt(V, R) result(Vrot)
    implicit none
    real*8, dimension(:,:) :: V, R
    real*8, dimension(size(V,1),size(V,1)):: Vrot, RTboxR
    RTboxR = sbox(R)
    Vrot = matmul(matmul(transpose(RTboxR), V), RTboxR)
    ! Vrot = T2V(rotate4(V2T(V),R))
  end function

!----------------------------------------------------------------------------------------------------------------------------------!
!--[VECTOR/MATRIX/TENSOR ROUTINES]-------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!
  
  ! calculates the angle between two vectors
  function getAngle(V1,V2) result(a)
    implicit none
    real*8, dimension(3) :: V1, V2
    real*8 :: a
    a = acos(dot_product(V1,V2)/(norm2(V1)*norm2(V2)))    
  end function

  function getDistancePP(P1,P2) result(d)
    implicit none
    real*8, dimension(3) :: P1,P2,d
    d = P2-P1
  end function 

  function getDistancePPV(P1,P2,V1) result(d)
    implicit none
    real*8, dimension(3) :: P1,P2,V1,d
    d = P2-P1-dot_product(P2-P1,V1/norm2(V1))*V1/norm2(V1)
  end function

  function getDistancePPVV(P1,P2,V1,V2) result(d)
    implicit none
    real*8, dimension(3) :: P1,P2,V1,V2,d
    d = cross(V1,V2)
    if(norm2(d).gt.0.d0) then
      d = d / norm2(d)
      d = dot_product(P2-P1,d) * d
    else
      d = getDistancePPV(P1,P2,V1)
    endif
  end function

  ! subroutine getIntersection(P1,P2,V1,V2,ip,a,b,dist,info)
  !   implicit none
  !   real*8, dimension(3) :: P1, P2, V1, V2, ip, distVec
  !   real*8 :: a,b,dist
  !   real*8, dimension(2,2) :: M
  !   integer :: info, ipiv(2)
  !   logical :: parallel

  !   info = 0

  !   ! get the distance between the vector V1 going through point P1 and the vector V2 going through point P2
  !     distVec = getDistance(P1,P2,V1,V2)
  !     dist = norm2(distVec)
  !   ! --

  !   ! set up matrix for equation system to find intersection P1+a*V1 = P2+b*V2 => (V1 -V2) (a b)^T = P2-P1
  !     M(:,1) = V1(1:2)
  !     M(:,2) = -V2(1:2)
  !     if(det(M).eq.0.d0) then
  !       M(:,1) = V1(2:3)
  !       M(:,2) = -V2(2:3)
  !       if(det(M).eq.0.d0) then
  !         M(:,1) = (/ V1(1), V1(3) /)
  !         M(:,2) = (/ V2(1), V2(3) /)
  !         if(det(M).eq.0.d0) then
  !           info = 11
  !         endif
  !       endif
  !     endif
  !   ! --

  !   ! if V1 and V2 are parallel
  !     if(info.eq.11) then
  !       ip = 0.5d0*(P2+P1)
  !       a = dot_product(ip-P1,V1/norm2(V1))
  !       b = dot_product(ip-P2,V2/norm2(V2))
  !   ! --

  !   ! if V1 and V2 are not parallel
  !     else
  !       ip = P2-dist-P1
  !       call dgesv(2,1,M,2,ipiv,ip(1:2),2,info)
  !       a = ip(1)
  !       b = ip(2)
  !       ip = P1 + a*V1 + 0.5d0*distVec
  !     endif
  !   ! --
  ! end subroutine


  function cross(V,W) result (c)
    implicit none
    real*8, dimension(3) :: V,W,c
    c(1) = V(2)*W(3) - V(3)*W(2)
    c(2) = V(3)*W(1) - V(1)*W(3)
    c(3) = V(1)*W(2) - V(2)*W(1)
  end function

  ! calculates double contraction (numbers give the order of tensor)
  function ddot22(a,b) result (dd)
    implicit none
    real*8, dimension(:,:) :: a,b
    real*8 :: dd
    integer :: i,j
    if(size(a,1).eq.size(b,1).and.size(a,2).eq.size(b,2)) then
      dd = 0.d0
      do i=1, size(a,1)
        do j=1, size(a,2)
          dd = dd + a(i,j)*b(i,j)
        enddo
      enddo
    else
      print*, 'ERROR STOP - ddot22: Matrix dimensions not equal'
      ERROR STOP
    endif
  end function

  function ddot32(a,b) result(dd)
    implicit none
    real*8, dimension(:,:,:) :: a
    real*8, dimension(size(a,2), size(a,3)) :: b
    real*8, dimension(size(a,1)) :: dd
    integer :: i,j,k
    dd = 0.d0
    do i=1, size(a,1)
      do j=1, size(a,2)
        do k=1, size(a,3)
          dd(i) = dd(i) + a(i,j,k)*b(j,k)
        enddo
      enddo
    enddo
  end function

  function ddot23(a,b) result(dd)
    implicit none
    real*8, dimension(:,:,:) :: b
    real*8, dimension(size(b,1),size(b,2)) :: a
    real*8, dimension(size(b,3)) :: dd
    integer :: i,j,k
    dd = 0.d0
    do i=1, size(b,1)
      do j=1, size(b,2)
        do k=1, size(b,3)
          dd(k) = dd(k) + a(i,j)*b(i,j,k)
        enddo
      enddo
    enddo
  end function

  function ddot24(a,b) result (dd)
    implicit none
    real*8, dimension(:,:) :: a
    real*8, dimension(:,:,:,:) :: b
    real*8, dimension(size(b,3), size(b,4)) :: dd
    integer :: i,j,k,l
    if(size(a,1).eq.size(b,1).and.size(a,2).eq.size(b,2)) then
      dd = 0.d0
      do i=1, size(b,1)
        do j=1, size(b,2)
          do k=1, size(b,3)
            do l=1, size(b,4)
              dd(k,l) = dd(k,l) + a(i,j)*b(i,j,k,l)
            enddo
          enddo
        enddo
      enddo
    else
      print*, 'ERROR STOP - ddot24: Matrix and tensor dimensions do not match'
      ERROR STOP
    endif
  end function

  function ddot42(a,b) result (dd)
    implicit none
    real*8, dimension(:,:,:,:) :: a
    real*8, dimension(:,:) :: b
    real*8, dimension(size(a,1), size(a,2)) :: dd
    integer :: i,j,k,l
    ! if(size(a,1).eq.size(b,1).and.size(a,2).eq.size(b,2)) then
      dd = 0.d0
      do i=1, size(a,1)
        do j=1, size(a,2)
          do k=1, size(a,3)
            do l=1, size(a,4)
              dd(i,j) = dd(i,j) + a(i,j,k,l)*b(k,l)
            enddo
          enddo
        enddo
      enddo
    ! else
    !   print*, 'ERROR STOP - ddot42: Matrix and tensor dimensions do not match'
    !   ERROR STOP
    ! endif
  end function

  function ddot44(a,b) result (dd)
    implicit none
    real*8, dimension(:,:,:,:) :: a,b
    real*8, dimension(size(a,1), size(a,2), size(a,3), size(a,4)) :: dd
    integer :: i,j,k,l,m,n, ndm
    ndm = size(a,1)
    ! if(size(a,1).eq.size(b,1).and.size(a,2).eq.size(b,2)) then
      dd = 0.d0
      do i=1, ndm
        do j=1, ndm
          do k=1, ndm
            do l=1, ndm
              do m=1, ndm
                do n=1, ndm
                  dd(i,j,m,n) = dd(i,j,m,n) + a(i,j,k,l)*b(k,l,m,n)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    ! else
    !   print*, 'ERROR STOP - ddot42: Matrix and tensor dimensions do not match'
    !   ERROR STOP
    ! endif
  end function

  !------------------------------------------------------------------------------------------------------------------------------!

  ! calculates determinant of matrix M
  function det(M)
    implicit none
    real*8 :: det
    real*8, dimension(:,:) :: M
    integer :: ndm
    if (size(M,1).eq.size(M,2)) then
      ndm = size(M,1)
      if(ndm.eq.2) then
        det = M(1,1)*M(2,2) - M(1,2)*M(2,1)
      elseif(ndm.eq.3) then
        det = M(1,1)*M(2,2)*M(3,3) - M(1,1)*M(2,3)*M(3,2) &
            + M(1,2)*M(2,3)*M(3,1) - M(1,2)*M(2,1)*M(3,3) & 
            + M(1,3)*M(2,1)*M(3,2) - M(1,3)*M(2,2)*M(3,1)
      else
        det = detLU(ndm, M)
        ! print*, 'ERROR STOP - det: not implemented for ndm=', ndm
      endif
    else
      print*, 'ERROR STOP - det: Matrix not quadratic'
      ERROR STOP
    endif
  end function

  function detLU(N, mat) 
    implicit none 
    real*8, dimension(:,:) :: mat
    real*8 :: detLU, sgn
    integer :: N
    integer, allocatable :: ipiv(:)
    integer :: i, info
    
    N = N
    allocate(ipiv(N))
    ipiv = 0
    
    call dgetrf(N, N, mat, N, ipiv, info)
    
    detLU = 1
    
    do i=1, N
       detLU = detLU*mat(i,i)
    end do  
    
    sgn = 1
    
    do i=1, N
    
       if(ipiv(i) /= i) then
        sgn=-sgn
       else
        sgn=sgn
       endif
    end do 
    
    detLU=sgn*detLU
      
  end function

  !------------------------------------------------------------------------------------------------------------------------------!

  ! calculates deviatoric part of matrix M
  function dev(M)
    implicit none
    real*8, dimension(:,:) :: M
    real*8, dimension(size(M,1),size(M,1)) :: dev
    real*8 :: tr
    integer :: ndm
    ndm = size(M,1)
    tr = trace(M)
    dev = M - 1.d0/3.d0*tr*one(ndm)
  end function

  !------------------------------------------------------------------------------------------------------------------------------!

  ! calculates dyadic product of two vectors
  function dyadic11(a,b) result (d)
    implicit none
    real*8, dimension(:) :: a,b
    real*8, dimension(size(a),size(b)) :: d
    integer :: i,j
    do i=1, size(a)
      do j=1, size(b)
        d(i,j) = a(i)*b(j)
      enddo
    enddo
  end function

  ! calculates dyadic product of matrix and vector
  function dyadic21(a,b) result (d)
    implicit none
    real*8, dimension(:,:) :: a
    real*8, dimension(:) :: b
    real*8, dimension(size(a,1),size(a,2),size(b)) :: d
    integer :: i,j,k
    do i=1, size(a,1)
      do j=1, size(a,2)
        do k=1, size(b)
          d(i,j,k) = a(i,j)*b(k)
        enddo
      enddo
    enddo
  end function

  ! calculates dyadic product of two matrices
  function dyadic22(a,b) result(d)
    implicit none
    real*8, dimension(:,:) :: a,b
    real*8, dimension(size(a,1),size(a,2),size(b,1),size(b,2)) :: d
    integer :: i,j,k,l
    do i=1, size(a,1)
      do j=1, size(a,2)
        do k=1, size(b,1)
          do l=1, size(b,2)
            d(i,j,k,l) = a(i,j)*b(k,l)
          enddo
        enddo
      enddo
    enddo
  end function

  !------------------------------------------------------------------------------------------------------------------------------!

  ! calculates the inverse of matrix a
  function inv(a)
    implicit none
    real*8, dimension(:,:) :: a
    real*8, dimension(size(a,1), size(a,1)) :: inv
    integer, dimension(size(a,1)) :: ipiv
    real*8, dimension(size(a,1)) :: work
    integer :: ndm, info
    inv = a
    if(size(a,1).eq.size(a,2)) then
      ndm = size(a,1)
      call dgetrf(ndm, ndm, inv, ndm, ipiv, info)
      if(info.ne.0) then
        print*, 'ERROR STOP - inv: dgetrf with problems'
        ERROR STOP
      else
        call dgetri(ndm, inv, ndm, ipiv, work, ndm, info)
        if(info.ne.0) then
          print*, 'ERROR STOP - inv: dgetri with problems'
          ERROR STOP
        endif
      endif
    else
      print*, 'ERROR STOP - inv: Matrix dimensions not equal'
      ERROR STOP
    endif
  end function

  !------------------------------------------------------------------------------------------------------------------------------!

  ! calculates the symmetric part of matrix M
  function sym(M)
    implicit none
    real*8, dimension(:,:) :: M
    real*8, dimension(size(M,1), size(M,2)) :: sym
    if(size(M,1).eq.size(M,2)) then
      sym = 0.5d0*(M+transpose(M))
    else
      print*, 'ERROR STOP - sym: Matrix dimensions not equal'
      ERROR STOP
    endif
  end function

  !------------------------------------------------------------------------------------------------------------------------------!

  ! calculates the trace of matrix a
  function trace(a)
    real*8 :: trace
    real*8, dimension(:,:) :: a
    integer :: i
    if(size(a,1).eq.size(a,2)) then
      trace = 0.d0
      do i=1, size(a,1)
        trace = trace + a(i,i)
      enddo 
    else
      print*, 'ERROR STOP - trace: Matrix not quadratic'
      ERROR STOP
    endif
  end function

!----------------------------------------------------------------------------------------------------------------------------------!
!--[VOIGT NOTATION ROUTINES]-------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!

  subroutine getEffectiveStress(sig, sigEff, ndm)
    implicit none
    integer :: ndm
    real*8 :: sigEff
    real*8, dimension(ndm, ndm) :: sig, devSig
    devSig = dev(sig)
    sigEff = sqrt(ddot(devSig,devSig))
  end subroutine

  subroutine getEffectiveStressVgt(sig, sigEff, ndm, nTens)
    implicit none
    integer :: ndm, nTens
    real*8 :: sigEff
    real*8, dimension(nTens) :: sig
    real*8, dimension(ndm, ndm) :: devSig
    devSig = dev(V2M(sig))
    sigEff = sqrt(ddot(devSig,devSig))
  end subroutine

  subroutine getVonMisesStress(sig, sigVm, ndm)
    implicit none
    integer :: ndm
    real*8 :: sigVM
    real*8, dimension(ndm, ndm) :: sig, devSig
    devSig = dev(sig)
    sigVM = sqrt(1.5d0*ddot(devSig,devSig))
  end subroutine

  subroutine getVonMisesStressVgt(sig, sigVm, ndm, nTens)
    implicit none
    integer :: ndm, nTens
    real*8 :: sigVM
    real*8, dimension(nTens) :: sig
    real*8, dimension(ndm, ndm) :: devSig
    devSig = dev(V2M(sig))
    sigVM = sqrt(1.5d0*ddot(devSig,devSig))
  end subroutine

!----------------------------------------------------------------------------------------------------------------------------------!
!--[VOIGT NOTATION ROUTINES]-------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!

  function M2V(M)
    implicit none
    real*8, dimension(:,:) :: M
    real*8, dimension(size(M,1)*3-3) :: M2V
    real*8, parameter :: tol=1.d-14
    integer :: ndm
    if(size(M,1).eq.size(M,2)) then
      if(norm2(M-transpose(M)).gt.tol) then
        ! print*, 'M2V: Matrix not symmetric?'
        ! call writeMat(M)
        ! M=sym(M)
      endif
      ndm = size(M,1)
      if(ndm.eq.2) then
        M2V(1) = M(1,1)
        M2V(2) = M(2,2)
        M2V(3) = sqrt2*M(2,1)
      elseif(ndm.eq.3) then
        M2V(1) = M(1,1)
        M2V(2) = M(2,2)
        M2V(3) = M(3,3)
        M2V(4) = sqrt2*M(3,2)
        M2V(5) = sqrt2*M(3,1)
        M2V(6) = sqrt2*M(2,1)
      else
        print*, 'ERROR STOP - M2V: not implemented for ndm=', ndm
        ERROR STOP
      endif
    else
      print*, 'ERROR STOP - M2V: Matrix dimensions not equal'
      ERROR STOP
    endif
  end function

  !------------------------------------------------------------------------------------------------------------------------------!

  function V2M(V)
    implicit none
    real*8, dimension(:) :: V
    real*8, dimension(size(V)/3+1,size(V)/3+1) :: V2M
    real*8, parameter :: tol=1.d-15
    integer :: nTens
    nTens = size(V)
    if(nTens.eq.3) then
      V2M(1,:) = (/ V(1), invsqrt2*V(3) /)
      V2M(2,:) = (/ invsqrt2*V(3), V(2) /)
    elseif(nTens.eq.6) then
      V2M(1,:) = (/ V(1), invsqrt2*V(6), invsqrt2*V(5) /)
      V2M(2,:) = (/ invsqrt2*V(6), V(2), invsqrt2*V(4) /)
      V2M(3,:) = (/ invsqrt2*V(5), invsqrt2*V(4), V(3) /)
    else
      print*, 'ERROR STOP - V2M: not implemented for nTens=', nTens
      ERROR STOP
    endif

  end function

  !------------------------------------------------------------------------------------------------------------------------------!

  function T2V(T)
    implicit none
    real*8, dimension(:,:,:,:) :: T
    real*8, dimension(size(T,1)*3-3, size(T,1)*3-3) :: T2V
    integer, dimension(6,2),parameter :: p = idxVgt3D
    real*8, dimension(6), parameter :: f = fctVgt3D
    integer :: ndm, i,j
    if(size(T,1).eq.size(T,2).and.size(T,1).eq.size(T,3).and.size(T,1).eq.size(T,4)) then
      ndm = size(T,1)
      if(ndm.eq.2) then
        T2V(1,:) = (/ T(1,1,1,1), T(1,1,2,2), sqrt2*T(1,1,2,1) /)
        T2V(2,:) = (/ T(2,2,1,1), T(2,2,2,2), sqrt2*T(2,2,2,1) /)
        T2V(3,:) = (/ sqrt2*T(2,1,1,1), sqrt2*T(2,1,2,2), 2.d0*T(2,1,2,1) /)
      elseif(ndm.eq.3) then
        do i=1, ndm*3-3
          do j=1, ndm*3-3
            T2V(i,j) = f(i)*f(j) * T(p(i,1),p(i,2),p(j,1),p(j,2))
          enddo
        enddo
      else
        print*, 'ERROR STOP - T2V: not implemented for ndm=', ndm
        ERROR STOP
      endif
    else
      print*, 'ERROR STOP - T2V: Tensor dimensions are not equal'
      ERROR STOP
    endif
  end function

  !------------------------------------------------------------------------------------------------------------------------------!

  function V2T(V)
    implicit none
    real*8, dimension(:,:) :: V
    real*8, dimension(size(V,1)/3+1, size(V,1)/3+1, size(V,1)/3+1, size(V,1)/3+1) :: V2T
    integer, dimension(6,2), parameter :: p = idxVgt3D
    real*8, dimension(6), parameter :: f = invFctVgt3D
    integer :: i,j,nTens
    nTens = size(V,1)
    if(nTens.eq.3) then
      V2T(1,1,1,:) = (/ V(1,1), V(1,3) /)
      V2T(1,1,2,:) = (/ V(1,3), V(1,2) /)
      V2T(2,2,1,:) = (/ V(2,1), V(2,3) /)
      V2T(2,2,2,:) = (/ V(2,3), V(2,2) /)
      V2T(1,2,1,:) = (/ V(3,1), V(3,3)/)
      V2T(1,2,1,:) = (/ V(3,1), V(3,3)/)
      V2T(2,1,:,:) = V2T(1,2,:,:)
    elseif(nTens.eq.6) then
      do i=1, nTens
        do j=1, nTens
          V2T(p(i,1),p(i,2),p(j,1),p(j,2)) = f(i)*f(j)*V(i,j)
          V2T(p(i,2),p(i,1),p(j,1),p(j,2)) = f(i)*f(j)*V(i,j)
          V2T(p(i,1),p(i,2),p(j,2),p(j,1)) = f(i)*f(j)*V(i,j)
          V2T(p(i,2),p(i,1),p(j,2),p(j,1)) = f(i)*f(j)*V(i,j)
        enddo
      enddo
    else
      print*, 'ERROR STOP - V2T: not implemented for nTens=', nTens
      ERROR STOP
    endif
  end function

  !------------------------------------------------------------------------------------------------------------------------------!

  function matmulVV2M(V,W) result(M)
    implicit none
    real*8, dimension(:) :: V, W
    real*8, dimension(size(V,1)/3+1, size(V,1)/3+1) :: M
    M = matmul(V2M(V), V2M(W))
  end function

!----------------------------------------------------------------------------------------------------------------------------------!
!--[BOX ROUTINES]-------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!

  function ATsBoxA(A) result (sbox)
    real*8, dimension(:,:) :: A
    real*8, dimension(size(A,1)*3-3, size(A,1)*3-3) :: sbox
    integer, dimension(6,2) :: p = idxVgt3D
    real*8, dimension(6) :: f = fctVgt3D
    integer :: i, j, ndm, nTens
    ndm = size(A,1); nTens=ndm*3-3
    if(ndm.eq.2) then
      sbox(1,:) = (/ A(1,1)*A(1,1), A(2,1)*A(2,1), sqrt2*A(1,1)*A(2,1)/)
      sbox(1,:) = (/ A(1,2)*A(1,2), A(2,2)*A(2,2), sqrt2*A(1,2)*A(2,2)/)
      sbox(1,:) = (/ sqrt2*A(1,1)*A(1,2), sqrt2*A(2,1)*A(2,2), A(1,1)*A(2,2)+A(1,2)*A(2,1)/)
    elseif(ndm.eq.3) then
      do i=1, nTens
        do j=1, nTens
          sbox(i,j) = f(i)*f(j)*0.5d0*( A(p(j,1),p(i,1))*A(p(j,2),p(i,2)) + A(p(j,1),p(i,2))*A(p(j,2),p(i,1)))
        enddo
      enddo
    else
      print*, 'ERROR STOP - ATsBoxA: not implemented for ndm=', ndm
      ERROR STOP
    endif    
  end function
  
  function AsboxsI(aa)        
    implicit none
    integer :: ndm
    real*8 :: aa(:,:), AsboxsI(3*(size(aa,1)-1),3*(size(aa,1)-1))
    real*8, parameter :: w2_2 = sqrt(2.d0)/2.d0
    ndm = size(aa,1)
    if (ndm == 2) then
        AsboxsI(1,:) = (/ aa(1,1), 0.d0, aa(1,2)*w2_2 /)
        AsboxsI(2,:) = (/ 0.d0, aa(2,2), aa(1,2)*w2_2 /)
        AsboxsI(3,:) = (/ aa(1,2) * w2_2, aa(1,2) * w2_2, 0.5d0*( aa(1,1)+aa(2,2) ) /)
    elseif (ndm == 3) then
        AsboxsI(1,:) = (/ aa(1,1), 0.d0, 0.d0, 0.d0, w2_2*aa(1,3), w2_2*aa(1,2) /)
        AsboxsI(2,:) = (/ 0.d0, aa(2,2), 0.d0, w2_2*aa(2,3), 0.d0, w2_2*aa(1,2) /)
        AsboxsI(3,:) = (/ 0.d0, 0.d0, aa(3,3), w2_2*aa(2,3), w2_2*aa(1,3), 0.d0 /)
        AsboxsI(4,:) = (/ 0.d0, w2_2*aa(2,3), w2_2*aa(2,3), 0.5d0*(aa(2,2)+aa(3,3)), 0.5d0*aa(1,2), 0.5d0*aa(1,3) /)
        AsboxsI(5,:) = (/ w2_2*aa(1,3), 0.d0, w2_2*aa(1,3), 0.5d0*aa(1,2), 0.5d0*(aa(1,1)+aa(3,3)), 0.5d0*aa(2,3) /)
        AsboxsI(6,:) = (/ w2_2*aa(1,2), w2_2*aa(1,2), 0.d0, 0.5d0*aa(1,3), 0.5d0*aa(2,3), 0.5d0*(aa(1,1)+aa(2,2)) /)
    else
      error stop
    endif
  end function

!----------------------------------------------------------------------------------------------------------------------------------!
!--[RANDOM ROUTINES]---------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!

  subroutine randZYZangle(zyz)
    implicit none
    real*8 :: dir(3), zyz(3)
    real*8, parameter :: PI = 2.d0 * asin(1.d0)
  
    dir = randDir3D()
    call random_number(zyz(1))
    zyz(1) = zyz(1) * 360.d0
    zyz(2) = (90.d0 - dir(2)/abs(dir(2)) * acos(norm2(dir(1:2))/norm2(dir)) / PI * 180.d0)
    zyz(3) = acos(dir(1)/norm2(dir(1:2))) /PI * 180.d0
  end subroutine

  subroutine randZYZangleMax(zyz,zyzmax)
    implicit none
    real*8 :: dir(3), zyz(3), zyzmax(3)
    real*8, parameter :: PI = 2.d0 * asin(1.d0)

    dir = randDir3D()
    call random_number(zyz(1))
    zyz(1) = zyz(1) * zyzmax(1)
    zyz(2) = (90.d0 - dir(2)/abs(dir(2)) * acos(norm2(dir(1:2))/norm2(dir)) / PI * 180.d0)
    zyz(3) = acos(dir(1)/norm2(dir(1:2))) /PI * zyzmax(3)
  end subroutine

  function randDir3D() result(dir)
    implicit none
    real*8 :: dir(3)
    do 
      call random_stdnormal(dir(1))
      call random_stdnormal(dir(2))
      call random_stdnormal(dir(3))
      if(norm2(dir).gt.1.d-10) exit
    enddo
    dir = dir/norm2(dir)
  end function

  subroutine randn(x,mean,sd)
    implicit none
    real*8 :: x, mean, sd
    integer :: i
    call random_stdnormal(x)
    x = sd * x + mean
  end subroutine
  
  subroutine random_stdnormal(x)
    implicit none
    real*8,intent(out) :: x
    real*8,parameter :: pi=3.14159265
    real*8 :: u1,u2
    call random_stduniform(u1)
    call random_stduniform(u2)
    x = dsqrt(-2*log(u1))*cos(2*pi*u2)
  end subroutine random_stdnormal
  
  subroutine random_stduniform(u)
    implicit none
    real*8,intent(out) :: u
    real*8 :: r
    call random_number(r)
    u = 1 - r
  end subroutine random_stduniform

!----------------------------------------------------------------------------------------------------------------------------------!
!--[OUTPUT ROUTINES]---------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!

  subroutine writeMat(M)
    implicit none
    real*8 :: M(:,:)
    integer :: i
    do i=1, size(M,1)
      print*, M(i,:)
    enddo
  end subroutine

  subroutine writeMatFormatted(M, fstr)
    implicit none
    real*8 :: M(:,:)
    character(*) :: fstr
    character*128 :: fmtstr
    integer :: i
    write(fmtstr,*) '(', size(M,2), trim(fstr), ')'
    do i=1, size(M,1)
      print trim(fmtstr), M(i,:)
    enddo
  end subroutine


  !------------------------------------------------------------------------------------------------------------------------------!

  subroutine writeMat2File(M,u,fstr)
    implicit none
    real*8 :: M(:,:)
    integer :: u,i
    character(*) :: fstr
    do i=1, size(M,1)
      write(u,fstr) M(i,:)
    enddo
  end subroutine

  !------------------------------------------------------------------------------------------------------------------------------!

  subroutine writeTensFormatted(T, fstr)
    implicit none
    character(*) :: fstr
    character*128 :: str, fmtstr
    real*8 :: T(:,:,:,:)
    integer :: i,j,k
    write(fmtstr,*) '(A,', size(T,1), trim(fstr), ',A)'
    do i=1, size(T,1)
      do k=1, size(T,3)
        str = ''
        do j=1, size(T,2)
          write(str,trim(fmtstr)) trim(str), T(i,j,k,:), ' | '
        enddo
        print*, str
      enddo
      print*, ''
    enddo
  end subroutine

  subroutine writeTens(T)
    implicit none
    character*256 :: str
    real*8 :: T(:,:,:,:)
    integer :: i,j,k
    do i=1, size(T,1)
      do k=1, size(T,3)
        str = ''
        do j=1, size(T,2)
          write(str,*) trim(str), T(i,j,k,:), '|'
        enddo
        print*, str
      enddo
      print*, ''
    enddo
  end subroutine

end module