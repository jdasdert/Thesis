program createMesh
  use tensorLib
  use sortLib
  implicit none

  !________________________________________________________________________________________________________________________________!
  !--[VARIABLES]--------------------------------------------------------------------------------------------------------------------

    real*8, parameter :: third = 1.d0/3.d0

    real*8, parameter :: alMax = 20.d0, alMin = 10.d0, arMax=1.d0, arMin=0.25d0, rhoZnO=5.61d-12
    real*8 :: Vtpod, tpoddim(2), boxSize(3), Vsample, rhoSample, rxIn(2), rxOut(2), centerIn(2), centerOut(2), centerShift(2)
    real*8 :: VtpodFact, ryIn(2), ryOut(2)
    integer :: nParticles, shape, nBoxesX,nBoxesY,nBoxesZ, nBoxes
    real*8, allocatable :: pos(:,:), nodes(:,:), al(:), ar(:), dd(:,:,:), ang(:,:), er(:), erRed(:)
    integer,allocatable :: nodeCM(:), elmtCM(:,:), elmtRedCM(:,:), poreElmtCM(:,:)
    real*8, allocatable :: nodesRed(:,:)
    integer :: nNodes, nNodesUnconnected, nPores, nPoresRed, nNodesFinal, nElmtFinal
    character :: filename*128, str1*32, arg*32, str*128, dl*128, lf=char(10)
    integer :: i,j,k, info, nargs, nInter, nElmt, nElmtRed, nNodesRed, outpt
    integer, allocatable :: boxBoxCM(:,:),partBoxCM(:),boxPartCM(:,:),nPartInBox(:),nodeBoxCM(:),boxNodeCM(:,:),nNodeInBox(:)
    integer, allocatable :: boxElmtCM(:,:),nElmtInBox(:), nElmtInPore(:), nodePoreCM(:)
    real*8 :: sd_volume, min_volume, lr_ratio, sd
    dl(:) = ' '
    ! --
  
  !________________________________________________________________________________________________________________________________!
  !--[CODE]-------------------------------------------------------------------------------------------------------------------------

  ! define default values for sample geometry, number of particles and tetrapod dimensions
      shape = 1
      rxIn = (/ 10.d0, 10.d0 /) ! inner radius (top,bottom) [µm]
      ryIn = (/ 10.d0, 10.d0 /)
      rxOut = (/ 100.d0, 100.d0 /) ! outer radius (top,bottom) [µm]
      ryOut = (/ 100.d0, 100.d0 /)
      boxSize(3) = 50.d0
      tpoddim = (/ 20.d0, 1.333d0 /)
      VtpodFact = 1.d0
      rhoSample = 0.25d-12
      centerShift = (/ 0.d0, 0.d0 /)
      sd = 0.0d0
      outpt = 1
  ! --

  call readInputArgs()
 
  filename = 'mesh'

  print*, ''
  print*, '--[ SAMPLE ]------------------------------------------------------------------------------------------------------------'

  call createPosByRejection()
  call createTpods()
  call writeTetrapods2VTK()
  nBoxesX = int(boxSize(1)/(2*alMax))+2
  nBoxesY = int(boxSize(2)/(2*alMax))+2
  nBoxesZ = int(boxSize(3)/(2*alMax))+2
  nBoxes = nBoxesX*nBoxesY*nBoxesZ
  call createBoxBoxCM()
  call createBoxPartCM()

  print*, ''
  print*, '--[ MESH ]--------------------------------------------------------------------------------------------------------------'

  call findIntersectionsAndCreateNodesAndElements()
  call createBoxNodeCM(nodes(1:nNodes,1:3),nNodes)
  call connectUnconnectedNodes()
  call writeNodesAndElements2VTK()

  print*, ''
  print*, '--[ REDUCED MESH ]------------------------------------------------------------------------------------------------------'

  call combineNodeClusters()
  deallocate(boxNodeCM)
  deallocate(nNodeInBox)
  deallocate(nodeBoxCM)
  call createBoxNodeCM(nodesRed(1:nNodesRed,1:3),nNodesRed)
  call reduceElements()
  call writeReducedNodesAndElements2VTK()


  print*, ''
  print*, '--[ FINAL MESH, REDUCED PORES ]-----------------------------------------------------------------------------------------'
  
  call createPores()
  call connectUnconnectedPores()
  ! call removeUnconnectedPores()
  if(nNodesFinal.gt.0) then
    call writeFinalNodesAndElements2VTK()
  else
    ERROR STOP 'No pores connecting inlet and outlet!'
  endif
  ! endif
  
  contains

  subroutine createPosByRejection()
    integer :: nptclTemp,n
    real*8, allocatable :: posTemp(:,:)
    VSample = boxSize(1)*boxSize(2)*boxSize(3)
    nptclTemp = int((Vsample * rhoSample) / (Vtpod*rhoZnO))
    allocate(posTemp(nptclTemp,3))
    call random_number(posTemp)
    posTemp(:,1) = posTemp(:,1) * boxSize(1)
    posTemp(:,2) = posTemp(:,2) * boxSize(2)
    posTemp(:,3) = posTemp(:,3) * boxSize(3)
    nparticles = 0
    do i=1,nptclTemp
      if((shape.eq.1 .or. shape.eq.2) .and. acceptPosInvertedConeAxisZ(posTemp(i,:),rxIn,rxOut,ryIn,ryOut,centerIn,centerOut)) then
        nparticles = nparticles +1
        posTemp(nparticles,:) = posTemp(i,:)
      endif
    enddo
    allocate(pos(nparticles,3))
    pos = posTemp(1:nparticles,:)
    deallocate(posTemp)

    write(str1,*) 'Number of particles'; print'(2A,3X,I0,X,I0)', str1, '=', nParticles
  end subroutine

  subroutine createTpods()
    implicit none
    real*8 :: volumeTemp, lrr, min_lrr=5
    allocate(ang(nParticles,3))
    allocate(al(nParticles))
    allocate(ar(nParticles))
    do i=1,nparticles
      call randZYZangle(ang(i,1:3))
      do 
        call randn(volumeTemp,Vtpod,sd_volume)
        if(volumeTemp .gt. min_volume) exit
      enddo
      do 
        call randn(lrr,lr_ratio,sd*lr_ratio)
        if(lrr .gt. min_lrr) exit
      enddo
      ! al(i) = (0.25d0*lr_ratio**2*volumeTemp / PI)**(third) 
      ! ar(i) = al(i) / lr_ratio
      al(i) = (0.25d0*lrr**2*volumeTemp / PI)**(third) 
      ar(i) = al(i) / lrr
    enddo

    ! al(:) = tpoddim(1)
    ! ar(:) = tpoddim(2)
    ! --

    ! allocate and set central axes (dd) of tetrapod arms
      allocate(dd(0:nParticles,4,3))

      ! tetrapod directions with length 1 and angles 0
      dd(0,1,:) = (/ -0.5d0, -0.5d0, -0.5d0 /)
      dd(0,2,:) = (/  0.5d0,  0.5d0, -0.5d0 /)
      dd(0,3,:) = (/  0.5d0, -0.5d0,  0.5d0 /)
      dd(0,4,:) = (/ -0.5d0,  0.5d0,  0.5d0 /)
      dd(0,:,:) = dd(0,:,:) / norm2(dd(0,1,:))
      ! --

      ! rotation and scaling of tetrapod arms
        do i=1, nParticles
          do j=1,4
            dd(i,j,:) = al(i) * matmul(matmul(matmul(rotMatZ(ang(i,2)),rotMatY(ang(i,1))),rotMatZ(ang(i,3))),dd(0,j,:))
          enddo
        enddo
      ! --

  end subroutine

  subroutine getTubeIntersection(p1,p2,v1,v2,r1,r2,ip,a,b,info)
    use tensorLib
    implicit none
    real*8,dimension(3) :: p1,p2,v1,v2,distVec,ip
    real*8 :: r1,r2, dist, a, b, c
    real*8, dimension(2,2) :: M
    integer, dimension(2) :: ipiv
    integer :: info
    logical :: parallel
  
    info = 0
    parallel = .false.
  
    ! get the distance between the vector V1 going through point P1 and the vector V2 going through point P2
      distVec = getDistance(P1,P2,V1,V2)
      dist = norm2(distVec)
    ! --
  
    if(dist.le.r1+r2) then
      ! set up matrix for equation system to find intersection P1+a*V1 = P2+b*V2 => (V1 -V2) (a b)^T = P2-P1
        M(:,1) =  V1(1:2)
        M(:,2) = -V2(1:2)
        if(det(M).eq.0.d0) then
          M(:,1) =  V1(2:3)
          M(:,2) = -V2(2:3)
          if(det(M).eq.0.d0) then
            M(:,1) =  (/ V1(1), V1(3) /)
            M(:,2) = -(/ V2(1), V2(3) /)
            if(det(M).eq.0.d0) then
              parallel = .true.
            endif
          endif
        endif
      ! --
  
      ! if V1 and V2 are parallel
        if(parallel) then
          ! ip = 0.5d0*(P2+P1)
          ! a = dot_product(ip-P1,V1/norm2(V1))
          ! b = dot_product(ip-P2,V2/norm2(V2))
  
          ! if(a.ge.0.d0.and.a.le.1.d0.and.b.ge.0.d0.and.b.le.1.d0 &
          !             .or. a.ge.0.d0.and.a.le.0.5d0 .or. b.ge.0.d0.and.b.le.0.5d0) then
          !   info = 1
          ! endif
      ! --
  
      ! if V1 and V2 are not parallel
        else
          ip = P2-distVec-P1
          call dgesv(2,1,M,2,ipiv,ip(1:2),2,info)
          a = ip(1)
          b = ip(2)
          ip = P1 + a*V1 + 0.5d0*distVec
  
          if(a.ge.0.d0.and.a.le.1.d0.and.b.ge.0.d0.and.b.le.1.d0) then
            info = 2
          ! else
          !   c = (r1+r2)/norm2(V1-V2)
          !   a = a-c
          !   b = b-c
          !   if(a.ge.0d0.and.a.le.1.d0.and.b.ge.0.d0.and.b.le.1.d0) then
          !     ip = 0.5d0*(P1 + a*V1 + P2 + b*V2)
          !     info = 3
          !   else
          !     a = a+c
          !     b = b+c
          !     info=4
          !   endif
          endif
        endif
      ! --
    else
      info = 0
    endif
  end subroutine

  logical function acceptPosInvertedConeAxisZ(pos,rxIn,rxOut,ryIn,ryOut,centerIn,centerOut) result (accept)
    implicit none
    real*8 :: pos(3), rxIn(2), rxOut(2), ryIn(2), ryOut(2), aI, aO, bI, bO, rI, rO, lI, lO, h, centerIn(2), centerOut(2), phiIn
    real*8 :: phiOut, pIn(2), pOut(2)
    accept = .false.
    h = pos(3)/boxSize(3)
    pIn = pos(1:2) - centerIn
    pOut = pos(1:2) - centerOut
    aI = rxIn(1)  * h + rxIn(2)  * (1.d0-h)
    aO = rxOut(1) * h + rxOut(2) * (1.d0-h)
    bI = ryIn(1)  * h + ryIn(2)  * (1.d0-h)
    bO = ryOut(1) * h + ryOut(2) * (1.d0-h)
    if(abs(pIn(1)).gt.1.d-10) then
      phiIn = atan(pIn(2)/pIn(1))
    else
      phiIn = 0.5d0*pi
    endif
    if(abs(pOut(1)).gt.1.d-10) then
      phiOut = atan(pOut(2)/pOut(1))
    else
      phiOut = 0.5d0*pi
    endif
    rI = aI*bI / dsqrt((aI*sin(phiIn))**2 + (bI*cos(phiIn))**2)
    rO = aO*bO / dsqrt((aO*sin(phiOut))**2 + (bO*cos(phiOut))**2)
    lI = norm2(pos(1:2)-centerIn)
    lO = norm2(pos(1:2)-centerOut)
    if(lI.ge.rI .and. lO.le.rO) accept=.true.
  end function

  logical function isAtSurfaceInvertedConeAxisZ(pos,rx,ry,tol) result(isSurfIn)
    implicit none
    real*8 :: pos(3), rx(2), ry(2), tol, r, h, phi, a, b, l
    isSurfIn = .false.
    h = pos(3)/boxSize(3)
    if(pos(1).gt.abs(1.d-10)) then
      phi = atan(pos(2)/pos(1))
    else
      phi = 0.5d0*pi
    endif
    a = rx(1)  * h + rx(2)  * (1.d0-h)
    b = ry(1)  * h + ry(2)  * (1.d0-h)
    r = a*b / dsqrt((a*sin(phi))**2+(b*cos(phi))**2)
    l = norm2(pos(1:2)-centerIn)
    if(abs(l-r).lt.tol) isSurfIn = .true.
  end function

  subroutine findIntersectionsAndCreateNodesAndElements()
    implicit none
    integer :: a,aa,b,bb,i,p,pp, iNr(nParticles,4), elmtCMtemp(30*nparticles,2)
    real*8 :: ip(3), vs1, vs2, erTemp(30*nparticles)
    real*8 :: inter(nParticles*30,3), interpos(nParticles,4,30,2)

    interpos = 0.d0
    iNr = 0
    nInter = 0

    do p=1,nparticles
      if(outpt .gt. 0) write(*,'(A,I0,A,I0)',ADVANCE='NO') achar(13)//'connecting particle ',p,' of ', nParticles
      b = partBoxCM(p)
      ! do a=1,4
      !   nInter = nInter + 1
      !   iNr(p,a) = iNr(p,a) + 1
      !   inter(nInter,:) = pos(p,:) + dd(p,a,:)
      !   interpos(p,a,iNr(p,a),:) = (/ 1.d0, dble(nInter+nParticles)/)
      ! enddo
      do i=1,27
        bb = boxBoxCM(b,i)
        do j=1,nPartInBox(bb)
          pp = boxPartCM(bb,j)
          if (norm2(pos(p,:)-pos(pp,:)).le.al(p)+al(pp) .and. p.lt.pp) then
            do a=1,4
              do aa=1,4
                call getTubeIntersection(pos(p,:),pos(pp,:),dd(p,a,:),dd(pp,aa,:),ar(p),ar(pp),ip,vs1,vs2,info)
                select case(info)
                  case(0); !print*, 'no connection'
                  case(1); print*, 'direct connection'
                  case(2); !print*, 'case2 intersection at ip=', ip
                    nInter = nInter + 1
                    iNr(p,a) = iNr(p,a) + 1
                    iNr(pp,aa) = iNr(pp,aa) + 1
                    if(iNr(p,a).gt.30) then; print*, 'Need to increase array size'; ERROR STOP; endif
                    if(iNr(pp,aa).gt.30) then; print*, 'Need to increase array size'; ERROR STOP; endif
                    inter(nInter,:) = ip
                    interpos(p,a,iNr(p,a),:) = (/ vs1, dble(nInter+nParticles) /)
                    interpos(pp,aa,iNr(pp,aa),:) = (/ vs2, dble(nInter+nParticles) /)
                  case(3); !print*, 'case3 intersection at ip=', ip
                    ! nInter = nInter + 1
                    ! iNr(i,j) = iNr(i,j) + 1
                    ! iNr(k,l) = iNr(k,l) + 1
                    ! if(iNr(i,j).gt.10) then; print*, 'Need to increase array size'; ERROR STOP; endif
                    ! if(iNr(k,l).gt.10) then; print*, 'Need to increase array size'; ERROR STOP; endif
                    ! inter(nInter,:) = ip
                    ! interpos(i,j,iNr(i,j),:) = (/ vs1, dble(nInter+nParticles) /)
                    ! interpos(k,l,iNr(k,l),:) = (/ vs2, dble(nInter+nParticles) /)
                  case(4);
                end select
              enddo
            enddo
          endif
        enddo
      enddo
    enddo
    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') achar(13)
    write(str1,*) 'Number of intersection points'; print'(2A,G15.7)', str1, '=', nInter

    nNodes = nParticles + nInter
    write(str1,*) 'Number of nodes'; print'(2A,G15.7)', str1, '=', nNodes
    allocate(nodes(nNodes,3))
    nodes = 0.d0
    nodes(1:nParticles,1:3) = pos(1:nparticles,1:3)
    if(nInter.gt.0) nodes(nParticles+1:nNodes,1:3) = inter(1:nInter,1:3)

    nElmt = 0
    do i=1,nParticles
      do j=1,4
        if(iNr(i,j).gt.0) then       
          call bubbleSort(transpose(interpos(i,j,1:iNr(i,j),:)),1)
          nElmt= nElmt+1
          elmtCMtemp(nElmt,:) = (/ i, int(interpos(i,j,1,2)) /)
          erTemp(nElmt) = ar(i)
          do k=1,iNr(i,j)-1
            nElmt = nElmt+1
            elmtCMtemp(nElmt,:) = (/ int(interpos(i,j,k,2)), int(interpos(i,j,k+1,2)) /)
            erTemp(nElmt) = ar(i)
          enddo
        endif
      enddo
    enddo

    ! label and count unconnected nodes
      allocate(nodeCM(nNodes))
      nNodesRed = 0
      nodeCM = -1
      do i=1, nElmt
        nodeCM(elmtCMtemp(i,1)) = 0
        nodeCM(elmtCMtemp(i,2)) = 0
      enddo
      nNodesUnconnected = 0
      do i=1,nNodes
        if(nodeCM(i).eq.-1) nNodesUnconnected = nNodesUnconnected + 1 
        ! if(nodeCM(i).eq.-1) call connectNodeToClosestNode(i)
      enddo
      write(str1,*) 'Number of unconnected nodes'; print'(2A,G15.7)', str1, '=', nNodesUnconnected
    ! -- 

    write(str1,*) 'Number of elements'; print'(2A,G15.7)', str1, '=', nElmt

    allocate(elmtCM(nElmt+nNodesUnconnected,2))
    allocate(er(nElmt+nNodesUnconnected))
    elmtCM(1:nElmt,1:2) = elmtCMtemp(1:nElmt,1:2)
    er(1:nElmt) = erTemp(1:nElmt)

  end subroutine

  subroutine combineNodeClusters()
    implicit none
    integer :: b,n,nn,nCluster
    real*8 :: nodeTemp(3)

    allocate(nodesRed(nNodes,3))
    nNodesRed = 0
    do n=1, nNodes
      if(outpt .gt. 0) write(*,'(A,I0,A,I0)',ADVANCE='NO') achar(13)//'Check for cluster node ',n,' of ', nNodes
      if(nodeCM(n).eq.0) then
        nNodesRed = nNodesRed + 1
        nodeCM(n) = nNodesRed
        nodeTemp = nodes(n,:)
        nCluster = 1
        b = nodeBoxCM(n)
        do i=1,nNodeInBox(b)
          nn = boxNodeCM(b,i)
          if(n.lt.nn .and. nodeCM(nn).eq.0 .and. norm2(nodes(n,:)-nodes(nn,:)).lt.2.d0*ar(1)) then
            nodeTemp = nodeTemp + nodes(nn,:)
            nodeCM(nn) = nNodesRed
            nCluster = nCluster+1
          endif
        enddo
        nodesRed(nNodesRed,:) = nodeTemp / dble(nCluster)
      endif
    enddo
    deallocate(nodes)
    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') achar(13)
    write(str1,*) 'Reduced number of nodes'; print'(2A,G15.7)', str1, '=', nNodesRed
  end subroutine

  subroutine reduceElements()
    implicit none
    integer :: b1, b2, n1, n2, nn1, nn2, e, ee, maxElmt, n(2), idx
    integer, allocatable :: elmtRedCMtemp(:,:), boxElmtCMtemp(:,:)
    logical :: newElmt, allNodes(nNodesRed)
    real*8 :: erRedTemp(nElmt)

    allocate(boxElmtCMtemp(nBoxes,nElmt*100/nBoxes))
    allocate(nElmtInBox(nBoxes))
    allocate(elmtRedCMtemp(nElmt,3))
    elmtRedCMtemp = 0
    nElmtRed = 0
    nElmtInBox = 0
    do e=1,nElmt
      if(outpt .gt. 0) write(*,'(A,I0,A,I0)',ADVANCE='NO') achar(13)//'Check to remove elmt ',e,' of ', nElmt
      newElmt = .true.
      n = (/ nodeCM(elmtCM(e,1)), nodeCM(elmtCM(e,2)) /)
      n1 = minval(n)
      n2 = maxval(n)
      b1 = nodeBoxCM(n1)
      b2 = nodeBoxCM(n2)
      if(n1 .eq. n2) then
        newElmt = .false.
      else
        do i=1,nElmtInBox(b1)
          ee = boxElmtCMtemp(b1,i)
          nn1 = elmtRedCMtemp(ee,1)
          nn2 = elmtRedCMtemp(ee,2)
          if(n1.eq.nn1 .and. n2.eq.nn2) then
            newElmt = .false.
            exit
          endif
        enddo
        if(newElmt .and. b1.ne.b2) then
          do i=1,nElmtInBox(b2)
            ee = boxElmtCMtemp(b2,i)
            nn1 = elmtRedCMtemp(ee,1)
            nn2 = elmtRedCMtemp(ee,2)
            if(n1.eq.nn1 .and. n2.eq.nn2) then
              newElmt = .false.
              exit
            endif
          enddo
        endif

        if(newElmt) then
          nElmtRed = nElmtRed + 1
          elmtRedCMtemp(nElmtRed,1) = n1
          elmtRedCMtemp(nElmtRed,2) = n2
          erRedTemp(nElmtRed) = er(e)

          nElmtInBox(b1) = nElmtInBox(b1) + 1
          boxElmtCMtemp(b1,nElmtInBox(b1)) = nElmtRed
          if(b1.ne.b2) then
            nElmtInBox(b2) = nElmtInBox(b2) + 1
            boxElmtCMtemp(b2,nElmtInBox(b2)) = nElmtRed
          endif
        endif
      endif
    enddo

    allNodes = .false.
    do i=1,nElmtRed
      allNodes(elmtRedCMtemp(i,1)) = .true.
      allNodes(elmtRedCMtemp(i,2)) = .true.
    enddo
    do i=1,nNodesRed
      if(.not.allNodes(i)) then
        ! print*, 'Node ', i, 'not in elmtRedCM'
        call getClosestNode(i,idx,nodesRed(1:nNodesRed,1:3),nNodesRed)
        nElmtRed = nElmtRed + 1
        elmtRedCMtemp(nElmtRed,1:2) = (/ i, idx /)
        erRedTemp(nElmtRed) = tpoddim(2)
        nodeCM(i) = 0
      endif
    enddo

    deallocate(elmtCM)
    deallocate(er)
    allocate(elmtRedCM(nElmtRed,3))
    allocate(erRed(nElmtRed))
    elmtRedCM = elmtRedCMtemp(1:nElmtRed,:)
    erRed = erRedTemp(1:nElmtRed)
    deallocate(elmtRedCMtemp)
    maxElmt = maxval(nElmtInBox)
    allocate(boxElmtCM(nBoxes,maxElmt))
    boxElmtCM(1:nBoxes,1:maxElmt) = boxElmtCMtemp(1:nBoxes,1:maxElmt)
    deallocate(boxElmtCMtemp)
    deallocate(nodeCM)

    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') achar(13)
    write(str1,*) 'Reduced number of elements'; print'(2A,G15.7)', str1, '=', nElmtRed 
  end subroutine

  subroutine createPores()
    implicit none
    integer :: e, ee, n1, n2, b1, b2, nn1, nn2, maxElmt
    integer :: nElmtInPoreTemp(nElmtRed), porePoreCM(nElmtRed), ElmtElmtCM(nElmtRed,20), nElmtAtElmt(nElmtRed)
    integer :: minIdx
    logical :: repeat
        
    do e=1,nElmtRed; elmtRedCM(e,3) = e; enddo
    nElmtAtElmt = 0
    do e=1, nElmtRed
      if(outpt .gt. 0) write(*,'(A,I0,A,I0,A)',ADVANCE='NO') achar(13)//'Assigning elmt ',e,' of ', nElmtRed, ' to pore'

      n1 = elmtRedCM(e,1)
      n2 = elmtRedCM(e,2)
      b1 = nodeBoxCM(n1)
      b2 = nodeBoxCM(n2)
      do i=1, nElmtInBox(b1)
        ee = boxElmtCM(b1,i)
        if(ee.gt.nElmtREd) print*, 'ee > nElmtRed', ee, nElmtred
        if(e.lt.ee) then
          nn1 = elmtRedCM(ee,1)
          nn2 = elmtRedCM(ee,2)
          if(n1.eq.nn1 .or. n1.eq.nn2 .or. n2.eq.nn1 .or. n2.eq.nn2) then
            nElmtAtElmt(e) = nElmtAtElmt(e) + 1
            nElmtAtElmt(ee) = nElmtAtElmt(ee) + 1
            elmtElmtCM(e,nElmtAtElmt(e)) = ee
            elmtElmtCM(ee,nElmtAtElmt(ee)) = e
          endif
        endif
      enddo
      if(b1.ne.b2) then
        do i=1, nElmtInBox(b2)
          ee = boxElmtCM(b2,i)
          if(ee.gt.nElmtREd) print*, 'ee > nElmtRed', ee, nElmtred
          if(e.lt.ee) then
            nn1 = elmtRedCM(ee,1)
            nn2 = elmtRedCM(ee,2)
            if(n1.eq.nn1 .or. n1.eq.nn2 .or. n2.eq.nn1 .or. n2.eq.nn2) then
              nElmtAtElmt(e) = nElmtAtElmt(e) + 1
              nElmtAtElmt(ee) = nElmtAtElmt(ee) + 1
              elmtElmtCM(e,nElmtAtElmt(e)) = ee
              elmtElmtCM(ee,nElmtAtElmt(ee)) = e
            endif
          endif
        enddo
      endif
    enddo

    do
      repeat = .false. 
      do e=1,nElmtRed
        minIdx = elmtRedCM(e,3)
        do i=1,nElmtAtElmt(e)
          if(elmtRedCM(elmtElmtCM(e,i),3).lt.minIdx) then
            minIdx = elmtRedCM(elmtElmtCM(e,i),3)
            repeat = .true.
          endif
        enddo
        elmtRedCM(e,3) = minIdx
        do i=1,nElmtAtElmt(e)
          elmtRedCM(elmtElmtCM(e,i),3) = minIdx
        enddo
      enddo
      if(.not.repeat) exit
    enddo

    porePoreCM = 0
    ! do i=1,nPores; porePoreCM(i) = i; enddo
    nElmtInPoreTemp = 0
    do e=1, nElmtRed
      nElmtInPoreTemp(elmtRedCM(e,3))  = nElmtInPoreTemp(elmtRedCM(e,3)) + 1
    enddo
    nPores = 0
    do e=1, nElmtRed
      if(nElmtInPoreTemp(e).gt.0)then
        nPores = nPores+1
        porePoreCM(e) = nPores
      endif
    enddo
    maxElmt = maxval(nElmtInPoreTemp)
    allocate(poreElmtCM(nPores,maxElmt))
    allocate(nodePoreCM(nNodesRed))
    allocate(nElmtInPore(nPores))
    nElmtInPore = 0
    do e=1, nElmtRed
      elmtRedCM(e,3) = porePoreCM(elmtRedCM(e,3))
      nElmtInPore(elmtRedCM(e,3)) = nElmtInPore(elmtRedCM(e,3)) +1 
      poreElmtCM(elmtRedCM(e,3),nElmtInPore(elmtRedCM(e,3))) = e
      nodePoreCM(elmtRedCM(e,1)) = elmtRedCM(e,3)
      nodePoreCM(elmtRedCM(e,2)) = elmtRedCM(e,3)
    enddo

    do i=1, nNodesREd
      if(nodePoreCM(i).gt.nPores .or. nodePoreCM(i).lt.1) then
        print*, 'Node', i, 'not correctly connected', nodePoreCM(i)
        read(*,*)
      endif
    enddo

    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') achar(13)//dl
    write(str1,*) 'Number of pores'; print'(2A,G15.7)', achar(13)//str1, '=', nPores

  end subroutine

  subroutine connectUnconnectedPores()
    implicit none
    integer :: i,j, idx1,idx2,nElmt2Add,elmtAddCM(nPores,2),nNodesSurf(2)
    logical :: poreSurf(nPores,2)

    poreSurf = .false.

    nNodesSurf = 0
    do i=1, nNodesRed
      if(outpt .gt. 0) write(*,'(A,I0,A,I0)',ADVANCE='NO') achar(13)//'Finding surface nodes ',i,' of ', nNodesRed
      if((shape.eq.1 .or. shape.eq.2) .and. isAtSurfaceInvertedConeAxisZ(nodesRed(i,:),rxIn,ryIn,tpoddim(1))) then
        poreSurf(nodePoreCM(i),1) = .true.
        nNodesSurf(1) = nNodesSurf(1) + 1
      elseif((shape.eq.1 .or. shape.eq.2) .and. isAtSurfaceInvertedConeAxisZ(nodesRed(i,:),rxOut,ryOut,tpoddim(1))) then  
        poreSurf(nodePoreCM(i),2) = .true.
        nNodesSurf(2) = nNodesSurf(2) + 1
      endif
    enddo

    nElmt2Add = 0
    write(str1,*) 'Number of surface nodes'; print'(2A,2G15.7)', str1, '=', nNodesSurf(1), nNodesSurf(2)
    do i=1, nPores
      if(poreSurf(i,1).and.poreSurf(i,2)) then
        nPoresRed = nPoresRed +1 
      else
        idx1 = elmtRedCM(poreElmtCM(i,1),1)
        call getClosestNodeOfDifferentPore(idx1,idx2)
        nElmt2Add = nElmt2Add+1
        elmtAddCM(nElmt2Add,1:2) = (/ idx1,idx2 /)
        
        poreSurf(i,1) = poreSurf(i,1) .or.  poreSurf(nodePoreCM(idx2),1)
        poreSurf(i,2) = poreSurf(i,2) .or.  poreSurf(nodePoreCM(idx2),2)
        poreSurf(nodePoreCM(idx2),1) = poreSurf(i,1) .or.  poreSurf(nodePoreCM(idx2),1)
        poreSurf(nodePoreCM(idx2),2) = poreSurf(i,2) .or.  poreSurf(nodePoreCM(idx2),2)
        do j=1,nElmtInPore(i)
          nodePoreCM(elmtRedCM(poreElmtCM(i,j),1)) = nodePoreCM(idx2)
          nodePoreCM(elmtRedCM(poreElmtCM(i,j),2)) = nodePoreCM(idx2)
          elmtRedCM(poreElmtCM(i,j),3) = nodePoreCM(idx2)
        enddo
      endif
    enddo
    do i=1,nPores
      ! print*, poreSurf(i,:)
    enddo

    
    nNodesFinal = nNodesRed
    nElmtFinal = nElmtRed+nElmt2Add

    allocate(nodes(nNodesFinal,3))
    allocate(elmtCM(nElmtFinal,2))
    allocate(er(nElmtFinal))


    elmtCM(1:nElmtRed,1:2) = elmtRedCM(1:nElmtRed,1:2)
    elmtCM(nElmtRed+1:nElmtFinal,1:2) = elmtAddCM(1:nElmt2Add,1:2)
    er(1:nElmtRed) = erRed(1:nElmtRed)
    er(nElmtRed+1:nElmtFinal) = tpoddim(2)
    nodes = nodesRed

    write(str1,*) 'Number of final nodes'; print'(2A,G15.7)', str1, '=', nNodesFinal
    write(str1,*) 'Number of final elements'; print'(2A,G15.7)', str1, '=', nElmtFinal


  end subroutine

  subroutine removeUnconnectedPores()
    implicit none
    integer :: idxs(nPores), idxs2(nElmtRed), idxs3(2*nElmtRed), nPores2Remove, nNodes2Remove, nElmt2Remove
    integer :: i,j,k, nNodesSurf(2)
    logical :: poreSurf(nPores,2)

    poreSurf = .false.
    nNodesSurf = 0
    do i=1, nNodesRed
      if(outpt .gt. 0) write(*,'(A,I0,A,I0)',ADVANCE='NO') achar(13)//'Finding surface nodes ',i,' of ', nNodesRed
      if((shape.eq.1 .or. shape.eq.2) .and. isAtSurfaceInvertedConeAxisZ(nodesRed(i,:),rxIn,ryIn,tpoddim(1))) then
        poreSurf(nodePoreCM(i),1) = .true.
        nNodesSurf(1) = nNodesSurf(1) + 1
      elseif((shape.eq.1 .or. shape.eq.2) .and. isAtSurfaceInvertedConeAxisZ(nodesRed(i,:),rxOut,ryOut,tpoddim(1))) then  
        poreSurf(nodePoreCM(i),2) = .true.
        nNodesSurf(2) = nNodesSurf(2) + 1
      endif
    enddo

    nPoresRed = 0
    nPores2Remove = 0
    idxs = 0
    do i=1, nPores
      if(poreSurf(i,1).and.poreSurf(i,2)) then
        nPoresRed = nPoresRed +1 
      else
        nPores2Remove = nPores2Remove+1
        idxs(nPores2Remove) = i
      endif
    enddo

    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') achar(13)
    write(str1,*) 'Number of surface nodes'; print'(2A,2(I0,2X))', str1, '=', nNodesSurf
    write(str1,*) 'Reduced number of pores'; print'(2A,G15.7)', str1, '=', nPoresRed

    nElmt2Remove = 0
    idxs2 = 0 
    idxs3 = 0
    do i=1, nPores2Remove
      do j=1, nElmtInPore(idxs(i))
        nElmt2Remove = nElmt2Remove+1
        idxs2(nElmt2Remove) = poreElmtCM(idxs(i),j)
        idxs3(nElmt2Remove*2-1:nElmt2Remove*2) = elmtRedCM(idxs2(nElmt2Remove),1:2)
      enddo
    enddo

    call quicksortInt(idxs2,1,nElmt2Remove)
    call quicksortInt(idxs3,1,2*nElmt2Remove)

    nNodes2Remove = 1
    do i=2,2*nElmt2Remove
      if(idxs3(i-1).lt.idxs3(i)) then
        nNodes2Remove = nNodes2Remove + 1
        idxs3(nNodes2Remove) = idxs3(i)
      endif
    enddo

    nNodesFinal = nNodesRed-nNodes2Remove
    nElmtFinal = nElmtRed-nElmt2Remove

    allocate(nodes(nNodesFinal,3))
    allocate(elmtCM(nElmtFinal,2))
    allocate(nodeCM(nNodesRed))

    k=0
    j=1
    nodeCM = 0
    do i=1,nNodesRed
      if(i.ne.idxs3(j)) then
        k = k+1
        nodes(k,1:3) = nodesRed(i,1:3)
        nodeCM(i) = k
      else
        j = j+1
      endif
    enddo

    k=0
    j=1
    do i=1,nElmtRed
      if(i.ne.idxs2(j)) then
        k = k+1
        elmtCM(k,1) = nodeCM(elmtRedCM(i,1))
        elmtCM(k,2) = nodeCM(elmtRedCM(i,2))
        if(nodeCM(elmtRedCM(i,1)).le.0 .or. nodeCM(elmtRedCM(i,1)).gt.nNodesFinal .or. &
           nodeCM(elmtRedCM(i,1)).le.0 .or. nodeCM(elmtRedCM(i,1)).gt.nNodesFinal) then
          print*, 'nodeCM not correct'
        endif
      else
        j = j+1
      endif
    enddo

    write(str1,*) 'Number of final nodes'; print'(2A,G15.7)', str1, '=', nNodesFinal
    write(str1,*) 'Number of final elements'; print'(2A,G15.7)', str1, '=', nElmtFinal

  end subroutine

  subroutine connectUnconnectedNodes()
    implicit none
    integer :: i,idx
    do i=1,nNodes
      ! if(nodeCM(i).eq.-1) nNodesUnconnected = nNodesUnconnected + 1 
      if(nodeCM(i).eq.-1) then
        call getClosestNode(i,idx,nodes(1:nNodes,1:3),nNodes)
        nElmt = nElmt + 1
        elmtCM(nElmt,1:2) = (/ i, idx /)
        er(nElmt) = tpoddim(2)
        nodeCM(i) = 0
      endif
    enddo
  end subroutine

  subroutine getClosestNode(idx1,minDistIdx, nodes, nNodes)
    implicit none
    integer :: i,j,idx1,idx2,bb,b,minDistIdx, nNodes
    real*8 :: minDist, dmmy, nodes(nNodes,3)
    minDist = boxSize(1) * boxSize(2) * boxSize(3)
    b = getBoxNrByPos(nodes(idx1,:))
    do i=1,27
      bb = boxBoxCM(b,i)
      do j=1,nNodeInBox(bb)
        idx2 = boxNodeCM(bb,j)
        dmmy = norm2(nodes(idx1,:)-nodes(idx2,:))
        if(idx1.ne.idx2 .and. dmmy.lt.minDist) then
          minDist = dmmy
          minDistIdx = idx2
        endif
      enddo
    enddo
  end subroutine

  subroutine getClosestNodeOfDifferentPore(idx1,minDistIdx)
    implicit none
    integer :: i,j,idx1,idx2,minDistIdx
    real*8 :: minDist, dmmy
    minDist = boxSize(1) * boxSize(2) * boxSize(3)
    do i=1,nNodesRed
      idx2 = i
      dmmy = norm2(nodesRed(idx1,:)-nodesRed(idx2,:))
      if(nodePoreCM(idx1).ne.nodePoreCM(idx2) .and. dmmy.lt.minDist) then
        minDist = dmmy
        minDistIdx = idx2
      endif
    enddo
  end subroutine

  ! BOXING ROUTINES ----------------------------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------------------------------------------

  subroutine createBoxBoxCM()
    implicit none
    integer :: xyz(3), i,j,k,l
    allocate(boxBoxCM(nBoxes,27))
    boxBoxCM(:,:) = 0
    do i=1,nBoxes
      xyz = getBoxXYZ(i)
      if(isCenterBoxXYZ(xyz)) then
        do j=-1,1
          do k=-1,1
            do l=-1,1
              boxBoxCM(i,(j+2)+(k+1)*3+(l+1)*9) = getBoxNrByXYZ(xyz+(/j,k,l/))
            enddo
          enddo
        enddo
      endif
    enddo
  end subroutine

  subroutine createBoxPartCM()
    implicit none
    integer :: bn, maxPart
    integer, allocatable :: boxPartCMtemp(:,:)
    allocate(partBoxCM(nParticles))
    allocate(boxPartCMtemp(nBoxes,nParticles*30/nBoxes))
    allocate(nPartInBox(nBoxes))
    nPartInBox = 0
    do i=1,nparticles
      bn = getBoxNrByPos(pos(i,:))
      nPartInBox(bn) = nPartInBox(bn) + 1
      boxPartCMtemp(bn,nPartInBox(bn)) = i
      partBoxCM(i) = bn
    enddo
    maxPart = maxval(nPartInBox)
    allocate(boxPartCM(nBoxes,maxPart))
    boxPartCM = boxPartCMtemp(:,1:maxPart)
    deallocate(boxPartCMtemp)
  end subroutine

  subroutine createBoxNodeCM(nodes,nNodes)
    implicit none
    integer :: nNodes, bn, n, maxNode
    real*8 :: nodes(nNodes,3)
    integer, allocatable :: boxNodeCMtemp(:,:)
    ! print*, nBoxes, nNodes, size(nodes,1), nNodesRed
    maxNode = nNodes*100/nBoxes
    allocate(nodeBoxCM(nNodes))
    allocate(nNodeInBox(nBoxes))
    allocate(boxNodeCMtemp(nBoxes,maxNode))

    nNodeInBox = 0
    do n=1,nNodes
      bn = getBoxNrByPos(nodes(n,:))
      nNodeInBox(bn) = nNodeInBox(bn) + 1
      if(nNodeInBox(bn).gt.maxNode) then
        print*, 'nNodeInBox>maxNode', nNodeInBox(bn), maxNode
        ERROR STOP
      endif
      boxNodeCMtemp(bn,nNodeInBox(bn)) = n
      nodeBoxCM(n) = bn
    enddo
    maxNode = maxval(nNodeInBox)
    allocate(boxNodeCM(nBoxes,maxNode))
    boxNodeCM = boxNodeCMtemp(:,1:maxNode)
    deallocate(boxNodeCMtemp)
  end subroutine

  subroutine createBoxElmtCM()
    implicit none
    integer :: bn1, bn2, e, maxElmt
    integer, allocatable :: boxElmtCMtemp(:,:)
    maxElmt = nElmt*30/nBoxes
    allocate(boxElmtCMtemp(nBoxes,maxElmt))
    allocate(nElmtInBox(nBoxes))
    boxElmtCMtemp = 0
    nElmtInBox = 0
    do e=1, nElmt
      bn1 = nodeBoxCM(elmtCM(e,1))
      bn2 = nodeBoxCM(elmtCM(e,2))
      nElmtInBox(bn1) = nElmtInBox(bn1) + 1
      if(nElmtInBox(bn1).gt.maxElmt) then
        print*, 'nElmtInBox > maxElmt', nElmtInBox(bn2), maxElmt
        ERROR STOP
      endif
      boxElmtCMtemp(bn1,nElmtInBox(bn1)) = e
      if(bn2.ne.bn1) then
        nElmtInBox(bn2) = nElmtInBox(bn2) + 1
        if(nElmtInBox(bn2).gt.maxElmt) then
          print*, 'nElmtInBox > maxElmt', nElmtInBox(bn2), maxElmt
          ERROR STOP
        endif
        boxElmtCMtemp(bn2,nElmtInBox(bn2)) = e
      endif
    enddo
    maxElmt = maxval(nElmtInBox)
    allocate(boxElmtCM(nBoxes,maxElmt))
    boxElmtCM(1:nBoxes,1:maxElmt) = boxElmtCM(1:nBoxes,1:maxElmt)
    deallocate(boxElmtCMtemp)
  end subroutine

  integer function getBoxNrByPos(pos) result(bn)
    implicit none
    real*8 :: pos(3)
    integer :: x,y,z
    x = floor(pos(1)/boxSize(1)*(nBoxesX-2))+2
    y = floor(pos(2)/boxSize(2)*(nBoxesY-2))+2
    z = floor(pos(3)/boxSize(3)*(nBoxesZ-2))+2
    bn = x + nBoxesX * (y-1) + nBoxesX * nBoxesY * (z-1)
    if(bn .gt. nBoxes) then
      print*, 'Box number not correct: nBoxes=', nBoxes, 'bn=', bn
      print*, boxSize, nBoxesX,nBoxesY,nBoxesZ
      print*, pos, x, y, z
      ERROR STOP 
    endif
  end function

  integer function getBoxNrByXYZ(xyz) result(bn)
    implicit none
    integer :: xyz(3)
    bn = xyz(1) + nBoxesX * (xyz(2)-1) + nBoxesX * nBoxesY * (xyz(3)-1)
  end function

  function getBoxXYZ(bn) result(xyz)
    integer :: bn, xyz(3)
    xyz(3) = (bn-1) / (nBoxesX*nBoxesY)
    xyz(2) = (bn - xyz(3) * nBoxesX * nBoxesY - 1) / nBoxesX
    xyz(1) = bn - xyz(3) * nBoxesX * nBoxesY - xyz(2) * nBoxesX - 1
    xyz = xyz + 1
  end function

  logical function isCenterBoxXYZ(xyz)
    implicit none
    integer :: xyz(3),bn
    isCenterBoxXYZ = (xyz(1)-1)*(xyz(2)-1)*(xyz(3)-1).ne.0 .and. (xyz(1)-nBoxesX)*(xyz(2)-nBoxesY)*(xyz(3)-nBoxesZ).ne.0
  end function 


  ! SORTING ROUTINES ---------------------------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------------------------------------------

  recursive subroutine quicksort(a, first, last)
    implicit none
    real*8 :: a(*), x, t
    integer :: first, last
    integer :: i, j

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
      do while (a(i) < x)
          i=i+1
      end do
      do while (x < a(j))
          j=j-1
      end do
      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      i=i+1
      j=j-1
    end do
    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
  end subroutine quicksort

  recursive subroutine quicksortInt(a, first, last)
    implicit none
    integer :: a(*), x, t
    integer :: first, last
    integer :: i, j

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
      do while (a(i) < x)
          i=i+1
      end do
      do while (x < a(j))
          j=j-1
      end do
      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      i=i+1
      j=j-1
    end do
    if (first < i-1) call quicksortInt(a, first, i-1)
    if (j+1 < last)  call quicksortInt(a, j+1, last)
  end subroutine quicksortInt


  ! I/O ROUTINES -------------------------------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------------------------------------------

  subroutine readInputArgs()
    implicit none
    nargs = command_argument_count();
    i=0;
    do; i=i+1; if(i .gt. nargs) exit
      call get_command_argument(i,arg)
      select case(arg)
        case('-shape')
          i=i+1; call get_command_argument(i,arg); read(arg,*) shape
        case('-rxIn')
          do j=1,2; i=i+1; call get_command_argument(i,arg); read(arg,*) rxIn(j); enddo
        case('-rxOut')
          do j=1,2; i=i+1; call get_command_argument(i,arg); read(arg,*) rxOut(j); enddo
        case('-height')
          i=i+1; call get_command_argument(i,arg); read(arg,*) boxSize(3)
        case('-matlab')
          outpt = 0
        case('-tpoddim')
          do j=1,2; i=i+1; call get_command_argument(i,arg); read(arg,*) tpoddim(j); enddo
        case('-rhosample')
          i=i+1; call get_command_argument(i,arg); read(arg,*) rhoSample
        case('-centershift')
          do j=1,2; i=i+1; call get_command_argument(i,arg); read(arg,*) centerShift(j); enddo
        case('-vtpodfact')
          i=i+1; call get_command_argument(i,arg); read(arg,*) Vtpodfact
        case('-sd')
          i=i+1; call get_command_argument(i,arg); read(arg,*) sd
        case default
          call get_command_argument(i,arg)
          print*, 'Unknown argument: '//arg
          print*, 'To continue press any key'
          read(*,*)
        endselect
      enddo
      if(shape.eq.1) then
        boxSize(1:2) = (/ maxval(rxOut(1:2)), maxval(rxOut(1:2)) /)
        centerIn = (/ 0.d0, 0.d0 /) 
        centerOut = (/ 0.d0, 0.d0 /)
      elseif(shape.eq.2) then
        ryIn = rxIn
        ryOut = rxOut
        boxSize(1) = 2.d0 * maxval(rxOut(1:2))
        boxSize(2) = 2.d0 * maxval(ryOut(1:2))
        centerIn = (/ maxval(rxOut(1:2)), maxval(ryOut(1:2))/) + centerShift
        centerOut = (/ maxval(rxOut(1:2)), maxval(ryOut(1:2))/)
      endif
      ! Vtpod = 4.d0 * tpoddim(1) * tpoddim(2)**2 * PI * VtpodFact
      Vtpod = 12.d0 * cos(PI/6.d0) * tpoddim(1) * tpoddim(2)**2 * VtpodFact
      lr_ratio = tpoddim(1) / tpoddim(2)
      sd_volume = sd*Vtpod
      min_volume = 0.1d0*Vtpod
      ryIn = rxIn
      ryOut = rxOut
      ! print*, 'shape   = ', shape
      ! print*, 'rxIn    = ', rxIn
      ! print*, 'rxOut   = ', rxOut
      ! print*, 'height  = ', boxSize(3)
      ! print*, 'tpoddim = ', tpoddim
      ! print*, 'rho     = ', rhoSample
      ! print*, 'shift   = ', centershift
      ! print*, 'vfct    = ', vtpodfact
      ! print*, 'sd      = ', sd
      ! read(*,*)
  end subroutine

  subroutine writeTetrapods2VTK()
    implicit none
    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') 'Writing tetrapods to file'
    open(unit=2,file='tetrapods.vtk',access='stream',status='replace')
    write(str,'(A)') '# vtk DataFile Version 3.0'//lf; write(2) trim(str)
    write(str,'(A)') '1D Diffusion Mesh'//lf; write(2) trim(str)
    write(str,'(A)') 'ASCII'//lf//lf; write(2) trim(str)
    write(str,'(A)') 'DATASET POLYDATA'//lf; write(2)trim(str)
    write(str,'(A,I0,A)') 'POINTS ',5*nParticles,' DOUBLE'//lf; write(2)trim(str)
    do i=1,nParticles
      write(str,'(3F9.2,A)') pos(i,:),lf; write(2) trim(str)
      write(str,'(3F9.2,A)') pos(i,:)+dd(i,1,:),lf; write(2) trim(str)
      write(str,'(3F9.2,A)') pos(i,:)+dd(i,2,:),lf; write(2) trim(str)
      write(str,'(3F9.2,A)') pos(i,:)+dd(i,3,:),lf; write(2) trim(str)
      write(str,'(3F9.2,A)') pos(i,:)+dd(i,4,:),lf; write(2) trim(str)
    enddo
    write(str,'(A,2(I0,X),A)') lf//'LINES ',nParticles*4,nParticles*12,lf; write(2) trim(str)
    do i=1,nParticles
      write(str,'(3(I0,X),A)') 2, (i-1)*5, (i-1)*5+1, lf; write(2) trim(str)
      write(str,'(3(I0,X),A)') 2, (i-1)*5, (i-1)*5+2, lf; write(2) trim(str)
      write(str,'(3(I0,X),A)') 2, (i-1)*5, (i-1)*5+3, lf; write(2) trim(str)
      write(str,'(3(I0,X),A)') 2, (i-1)*5, (i-1)*5+4, lf; write(2) trim(str)
    enddo
    write(str,'(A,I0,A)') lf//'CELL_DATA ',nParticles*4,lf; write(2) trim(str)
    write(str,'(A)') 'SCALARS radius FLOAT'//lf; write(2) trim(str)
    write(str,'(A)') 'LOOKUP_TABLE default'//lf; write(2) trim(str)
    do i=1,nParticles
      write(str,'(F8.6,A)') ar(i), lf; write(2) trim(str)
      write(str,'(F8.6,A)') ar(i), lf; write(2) trim(str)
      write(str,'(F8.6,A)') ar(i), lf; write(2) trim(str)
      write(str,'(F8.6,A)') ar(i), lf; write(2) trim(str)
    enddo
    close(2)
    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') achar(13)//dl
    print*, achar(13)//' Created tetrapod file' 
  end subroutine

  subroutine writeNodesAndElements2VTK()
    implicit none
    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') 'Writing nodes and elements to vtk file'
    open(unit=1,file='mesh.vtk',access='stream',status='replace')
    write(str,'(A)') '# vtk DataFile Version 3.0'//lf; write(1) trim(str)
    write(str,'(A)') '1D Diffusion Mesh'//lf; write(1) trim(str)
    write(str,'(A)') 'ASCII'//lf//lf; write(1) trim(str)
    write(str,'(A)') 'DATASET POLYDATA'//lf; write(1)trim(str)
    write(str,'(A,I0,A)') 'POINTS ',nNodes,' DOUBLE'//lf; write(1)trim(str)
    do i=1,nNodes
      write(str,'(3F9.2,A)') nodes(i,1:3),lf; write(1) trim(str)
    enddo
    write(str,'(A,2(I0,X),A)') lf//'LINES ',nElmt,nElmt*3,lf; write(1) trim(str)
    do i=1,nElmt
      write(str,'(3(I0,X),A)') 2, elmtCM(i,1)-1, elmtCM(i,2)-1, lf; write(1) trim(str)
    enddo
    write(str,'(A,I0,A)') lf//'CELL_DATA ',nElmt,lf; write(1) trim(str)
    write(str,'(A)') 'SCALARS radius FLOAT'//lf; write(1) trim(str)
    write(str,'(A)') 'LOOKUP_TABLE default'//lf; write(1) trim(str)
    do i=1,nElmt
      write(str,'(F8.6,A)') er(i), lf; write(1) trim(str)
    enddo
    close(1)
    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') achar(13)//dl
    print*, achar(13)//' Created mesh file' 
  end subroutine

  subroutine writeReducedNodesAndElements2VTK()
    implicit none
    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') 'Writing reduced number of nodes and elements to vtk file'
    open(unit=1,file='meshReduced.vtk',access='stream',status='replace')
    write(str,'(A)') '# vtk DataFile Version 3.0'//lf; write(1) trim(str)
    write(str,'(A)') '1D Diffusion Mesh'//lf; write(1) trim(str)
    write(str,'(A)') 'ASCII'//lf//lf; write(1) trim(str)
    write(str,'(A)') 'DATASET POLYDATA'//lf; write(1)trim(str)
    write(str,'(A,I0,A)') 'POINTS ',nNodesRed,' DOUBLE'//lf; write(1)trim(str)
    do i=1,nNodesRed
      write(str,'(3F9.2,A)') nodesRed(i,1:3),lf; write(1) trim(str)
    enddo
    write(str,'(A,2(I0,X),A)') lf//'LINES ',nElmtRed,nElmtRed*3,lf; write(1) trim(str)
    do i=1,nElmtRed
      write(str,'(3(I0,X),A)') 2, elmtRedCM(i,1)-1, elmtRedCM(i,2)-1, lf; write(1) trim(str)
    enddo
    write(str,'(A,I0,A)') lf//'CELL_DATA ',nElmtRed,lf; write(1) trim(str)
    write(str,'(A)') 'SCALARS radius FLOAT'//lf; write(1) trim(str)
    write(str,'(A)') 'LOOKUP_TABLE default'//lf; write(1) trim(str)
    do i=1,nElmtRed
      write(str,'(F8.6,A)') erRed(i), lf; write(1) trim(str)
    enddo
    close(1)
    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') achar(13)//dl
    print*, achar(13)//' Created reduced mesh file' 
    print*, ''
  end subroutine

  subroutine writeFinalNodesAndElements2VTK()
    implicit none
    if(outpt .gt. 0) write(*,'(A)',ADVANCE='NO') 'Writing reduced number of nodes and elements to vtk file'
    open(unit=1,file='meshReduced2.vtk',access='stream',status='replace')
    write(str,'(A)') '# vtk DataFile Version 3.0'//lf; write(1) trim(str)
    write(str,'(A)') '1D Diffusion Mesh'//lf; write(1) trim(str)
    write(str,'(A)') 'ASCII'//lf//lf; write(1) trim(str)
    write(str,'(A)') 'DATASET POLYDATA'//lf; write(1)trim(str)
    write(str,'(A,I0,A)') 'POINTS ',nNodesFinal,' DOUBLE'//lf; write(1)trim(str)
    do i=1,nNodesFinal
      write(str,'(3F9.2,A)') nodes(i,1:3),lf; write(1) trim(str)
    enddo
    write(str,'(A,2(I0,X),A)') lf//'LINES ',nElmtFinal,nElmtFinal*3,lf; write(1) trim(str)
    do i=1,nElmtFinal
      write(str,'(3(I0,X),A)') 2, elmtCM(i,1)-1, elmtCM(i,2)-1, lf; write(1) trim(str)
    enddo
    write(str,'(A,I0,A)') lf//'CELL_DATA ',nElmtFinal,lf; write(1) trim(str)
    write(str,'(A)') 'SCALARS radius FLOAT'//lf; write(1) trim(str)
    write(str,'(A)') 'LOOKUP_TABLE default'//lf; write(1) trim(str)
    do i=1,nElmtFinal
      write(str,'(F15.7,A)') er(i), lf; write(1) trim(str)
    enddo
    close(1)
    write(*,'(A)',ADVANCE='NO') achar(13)//dl
    print*, achar(13)//' Created reduced mesh file' 
    print*, ''
  end subroutine
    
end program