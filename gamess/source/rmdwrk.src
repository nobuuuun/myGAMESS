C 18 Apr 16 - SRP - tweaks for RMD
! 15 May 14 - OS & SRP - Created source file for FMO/RMD multistate method
!***********************************************************************
!*MODULE RMDWRK  *DECK RMD_BUILD
!> @brief      Routine does main work for RMD method
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details    Routine builds and diagonalizes composite Hamiltonian
!>             Builds matrix of state gradients 
!>
!> @date : May 15, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>
!> @date : March 25, 2015 - Spencer R. Pruitt
!>         Added synchronization with DDI_BCAST (in other routines too)
!>
!> @param nstates : number of fragmentation states
!>
!> @param coord1d : 1D array of coordinates
!>
!> @param natoms : number of atoms
!>
!> @param atom_name : array of atom names
!>
!> @param fmo_energies : fmo energies
!>
!> @param fmo_gradients :  fmo gradients
!>
!> @param E : energies
!>
!> @param Gradients :  gradients
!***********************************************************************
      subroutine RMD_BUILD(nstates,coord1d,natoms,atom_name,
     &                     fmo_energies,fmo_gradients,E,Gradients)
       
      implicit none
       
      integer :: i
      integer :: ij
      integer :: nstates
      integer :: natoms
      integer :: IERR

      integer :: atom_name(natoms)
      integer :: iwrk(nstates)
      integer :: LEN_H1D
      
      integer, allocatable :: element(:)

      double precision, parameter :: ZERO=0.0D+00

      double precision :: E(1)
      double precision :: fmo_energies(30)
      double precision :: fmo_gradients(30*natoms*3)
      double precision :: Gradients(3*natoms)
      double precision :: SCR(nstates,8)
      double precision :: coord1d(3*natoms)
      double precision :: Htmp(nstates,nstates)
      double precision :: HXtmp(nstates,nstates,3*natoms)
      double precision :: VECtmp(nstates)

      double precision, allocatable :: coord(:)
      double precision, allocatable :: EIG(:)
      double precision, allocatable :: VEC(:,:)
      double precision, allocatable :: HX(:,:,:)
      double precision, allocatable :: H(:,:)
      double precision, allocatable :: H1D(:)

      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      !Length of H1D upper triangular matrix
      LEN_H1D = ((nstates+1)*nstates)/2

      if(.not.allocated(element)) allocate(element(1:natoms))
      if(.not.allocated(coord)) allocate(coord(1:natoms*3))
      if(.not.allocated(VEC)) allocate(VEC(nstates,nstates))
      if(.not.allocated(EIG)) allocate(EIG(nstates))
      if(.not.allocated(H)) allocate(H(nstates,nstates))
      if(.not.allocated(H1D)) allocate(H1D(LEN_H1D))
      if(.not.allocated(HX)) allocate(HX(nstates,nstates,3*natoms))

      !Return values to MD code
      Gradients(:)=ZERO
      E=ZERO
      
      !Passed arrays
      element(:) = atom_name(:)
      coord(:)   = coord1d(:)

      !Temporary arrays (not allocated)
      VECtmp(:)   = ZERO
      Htmp(:,:)   = ZERO
      HXtmp(:,:,:)= ZERO

      !Allocated arrays
      VEC(:,:) = ZERO
      EIG(:)   = ZERO
      H(:,:)   = ZERO
      H1D(:)   = ZERO
      HX(:,:,:)= ZERO

      do i = 1 ,nstates
         H(i,i) = fmo_energies(i)
      enddo

      !Compute couplings and repulsions between state energies
      call ComputeRepulsions(natoms,nstates,coord,element,Htmp)
      call ComputeCouplings(natoms,nstates,coord,element,Htmp)

      H(:,:) = H(:,:) + Htmp(:,:)

      !Transform H to a triangular matrix for GLDIAG 
      CALL CPYSQT(H,H1D,NSTATES,1)

      !CALL DDI_BCAST(1110,'F',H1D,LEN_H1D,MASTER)

      !Diagonalize H
      CALL GLDIAG(NSTATES,NSTATES,NSTATES,H1D,SCR,EIG,VEC,IERR,IWRK)

      !CALL DDI_BCAST(1109,'F',EIG,NSTATES,MASTER)

      E=EIG(1) !Pass the lowest ("ground state") energy back to the dynamics

      VECtmp(:)=VEC(:,1)

      do I=1,nstates
         do IJ=1,natoms*3
            HX(I,I,IJ)=fmo_gradients(IJ+natoms*3*(I-1))
         enddo
      enddo

      Hxtmp(:,:,:) = HX(:,:,:)

      !Compute coupling and repulsions between state gradients
      call computerepulsionsx(natoms,nstates,coord,element,HXtmp)
      call computecouplingsx(natoms,nstates,coord,element,HXtmp)

      HX(:,:,:) = HXtmp(:,:,:)

      !CALL DDI_BCAST(1108,'F',HX,nstates*nstates*3*natoms,MASTER)

      call ComputeHFGrad(nstates,natoms,Gradients,VECtmp,HX)

      !CALL DDI_BCAST(1107,'F',Gradients,natoms*3,MASTER)
      !CALL DDI_BCAST(1107,'F',E,1,MASTER)

      deallocate(element)
      deallocate(coord)
      deallocate(VEC)
      deallocate(EIG)
      deallocate(H)
      deallocate(H1D)
      deallocate(HX)
      
      return
      end
!***********************************************************************
!*MODULE RMDWRK  *DECK computecouplings
!> @brief      Compute the coupling energies
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details    
!>            
!> @date : May 18, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>      
!> @param natoms : number of atoms
!>
!> @param nstates : number of fragmentation states
!>
!> @param coord : coordinates
!>
!> @param element : element
!>
!> @param Htmp : temporary hamiltonian
!***********************************************************************
      subroutine computecouplings(natoms,nstates,coord,
     &                            element,Htmp)
      
      implicit none
      
      !Passed variables
      integer :: natoms
      integer :: nstates

      integer :: element(natoms)

      double precision :: coord(natoms*3)
      double precision :: Htmp(nstates,nstates)
            
      !Temporaray variables
      integer :: i
      integer :: j
      double precision :: ecouple
      
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      do i = 1, nstates
         do j = i+1, nstates 
            call computecoupling(natoms,coord,element,i,j,ecouple)
            Htmp(i,j) = ecouple
            Htmp(j,i) = ecouple
         enddo
      enddo

      !CALL DDI_BCAST(3006,'F',Htmp,nstates*nstates,MASTER)

      return
      end      
!***********************************************************************
!*MODULE RMDWRK  *DECK computecoupling
!> @brief      Compute the coupling energy b/w states I and J 
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details    
!>            
!> @date : May 18, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>            
!> @param natoms : number of atoms
!>
!> @param nstates : number of fragmentation states
!>
!> @param coord : coordinates
!>
!> @param element : element
!>
!> @param istate : state i   
!>
!> @param jstate : state j
!>
!> @param ecouple : coupling between states
!***********************************************************************
      subroutine computecoupling(natoms,coord,element,
     &                           istate,jstate,ecouple)

      use rmddat

      implicit none
      
      !Parameters

      double precision, parameter :: cut_OH  = 4.5353D+00 ! Distance cutoff in state search between O and H atoms
      double precision, parameter :: cut_OH2 = cut_OH * cut_OH !
      double precision, parameter :: AA = 0.01210450810D+00 ! Geometric mean in Bohr
      
      !Passed variables
      integer :: natoms
      integer :: istate
      integer :: jstate

      integer :: element(natoms)

      double precision :: ecouple

      double precision :: coord(natoms*3)
      
      !Temporary variables
      integer :: i
      integer :: j
      integer :: shared
      integer :: iproton

      double precision :: dx
      double precision :: dy
      double precision :: dz
      double precision :: dd
      double precision :: tmp
      double precision :: erepItoJ
      double precision :: erepJtoI
      
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      ecouple = 0.0D+00
      
      !If not shared atom between state I and J, then coupling is zero
      shared = 0
      iproton = -1 !index of shuttling proton between states I and J
      do i = 1, natoms
         if ( element(i)==1 .and. reactive((istate-1)*natoms+i)==1
     &   .and. reactive((jstate-1)*natoms+i)==1 ) then
          !Reactive fragment in state I shares at least one H atom with reactive fragment in state J
            shared = 1
            iproton = i
            exit
         endif
      enddo
      
      if (shared==1) then
         !Must check that I and J reactive fragments are not on the same oxygen,
         !as it is not possible for self-coupling
         !This can occur with bifurcated waters
         do i = 1, natoms
            if ( element(i)==8 .and. reactive((istate-1)*natoms+i)==1
     &      .and. reactive((jstate-1)*natoms+i)==1 ) then
               !Same oxygen. No coupling allowed.
               shared = 0
               exit
           endif
         enddo
      endif
      
      !If passed above, now actually compute the coupling
      if (shared==1) then
         !Compute interstate coupling repulsion for I to J
         erepItoJ = 0.0D+00
         do i = 1, natoms
            if ( element(i)==1 .and. reactive((istate-1)*natoms+i)==1
     &      .and. iproton==i ) then ! reactive H in state I
               do j = 1, natoms
                  if ( element(j)==8.and.
     &            reactive((jstate-1)*natoms+j)==1)then ! reactive O in state J

                  dx = coord(3*(i-1)+1) - coord(3*(j-1)+1)
                  dy = coord(3*(i-1)+2) - coord(3*(j-1)+2)
                  dz = coord(3*(i-1)+3) - coord(3*(j-1)+3)

                  dd =  dx*dx + dy*dy + dz*dz

                     if ( dd < cut_OH2 ) then
                        call Sr(dd, 0.0D+00, cut_OH2, tmp)
                        erepItoJ = erepItoJ + tmp*tmp
                     endif

                  endif
               enddo
            endif
         enddo
        !Compute interstate coupling repulsion for J to I
         erepJtoI = 0.0D+00
         do i = 1, natoms
            if ( element(i)==1.and.
     &      reactive((jstate-1)*natoms+i)==1.and.iproton==i)then ! reactive H in state J
               do j = 1, natoms
                  if ( element(j)==8.and.
     &            reactive((istate-1)*natoms+j)==1)then ! reactive O in state I

                     dx = coord(3*(i-1)+1) - coord(3*(j-1)+1)
                     dy = coord(3*(i-1)+2) - coord(3*(j-1)+2)
                     dz = coord(3*(i-1)+3) - coord(3*(j-1)+3)
 
                     dd =  dx*dx + dy*dy + dz*dz

                     if ( dd < cut_OH2 ) then
                        call Sr(dd, 0.0D+00, cut_OH2, tmp)
                        erepJtoI = erepJtoI + tmp*tmp
                     endif

                  endif
               enddo
            endif
         enddo
      
         !Geometric mean
         ecouple = AA * sqrt( erepItoJ  * erepJtoI )

      endif

      !CALL DDI_BCAST(3007,'F',ecouple,1,MASTER)

      return
      end
!***********************************************************************
!*MODULE RMDWRK  *DECK computecouplingsx
!> @brief      Compute the coupling gradients
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details    
!>            
!> @date : May 18, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>            
!> @param natoms : number of atoms
!>
!> @param nstates : number of fragmentation states
!>
!> @param coord : coordinates
!>
!> @param element : element
!>
!> @param HXtmp : temporary hamiltonian
!***********************************************************************
      subroutine computecouplingsx(natoms,nstates,coord,element,
     &                             HXtmp)
 
      use rmddat
     
      implicit none
      
      !Passed variables
      integer :: natoms
      integer :: nstates

      integer :: element(natoms)

      double precision :: coord(natoms*3)
      double precision :: HXtmp(nstates,nstates,3*natoms)
            
      !Temporaray variables
      integer :: i
      integer :: j
      
      do i = 1, nstates
         do j = i+1, nstates
            call computecouplingX(nstates,natoms,coord,
     &                            element,i,j,HXtmp)
         enddo
      enddo

      return  
      end
!***********************************************************************
!*MODULE RMDWRK  *DECK computecouplingX
!> @brief      Compute the coupling energy b/w states I and J
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details    
!>            
!> @date : May 18, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>      
!> @param natoms : number of atoms
!>
!> @param nstates : number of fragmentation states
!>
!> @param coord : coordinates
!>
!> @param element : element
!>
!> @param istate : state i   
!>
!> @param jstate : state j
!>
!> @param HXtmp : temporary hamiltonian
!***********************************************************************
      subroutine computecouplingX(nstates,natoms,coord,element,
     &                            istate,jstate,HXtmp)

      use rmddat      

      implicit none
      
      !Parameters
      double precision, parameter :: cut_OH  = 4.5353D+00 ! Distance cutoff in state search between O and H atoms
      double precision, parameter :: cut_OH2 = cut_OH * cut_OH !
      double precision, parameter :: AA = 0.01210450810D+00 ! Geometric mean in Bohr
      
      !Passed variables
      integer :: natoms
      integer :: nstates
      integer :: istate
      integer :: jstate

      integer :: element(natoms)

      double precision :: coord(natoms*3)
      double precision :: HXtmp(nstates,nstates,3*natoms)
      
      !Temporary variables
      integer :: i
      integer :: j
      integer :: shared
      integer :: iproton

      double precision :: erepItoJ
      double precision :: erepJtoI
      double precision :: d
      double precision :: dx
      double precision :: dy
      double precision :: dz
      double precision :: dd
      double precision :: fx
      double precision :: fy
      double precision :: fz
      double precision :: tmp
      double precision :: tmp1
      double precision :: tmp2
      double precision :: prefactor

      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      
      !If not shared atom between state I and J, then coupling is zero
      shared = 0
      iproton = -1 !index of shuttling proton between states I and J
      do i = 1, natoms
         if ( element(i)==1.and.reactive((istate-1)*natoms+i)==1
     &   .and.reactive((jstate-1)*natoms+i)==1) then
            !Reactive fragment in I shares at least one H atom with reactive fragment in J
            shared = 1
            iproton = i
            exit
         endif
      enddo
      
      if (shared==1) then
         !Must check that state I and state J reactive fragments are not on the same oxygen
         !as it is not possible for self-coupling
         !This can occur with bifurcated waters
         do i = 1, natoms
            if ( element(i)==8.and.reactive((istate-1)*natoms+i)==1
     &      .and. reactive((jstate-1)*natoms+i)==1) then
               !Same oxygen. No coupling allowed.
               shared = 0
               exit
            endif
         enddo
      endif
      
      ! ** If passed above, now actually compute the coupling **
      if (shared==1) then
         !Compute interstate coupling repulsion
         erepItoJ = 0.0D+00
         do i = 1, natoms
            if ( element(i)==1.and.reactive((istate-1)*natoms+i)==1
     &      .and.iproton==i) then ! reactive H in state I
               do j = 1, natoms
                  if ( element(j)==8.and.
     &            reactive((jstate-1)*natoms+j)==1) then ! reactive O in state J

                     dx = coord(3*(i-1)+1) - coord(3*(j-1)+1)
                     dy = coord(3*(i-1)+2) - coord(3*(j-1)+2)
                     dz = coord(3*(i-1)+3) - coord(3*(j-1)+3)

                     dd =  dx*dx + dy*dy + dz*dz

                     if ( dd < cut_OH2 ) then
                        call Sr(dd, 0.0D+00, cut_OH2, tmp)
                        erepItoJ = erepItoJ + tmp*tmp
                     endif

                  endif
               enddo
            endif
         enddo
      
         erepJtoI = 0.0D+00
         do i = 1, natoms
            if ( element(i)==1.and.reactive((jstate-1)*natoms+i)==1
     &      .and. iproton==i ) then ! reactive H in state J
               do j = 1, natoms
                  if ( element(j)==8.and.
     &            reactive((istate-1)*natoms+j)==1  ) then ! reactive O in state I

                     dx = coord(3*(i-1)+1) - coord(3*(j-1)+1)
                     dy = coord(3*(i-1)+2) - coord(3*(j-1)+2)
                     dz = coord(3*(i-1)+3) - coord(3*(j-1)+3)

                     dd =  dx*dx + dy*dy + dz*dz
                     if ( dd < cut_OH2 ) then
                        call Sr(dd, 0.0D+00, cut_OH2, tmp)
                        erepJtoI = erepJtoI + tmp*tmp
                     endif

                  endif
               enddo
            endif
         enddo

         prefactor = AA / ( 2.0D+00 * sqrt (erepItoJ * erepJtoI) )
      
        !Now get the gradients
         do i = 1, natoms
            if ( element(i)==1.and.reactive((istate-1)*natoms+i)==1
     &      .and. iproton==i ) then ! reactive H in state I
               do j = 1, natoms
                  if ( element(j)==8.and.
     &            reactive((jstate-1)*natoms+j)==1 ) then ! reactive O in state J

                     dx = coord(3*(i-1)+1) - coord(3*(j-1)+1)
                     dy = coord(3*(i-1)+2) - coord(3*(j-1)+2)
                     dz = coord(3*(i-1)+3) - coord(3*(j-1)+3)

                     dd =  dx*dx + dy*dy + dz*dz

                     if ( dd < cut_OH2 ) then
                        d = sqrt(dd)
                        call Sr(dd, 0.0D+00, cut_OH2, tmp1)
                        call dSrdr(dd, 0.0D+00, cut_OH2, tmp2)
                        tmp = 2.0D+00 * tmp1 * tmp2
                        tmp = tmp * prefactor * erepJtoI
      
                        fx = tmp * dx / d
                        fy = tmp * dy / d
                        fz = tmp * dz / d
      
      HXtmp(istate,jstate,3*(i-1)+1) = HXtmp(istate,jstate,3*(i-1)+1)+fx
      HXtmp(istate,jstate,3*(i-1)+2) = HXtmp(istate,jstate,3*(i-1)+2)+fy
      HXtmp(istate,jstate,3*(i-1)+3) = HXtmp(istate,jstate,3*(i-1)+3)+fz

      HXtmp(istate,jstate,3*(j-1)+1) = HXtmp(istate,jstate,3*(j-1)+1)-fx
      HXtmp(istate,jstate,3*(j-1)+2) = HXtmp(istate,jstate,3*(j-1)+2)-fy
      HXtmp(istate,jstate,3*(j-1)+3) = HXtmp(istate,jstate,3*(j-1)+3)-fz

                     endif
                  endif
               enddo
            endif
         enddo

         do i = 1, natoms
            if ( element(i)==1.and.reactive((jstate-1)*natoms+i)==1
     &      .and. iproton==i ) then ! reactive H in state J
               do j = 1, natoms
                  if ( element(j)==8.and.
     &            reactive((istate-1)*natoms+j)==1 ) then ! reactive O in state I

                     dx = coord(3*(i-1)+1) - coord(3*(j-1)+1)
                     dy = coord(3*(i-1)+2) - coord(3*(j-1)+2)
                     dz = coord(3*(i-1)+3) - coord(3*(j-1)+3)

                     dd =  dx*dx + dy*dy + dz*dz

                     if ( dd < cut_OH2 ) then
                        d = sqrt(dd)
                        call Sr(dd, 0.0D+00, cut_OH2, tmp1)
                        call dSrdr(dd, 0.0D+00, cut_OH2, tmp2)
                        tmp = 2.0D+00 * tmp1 * tmp2
                        tmp = tmp * prefactor * erepItoJ
      
                        fx = tmp * dx / d
                        fy = tmp * dy / d
                        fz = tmp * dz / d
                        
      HXtmp(istate,jstate,3*(i-1)+1) = HXtmp(istate,jstate,3*(i-1)+1)+fx
      HXtmp(istate,jstate,3*(i-1)+2) = HXtmp(istate,jstate,3*(i-1)+2)+fy
      HXtmp(istate,jstate,3*(i-1)+3) = HXtmp(istate,jstate,3*(i-1)+3)+fz
      
      HXtmp(istate,jstate,3*(j-1)+1) = HXtmp(istate,jstate,3*(j-1)+1)-fx
      HXtmp(istate,jstate,3*(j-1)+2) = HXtmp(istate,jstate,3*(j-1)+2)-fy
      HXtmp(istate,jstate,3*(j-1)+3) = HXtmp(istate,jstate,3*(j-1)+3)-fz
                     endif
                  endif
               enddo
            endif
         enddo
      
         !Symmetrize coupling gradient for convenience
         do i = 1, natoms
      HXtmp(jstate,istate,3*(i-1)+1) = HXtmp(istate,jstate,3*(i-1)+1)
      HXtmp(jstate,istate,3*(i-1)+2) = HXtmp(istate,jstate,3*(i-1)+2)
      HXtmp(jstate,istate,3*(i-1)+3) = HXtmp(istate,jstate,3*(i-1)+3)
         enddo

      endif

      !CALL DDI_BCAST(1005,'F',Hxtmp,nstates*nstates*3*natoms,MASTER)      

      return      
      end

!***********************************************************************
!*MODULE RMDWRK  *DECK computerepulsions
!> @brief      Compute the repulsions energies
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details
!>
!> @date : June 9, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>
!> @param natoms : number of atoms
!>
!> @param nstates : number of fragmentation states
!>
!> @param coord : coordinates
!>
!> @param element : element
!>
!> @param Htmp : temporary hamiltonian
!***********************************************************************
      subroutine computerepulsions(natoms,nstates,coord,
     &                            element,Htmp)

      use rmddat

      implicit none

      double precision, parameter :: ZERO=0.0D+00

      !Passed variables
      integer :: natoms
      integer :: nstates

      integer :: element(natoms)

      double precision :: coord(natoms*3)
      double precision :: Htmp(nstates,nstates)

      !Temporaray variables
      integer :: i
      double precision :: erepulsion, epenalty

      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      do i = 1, nstates
         epenalty = ZERO
         erepulsion = ZERO
         call computerepulsion(natoms,coord,element,i,erepulsion)
         if(penlty) call computepenalty(natoms,coord,element,i,epenalty)
         Htmp(i,i) = Htmp(i,i) + erepulsion + epenalty
      enddo

      !CALL DDI_BCAST(2005,'F',Htmp,nstates*nstates,MASTER)

      return
      end

!***********************************************************************
!*MODULE RMDWRK  *DECK computerepulsion
!> @brief      Compute the repulsion energy of states I
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details
!>
!> @date : June 9, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>
!> @param natoms : number of atoms
!>
!> @param coord : coordinates
!>
!> @param element : element
!>
!> @param istate : state i   
!>
!> @param erepulsion : repulsion energy
!***********************************************************************
      subroutine computerepulsion(natoms,coord,element,
     &                           istate,erepulsion)

      use rmddat

      implicit none

      !Parameters
      double precision, parameter :: ZERO=0.0D+00
      double precision, parameter :: cut_OH  = 4.5353D+00 ! Distance cutoff in state search between O and H atoms
      double precision, parameter :: cut_OH2 = cut_OH * cut_OH !
      double precision, parameter :: AA = 0.01210450810D+00 ! Geometric mean in Bohr

      !Passed variables
      integer :: natoms
      integer :: istate

      integer :: element(natoms)

      double precision :: erepulsion

      double precision :: coord(natoms*3)

      !Temporary variables
      integer :: i
      integer :: j

      double precision :: dx
      double precision :: dy
      double precision :: dz
      double precision :: dd
      double precision :: tmp

      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      erepulsion = ZERO

      !Let each hydronium H atom interact with all non-hydronium O atoms in state istate
      do i = 1, natoms
        if ( element(i)==1 .and. reactive((istate-1)*natoms+i)==1 ) then
          do j = 1, natoms
            if ( element(j)==8 .and.
     &      reactive((istate-1)*natoms+j)==0 ) then

              dx = coord(3*(i-1)+1) - coord(3*(j-1)+1)
              dy = coord(3*(i-1)+2) - coord(3*(j-1)+2)
              dz = coord(3*(i-1)+3) - coord(3*(j-1)+3)

              dd =  dx*dx + dy*dy + dz*dz

              if ( dd < cut_OH2 ) then
                call Sr(dd, ZERO, cut_OH2, tmp)
                erepulsion = erepulsion + tmp*tmp
              endif
            endif
          enddo
        endif
      enddo

      erepulsion = erepulsion * AA

      !CALL DDI_BCAST(2006,'F',erepulsion,1,MASTER)

      return
      end
!***********************************************************************
!*MODULE RMDWRK  *DECK computepenalty
!> @brief      Compute the internal penalty energy for state I
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2015
!>
!> @details
!>
!> @date : February 6, 2015 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>
!> @param natoms : number of atoms
!>
!> @param coord : coordinates
!>
!> @param element : element
!>
!> @param istate : state i   
!>
!> @param epenalty : energy penalty term
!***********************************************************************
      subroutine computepenalty (natoms,coord,element,istate,epenalty)

      use rmddat

      implicit none
       
      !Parameters
      double precision, parameter :: ZERO=0.0D+00
      double precision, parameter :: ONE=1.0D+00
      double precision, parameter :: cut_OH  = 4.5353D+00 ! Distance cutoff in state search between O and H atoms
      double precision, parameter :: cut_OH2 = cut_OH * cut_OH !
      double precision, parameter :: BB = 0.22D+00 ! hartree / Bohr**2
      double precision, parameter :: r0 = 2.65D+00
      double precision, parameter :: r02 = r0 * r0
       
      !Passed variables
      integer :: natoms
      integer :: istate
      integer :: element(natoms)
      double precision :: coord(natoms*3)
      double precision :: epenalty
       
      !Temporary variables
      double precision :: dx
      double precision :: dy
      double precision :: dz
      double precision :: dd
      double precision :: tmp
      double precision :: eint
      integer :: i
      integer :: j

      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      
      eint = ZERO
      epenalty = ZERO
      
      do i = 1, natoms
         if ( element(i)==1 ) then
            do j = 1, natoms
               if ( element(j)==8 ) then
                  if ( irmddat((istate-1)*natoms+i) ==
     &                 irmddat((istate-1)*natoms+j) )then
                     dx = coord(3*(i-1)+1) - coord(3*(j-1)+1)
                     dy = coord(3*(i-1)+2) - coord(3*(j-1)+2)
                     dz = coord(3*(i-1)+3) - coord(3*(j-1)+3)
                     dd =  dx*dx + dy*dy + dz*dz
                     if ( dd > r02 ) then
                        call Sr(dd, r02, cut_OH2, tmp)
                        eint = eint + ( ONE - tmp ) * ( ONE - tmp )
                     endif
                  endif
               endif
            enddo
         endif
      enddo
       
      epenalty = eint * BB

      !CALL DDI_BCAST(2007,'F',epenalty,1,MASTER)

      IF(MASWRK) WRITE(6,9001) ISTATE, epenalty

 9001 FORMAT(2X,' '/
     *       2X,49(1H-),/,
     *       2X,'Penalty (Hartrees) for state #',I2,' =',F15.10,/
     *       2X,49(1H-))
      return
      end
!***********************************************************************
!*MODULE RMDWRK  *DECK computerepulsionsx
!> @brief      Compute the repulsion gradients
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details
!>
!> @date : June 9, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>
!> @param natoms : number of atoms
!>
!> @param nstates : number of fragmentation states
!>
!> @param coord : coordinates
!>
!> @param element : element
!>
!> @param HXtmp : temporary hamiltonian
!***********************************************************************
      subroutine computerepulsionsx(natoms,nstates,
     &                              coord,element,HXtmp)

      use rmddat

      implicit none

      !Passed variables
      integer :: natoms
      integer :: nstates

      integer :: element(natoms)

      double precision :: coord(natoms*3)
      double precision :: HXtmp(nstates,nstates,3*natoms)

      !Temporaray variables
      integer :: i
      
      do i = 1, nstates
        call computerepulsionX(natoms,nstates,coord,element,i,HXtmp)
        if(penlty) call computepenaltyX (natoms,nstates,coord,
     &                                   element,i,HXtmp)
      enddo

      return
      end


!***********************************************************************
!*MODULE RMDWRK  *DECK computerepulsionX
!> @brief      Compute the repulsion gradients for state I
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details
!>
!> @date : June 9, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>
!> @param natoms : number of atoms
!>
!> @param nstates : number of fragmentation states
!>
!> @param coord : coordinates
!>
!> @param element : element
!>
!> @param HXtmp : temporary hamiltonian
!***********************************************************************
      subroutine computerepulsionX(natoms,nstates,coord,element,
     &                            istate,HXtmp)

      use rmddat

      implicit none

      !Parameters
      double precision, parameter :: ZERO=0.0D+00
      double precision, parameter :: cut_OH  = 4.5353D+00 ! Distance cutoff in state search between O and H atoms
      double precision, parameter :: cut_OH2 = cut_OH * cut_OH !
      double precision, parameter :: AA = 0.01210450810D+00 ! Geometric mean in Bohr

      !Passed variables
      integer :: natoms
      integer :: istate
      integer :: nstates

      integer :: element(natoms)

      double precision :: coord(natoms*3)
      double precision :: HXtmp(nstates,nstates,3*natoms)

      !Temporary variables
      integer :: i
      integer :: j

      double precision :: d
      double precision :: dx
      double precision :: dy
      double precision :: dz
      double precision :: dd
      double precision :: fx
      double precision :: fy
      double precision :: fz
      double precision :: tmp
      double precision :: tmp1
      double precision :: tmp2

      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      !Let each hydronium H atom interact with all non-hydronium O atoms in state istate
      do i = 1, natoms
        if ( element(i)==1 .and. reactive((istate-1)*natoms+i)==1 ) then
          do j = 1, natoms
            if ( element(j)==8 .and.
     &      reactive((istate-1)*natoms+j)==0 ) then

              dx = coord(3*(i-1)+1) - coord(3*(j-1)+1)
              dy = coord(3*(i-1)+2) - coord(3*(j-1)+2)
              dz = coord(3*(i-1)+3) - coord(3*(j-1)+3)

              dd =  dx*dx + dy*dy + dz*dz

              if ( dd < cut_OH2 ) then
                d = sqrt(dd)
                call Sr(dd, ZERO, cut_OH2, tmp1)
                call dSrdr(dd, ZERO, cut_OH2, tmp2)
                tmp = 2.0D+00 * tmp1 * tmp2
                tmp = tmp * AA
                !tmp = tmp * b2a

                fx = tmp * dx / d
                fy = tmp * dy / d
                fz = tmp * dz / d

      HXtmp(istate,istate,3*(i-1)+1)=HXtmp(istate,istate,3*(i-1)+1)+fx
      HXtmp(istate,istate,3*(i-1)+2)=HXtmp(istate,istate,3*(i-1)+2)+fy
      HXtmp(istate,istate,3*(i-1)+3)=HXtmp(istate,istate,3*(i-1)+3)+fz

      HXtmp(istate,istate,3*(j-1)+1)=HXtmp(istate,istate,3*(j-1)+1)-fx
      HXtmp(istate,istate,3*(j-1)+2)=HXtmp(istate,istate,3*(j-1)+2)-fy
      HXtmp(istate,istate,3*(j-1)+3)=HXtmp(istate,istate,3*(j-1)+3)-fz
              endif
            endif
          enddo
        endif
      enddo

      !CALL DDI_BCAST(1006,'F',Hxtmp,nstates*nstates*3*natoms,MASTER)

      return
      end
!***********************************************************************
!*MODULE RMDWRK  *DECK computepenaltyX
!> @brief      Compute the internal penalty gradient for state I
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2015
!>
!> @details
!>
!> @date : February 6, 2015 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>
!> @param natoms : number of atoms
!>
!> @param nstates : number of fragmentation states
!>
!> @param coord : coordinates
!>
!> @param element : element
!>
!> @param istate : state i   
!>
!> @param HXtmp : temporary hamiltonian
!***********************************************************************
      subroutine computepenaltyX
     &           (natoms,nstates,coord,element,istate,HXtmp)
      
      use rmddat

      implicit none
      
      !Parameters
      double precision, parameter :: ONE=1.0D+00
      double precision, parameter :: TWO=2.0D+00
      double precision, parameter :: cut_OH  = 4.5353D+00 ! Distance cutoff in state search between O and H atoms
      double precision, parameter :: cut_OH2 = cut_OH * cut_OH !
      double precision, parameter :: BB = 0.22D+00 ! hartree / Bohr**2
      double precision, parameter :: r0 = 2.65D+00
      double precision, parameter :: r02 = r0 * r0
      
      !Passed variables
      integer :: natoms
      integer :: nstates
      integer :: istate
      integer :: element(natoms)
      double precision :: coord(natoms*3)
      double precision :: HXtmp(nstates,nstates,3*natoms)
      
      !Temporary variables
      double precision :: d
      double precision :: dx
      double precision :: dy
      double precision :: dz
      double precision :: dd
      double precision :: fx
      double precision :: fy
      double precision :: fz
      double precision :: tmp
      double precision :: tmp1
      double precision :: tmp2
      integer :: i
      integer :: j
      
      do i = 1, natoms
         if ( element(i)==1 ) then
            do j = 1, natoms
               if ( element(j)==8 ) then
                  if ( irmddat((istate-1)*natoms+i) ==
     &                 irmddat((istate-1)*natoms+j) ) then
                     dx = coord(3*(i-1)+1) - coord(3*(j-1)+1)
                     dy = coord(3*(i-1)+2) - coord(3*(j-1)+2)
                     dz = coord(3*(i-1)+3) - coord(3*(j-1)+3)
                     dd =  dx*dx + dy*dy + dz*dz
      
                     if ( dd > r02 ) then
                        d = sqrt(dd)
                        call Sr(dd, r02, cut_OH2, tmp1)
                        call dSrdr(dd, r02, cut_OH2, tmp2)
                        tmp = TWO*(ONE-tmp1)*(-ONE*tmp2)
                        tmp = BB * tmp
      
      
                  !Convert from hartree/Angs to hartree/bohr
                        !tmp = tmp * b2a
      
                        fx = tmp * dx / d
                        fy = tmp * dy / d
                        fz = tmp * dz / d
                        
        HXtmp(istate,istate,3*(i-1)+1) = 
     $  (HXtmp(istate,istate,3*(i-1)+1) + fx)
        HXtmp(istate,istate,3*(i-1)+2) = 
     $  (HXtmp(istate,istate,3*(i-1)+2) + fy)
        HXtmp(istate,istate,3*(i-1)+3) = 
     $  (HXtmp(istate,istate,3*(i-1)+3) + fz)
                        
        HXtmp(istate,istate,3*(j-1)+1) = 
     $ (HXtmp(istate,istate,3*(j-1)+1) - fx)
        HXtmp(istate,istate,3*(j-1)+2) = 
     $ (HXtmp(istate,istate,3*(j-1)+2) - fy)
        HXtmp(istate,istate,3*(j-1)+3) = 
     $ (HXtmp(istate,istate,3*(j-1)+3) - fz)
      
                     endif
                  endif
               endif
            enddo
         endif
      enddo
      return
      end
!***********************************************************************
!*MODULE RMDWRK  *DECK Sr
!> @brief      Switching function
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details
!>
!> @date : May 18, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>
!> @param r2 : r squared
!>
!> @param rin2 : input r squared
!>
!> @param rout2 : output r squared
!>
!> @param output : output
!***********************************************************************
      subroutine Sr ( r2, rin2, rout2, output )
      
      implicit none
      
      double precision :: r2
      double precision :: rin2
      double precision :: rout2
      double precision :: output
      double precision :: denom
      double precision :: tmp

      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      
      if ( r2 < rin2 ) then
         output = 1.0D+00
      else if ( r2 > rout2 ) then
         output = 0.0D+00
      else
         denom = rout2 - rin2
         denom = denom * denom * denom
         tmp = rout2 - r2
         output = tmp*tmp* ( rout2 + 2.0D+00*r2 - 3.0D+00*rin2 ) / denom
      endif

      !CALL DDI_BCAST(2000,'F',output,1,MASTER)

      return
      end
!***********************************************************************
!*MODULE RMDWRK  *DECK dSrdr
!> @brief      Switching function derivative
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details    
!>            
!> @date : May 18, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>
!> @param r2 : r squared
!>
!> @param rin2 : input r squared
!>
!> @param rout2 : output r squared
!>
!> @param output : output
!***********************************************************************
      subroutine dSrdr ( r2, rin2, rout2, output )
      
      implicit none
      
      double precision :: r2
      double precision :: rin2
      double precision :: rout2
      double precision :: output
      double precision :: denom
      double precision :: r

      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      
      if ( r2 < rin2 ) then
         output = 0.0D+00
      else if ( r2 > rout2 ) then
         output = 0.0D+00
      else
         denom = rout2 - rin2
         denom = denom * denom * denom
         r = sqrt( r2 )
         output = -12.0D+00*r * (rout2 - r2) * (r2 -rin2) / denom
      endif

      !CALL DDI_BCAST(2001,'F',output,1,MASTER)

      return
      end
!***********************************************************************
!*MODULE RMDWRK  *DECK ComputeHFGrad
!> @brief      Routine computes composite gradient
!>
!> @author Olaseni Sode and Spencer R. Pruitt
!>         -2014
!>
!> @details    Routine computes composite gradient via Hellmann-Feynman 
!>             theorem
!>            
!> @date : May 15, 2014 - Olaseni Sode & Spencer R. Pruitt
!>         Created subroutine
!>
!> @param nstates : number of fragmentation states
!>
!> @param natoms : number of atoms
!>
!> @param GSGradients : final gradients
!>
!> @param GSCoeffs : final coefficients
!>
!> @param HX : final gradient hamiltonian
!***********************************************************************
      subroutine ComputeHFGrad(nstates,natoms,
     &           GSGradients,GSCoeffs,HX)
      
      implicit none

      !Passed variables      
      integer :: nstates
      integer :: natoms

      !Temporary variables      
      integer :: i
      integer :: istate
      integer :: jstate

      double precision :: cI
      double precision :: cJ
      double precision :: cc

      double precision :: HX(nstates,nstates,3*natoms)
      double precision :: GSGradients(natoms*3)
      double precision :: GSCoeffs(nstates)

      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      GSGradients(:)=0.0D+00

      do istate = 1, nstates
         cI = GSCoeffs(istate)
         do jstate = 1, nstates
            cJ = GSCoeffs(jstate)
            cc = cI*cJ
            do i = 1, natoms
               GSGradients(3*(i-1)+1)=GSGradients(3*(i-1)+1)+
     &                                cc*HX(istate,jstate,3*(i-1)+1)

               GSGradients(3*(i-1)+2)=GSGradients(3*(i-1)+2)+
     &                                cc*HX(istate,jstate,3*(i-1)+2)

               GSGradients(3*(i-1)+3)=GSGradients(3*(i-1)+3)+
     &                                cc*HX(istate,jstate,3*(i-1)+3)
            enddo
         enddo
      enddo

      !CALL DDI_BCAST(3000,'F',GSGradients,natoms*3,MASTER)

      return       
      end
