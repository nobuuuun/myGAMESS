      program Main
      implicit none
      integer, parameter :: MEM_REPL=1, MEM_DIST=1
      double precision dummy
      integer h(10),mw(10),mavail
      integer i,me,np
      
      call DDI_Init()
      call DDI_Nproc(np,me)
      write (6,*) 'My rank=',me,' PID=',getpid()
      flush 6
C      call sleep(20)
      
      call DDI_Memory(MEM_DIST,MEM_REPL,dummy)

      if (me.eq.0) call DDI_Memory_map_dbgprint()

      do i=1,6
         mw(i)=100
         call create_and_fill(h(i),mw(i))
      enddo
      call ddi_sync(100)

      do i=1,6
         call check_array(h(i),mw(i))
      enddo
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      call ddi_sync(200)

      if (me.eq.0) then
         write (6,*) 'Arrays are created and verified'
         write (6,'(/,''Test 1: simple out-of-order deallocation'')')
      endif
      flush 6
      call DDI_Destroy(h(5))
      h(5)=-1
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      flush 6

      if (me.eq.0) write (6,'(/,''Test 2: expanding down'')')
      flush 6
      call DDI_Destroy(h(4))
      h(4)=-1
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      flush 6

      if (me.eq.0)
     $     write (6,'(/,''Test 3: Deallocation at the bottom'')')
      flush 6
      call DDI_Destroy(h(1))
      h(1)=-1
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      flush 6

      if (me.eq.0) write (6,'(/,''Test 4: expanding up'')')
      flush 6
      call DDI_Destroy(h(2))
      h(2)=-1
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      flush 6

      if (me.eq.0) write (6,'(/,''Test 5: Coalescing'')')
      flush 6
      call DDI_Destroy(h(3))
      h(3)=-1
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      flush 6

      if (me.eq.0) write (6,'(/,''Test 6: collapsing'')')
      flush 6
      call DDI_Destroy(h(6))
      h(6)=-1
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      flush 6

      call ddi_sync(300)
      if (me.eq.0)
     $     write (6,
     $     '(/,''All arrays are destroyed, starting new test set'')')
      flush 6

      do i=1,3
         call create_and_fill(h(i),mw(i))
      enddo
      call ddi_sync(1100)

      do i=1,3
         call check_array(h(i),mw(i))
      enddo
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      call ddi_sync(1200)

      if (me.eq.0)
     $     write (6,'(/,''Arrays are created and verified'')')
      

      if (me.eq.0) write (6,'(/,''Test 1: exact fill'')')
      flush 6
      call DDI_Destroy(h(2))
      call create_and_fill(h(2),mw(i))
      do i=1,3
         call check_array(h(i),mw(i))
      enddo
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      flush 6

      if (me.eq.0) write (6,'(/,''Test 2: smaller array'')')
      flush 6
      call DDI_Destroy(h(2))
      mw(2)=mw(2)/4
      call create_and_fill(h(2),mw(2))
      do i=1,3
         call check_array(h(i),mw(i))
      enddo
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      flush 6

      if (me.eq.0) write (6,'(/,''Test 3: allocation of all space'')')
      flush 6
      call DDI_Memory_avail(mavail)
      write (6,'(''['',I2,''] available memory='',I12)') me,mavail
c$$$      call ddi_gsumi(2000,mavail,1)
c$$$      if (me.eq.0) write (6,*) 'Total available memory:',mavail
      flush 6
      mw(4)=mavail
      call create_and_fill(h(4),mw(4))
      do i=1,4
         call check_array(h(i),mw(i))
      enddo
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      flush 6

      if (me.eq.0) write (6,'(/,''Test 4: low-mem allocation'')')
      flush 6
      mw(5)=mw(2)
      call create_and_fill(h(5),mw(5))
      do i=1,5
         call check_array(h(i),mw(i))
      enddo
      if (me.eq.0) call DDI_Memory_map_dbgprint()
      flush 6

      if (me.eq.0) write (6,'(/,''Test 6: deallocation'')')
      flush 6
      do i=5,1,-1
         if (me.eq.0) write (6,*) 'Destroying h(',i,')=',h(i)
         if (h(i).ge.0) call DDI_Destroy(h(i))
         if (me.eq.0) call DDI_Memory_map_dbgprint()
      enddo
      flush 6

      call DDI_Finalize()
      end


C>    @brief Create and populate the array. This is a collective operation.
      subroutine create_and_fill(h,mwords)
      implicit none
      integer, intent(out) ::  h
      integer, intent(in) :: mwords
      integer, parameter :: NBUF=1000
      integer me,np
      integer j,jcol,jj
      integer mycol0
      double precision buf(NBUF),fillval

      call DDI_Nproc(np,me)

      call DDI_Create(1,np*mwords,h)
      mycol0=me*mwords
      do jcol=1,mwords,NBUF
         do j=1,min(NBUF,mwords-jcol+1)
            jj=jcol+j-1
            fillval = (mycol0+jj)+1000*h
            buf(j) = fillval
         enddo
         call DDI_put(h, 1,1, mycol0+jcol,
     $        mycol0+jcol+min(NBUF,mwords-jcol), buf)
      enddo
      write (6,'(''['',I2,''] Array '',I2,'' is created and filled.'')')
     $     me,h
      flush 6
      end
      

C>    Check the correctness of the array h of size mwords (across processes)     
      subroutine check_array(h,mwords)
      implicit none
      integer, intent(out) ::  h
      integer, intent(in) :: mwords
      integer, parameter :: NBUF=1000
      integer me,np
      integer j,jcol,jj
      integer mycol0
      double precision buf(NBUF), fillval, v

      call DDI_Nproc(np,me)
      
      if (me.eq.0) then
         mycol0=(np-1)*mwords
      else
         mycol0=(me-1)*mwords
      endif

      mycol0=me
      do jcol=1,mwords,NBUF
         call DDI_get(h, 1,1, mycol0+jcol,
     $        mycol0+jcol+min(NBUF,mwords-jcol), buf)
         do j=1,min(NBUF,mwords-jcol+1)
            jj=jcol+j-1
            fillval=(mycol0+jj)+1000*h
            v=buf(j)
            if (abs(v-fillval) .gt. 1E-4) then
               write (6,1000) me, h, jj, mycol0+jj, fillval, v
 1000          format ('Proc ',I3,' array ',I3,' at local col=',I5,
     $              ' glob col=',I5,' expected ',F10.1,' got ',F10.1)
            endif
         enddo
      enddo
      write (6,'(''['',I2,''] Array '',I2,'' is OK'')') me,h
      flush 6
      end
