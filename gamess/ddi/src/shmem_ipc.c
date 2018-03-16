/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 *
 * Wrapper subroutines for shared memory calls (SysV or POSIX)
 * Wrapper subroutines for semaphore calls     (SysV only)
 *
 * Author: Ryan M. Olson
 *  6 Feb 13 - Jason Rigby - retry if SEMOP fails due to interupt signal
 * 15 Dec 13 - Alex Gaenko - POSIX Shared memory implementation
\* -------------------------------------------------------------------- */

 # include "ddi_base.h"
 # if defined USE_SYSV

/* ------------------------------------------------------- *\
   Remove a shared-memory id from the system (either type)
\* ------------------------------------------------------- */
   void DDI_Shm_remove(int shmid) {
      Shmctl(shmid,IPC_RMID,NULL);
   }


/* ------------------------------- *\
   SYSV interface to shared memory
/* ------------------------------- */

# if defined USE_SHMEM_SYSV

/* --------------------------- *\
   Wrapper function for shmget
\* --------------------------- */
   int Shmget(key_t key, size_t size, int flag) {
     int shmid;
     if((shmid = shmget(key,size,flag)) < 0) {
       fprintf(stdout,"%s: shmget returned an error.\n",DDI_Id());
       switch(errno) {
         case ENOSPC:
           fprintf(stdout," Error ENOSPC: The number of shared memory identifiers for the\nnode has been reached.\nThis may be due to asking for too much shared memory!\n");
           break;
         case ENOMEM:
           fprintf(stdout," Error ENOMEM: Not enough memory to allocate %li bytes of shared memory.\n",(long)size);
           break;
         case EACCES:
           fprintf(stdout," Error EACCES\n");
           break;
         case EEXIST:
           fprintf(stdout," Error EEXIST\n");
           break;
         case EINVAL:
           fprintf(stdout," Error EINVAL: Attempting to create %lu bytes of shared memory.\n",(unsigned long) size);
           fprintf(stdout," Check system limits on the size of SysV shared memory segments.\n");
           fprintf(stdout,"\n");
           fprintf(stdout," The file ~/gamess/ddi/readme.ddi contains information on how to display\n");
           fprintf(stdout," the current SystemV memory settings, and how to increase their sizes.\n");
           fprintf(stdout," Increasing the setting requires the root password, and usually a sytem reboot.\n");
           fprintf(stdout,"\n");
           fflush(stdout);
           break;
         case ENOENT:
           fprintf(stdout," Error ENOENT\n");
           break;
         default:
           fprintf(stdout," unusual shmget errno=%d\n",errno); break;
           break;
       }
       Fatal_error(911);
     }

     return shmid;
   }


/* -------------------------- *\
   Wrapper function for shmat
\* -------------------------- */
   void *Shmat(int shmid, void *addr, int flag) {
     void *shmaddr = NULL;
     void *error   = (void *) -1;
     const DDI_Comm *comm = (DDI_Comm *) &gv(ddi_base_comm);


     if((shmaddr = shmat(shmid,addr,flag)) == error) {
       fprintf(stdout," DDI Process %i: shmat returned an error.\n",comm->me);
       switch(errno) {
         case EINVAL:
           fprintf(stdout," Error EINVAL: shmid=%i is not a valid shared memory identifier.\n",shmid);
           break;
         case ENOMEM:
           fprintf(stdout," Error ENOMEM: Can't attach shared memory segment due to data space limits.\n");
           break;
         case EACCES:
           fprintf(stdout," Error EACCES: Permission to attach to the shared memory segment denied.\n");
           break;
         case EMFILE:
           fprintf(stdout," Error EMFILE: no. of shared memory segments exceeds system-imposed limit.\n");
           break;
         default:
           fprintf(stdout," unusual shmat errno=%d\n",errno); break;
           break;
       }
       fflush(stdout);
       Fatal_error(911);
     }

     return shmaddr;
   }


/* --------------------------- *\
   Wrapper function for shmctl
\* --------------------------- */
   int Shmctl(int shmid, int cmd, struct shmid_ds *buff) {
     int ret;
     const DDI_Comm *comm = (DDI_Comm *) &gv(ddi_base_comm);

     if((ret = shmctl(shmid,cmd,buff)) != 0) {
       fprintf(stdout," DDI Process %i: shmctl return an error.\n",comm->me);
       switch(errno) {
         case EINVAL:
           fprintf(stdout," Error EINVAL: shmid is not a valid shared memory segment.\n");
           break;
         case EFAULT:
           fprintf(stdout," Error EFAULT: argument 3 is not a valid struct shmid_ds.\n");
           break;
         case EPERM:
           fprintf(stdout," Error EPREM: permission to access/change shared mem segment denied.\n");
           break;
         default:
           fprintf(stdout," unusual shmctl errno=%d\n",errno); break;
           break;
       }
       Fatal_error(911);
     }

     return ret;

   }

# endif
     /*   the line above closes the SYSV shared memory implementation   */

/* -------------------------------- *\
   POSIX interface to shared memory
/* -------------------------------- */

# ifdef USE_SHMEM_POSIX

#include <sys/mman.h>
#include <sys/stat.h>

static char* number2shmemname(char* buf, int n)

/*       number2shmemnam:
    @brief Makes a shmem name from a number.
    @param{buf} Buffer to hold the name (must be of sufficient size!)
    @param{n} The number
    @return The pointer to the buffer.
*/
{
  sprintf(buf,"/ddi-shmem-%x",n);
  return buf;
}

/* ------------------------------------------------------ *\
   Wrapper function for get a POSIX shared memory segment
\* ------------------------------------------------------ */
int Shmget(key_t key, size_t sz, int flag)

/*  note that only the 2nd argument is used below  */

{
  char namebuf[NAME_MAX];
  int memid, fd;
  unsigned int i;
  const unsigned int ntries=1000000; /* no. tries to generate a unique name */

  for (i=ntries; i!=0; i--) {
    memid=rand();
    number2shmemname(namebuf,memid);

/*
    The 'shm_open' creates and opens a new, or opens an existing, POSIX
    shared memory object. A POSIX shared memory object is in effect a
    handle which can be used by unrelated processes to mmap(2) the same
    region of shared memory.
*/
    fd=shm_open(namebuf,O_CREAT|O_EXCL|O_RDWR,0666);
    if (fd>=0) break; /* successfully created a named shared memory segment */
    if (errno==EEXIST) continue; /* try again with a different name */
    perror("shm_open() error");
    printf("Cannot create a POSIX shared memory segment '%s', errno=%d\n",namebuf,errno);
    Fatal_error(911);
    return -1; /* never returns */
  }
  if (i==0) {
    printf("Cannot create a unique POSIX shared memory segment name even after %d attempts\n",ntries);
    Fatal_error(911);
    return -1;  /* never returns */
  }

  /* The number and the name are chosen, now set the size */
  if (ftruncate(fd,sz)==-1) {
    perror("ftruncate() error");
    printf("Cannot set size (%ld) for a POSIX shared memory segment '%s', errno=%d\n",(long int)sz,namebuf,errno);
    Fatal_error(911);
    return -1; /* never returns */
  }

  /* Close the descriptor */
  close(fd);

  return memid;
}

/* -------------------------------------------------------- *\
   Wrapper function to attach a POSIX shared memory segment
\* -------------------------------------------------------- */
void *Shmat(int memid, void *dummy1, int dummy2)

/*  note that only the 1st argument is used below  */

{
  char namebuf[NAME_MAX];
  struct stat statbuf;
  int fd;
  void* addr;

  number2shmemname(namebuf,memid);

  fd=shm_open(namebuf,O_RDWR,0);
  if (fd==-1) {
    perror("shm_open() error");
    printf("Cannot open a POSIX shared memory segment '%s', errno=%d\n",
             namebuf,errno);
    Fatal_error(911);
    return 0; /* never returns */
  }
  if (fstat(fd,&statbuf)!=0) {
    perror("fstat() error");
    printf("Cannot fstat() a POSIX shared memory segment '%s', errno=%d\n",namebuf,errno);
    Fatal_error(911);
    return 0; /* never returns */
  }
  if (statbuf.st_size==0) {
    printf("fstat() on a POSIX shared memory segment '%s' returned 0 (should not happen).\n",namebuf);
    Fatal_error(911);
    return 0; /* never returns */
  }

  addr=mmap(NULL,statbuf.st_size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);

  if (addr==MAP_FAILED) {
    perror("mmap() error");
    printf("Cannot mmap() a POSIX shared memory segment '%s' of size %lu, errno=%d\n",
           namebuf,(unsigned long)statbuf.st_size,errno);
    Fatal_error(911);
    return 0; /* never returns */
  }

  close(fd);

  return addr;
}

/* ------------------------------------------------------------ *\
   Wrapper function to pre-delete a POSIX shared memory segment
\* ------------------------------------------------------------ */
int Shmctl(int memid, int cmd, struct shmid_ds *buff)

/*  note that only the 1st argument is used below  */

{
  int iret;
  char namebuf[NAME_MAX];

  number2shmemname(namebuf,memid);

/*
     The 'shm_unlink' removes a shared memory object name, and, once all
     processes have unmapped the object, de-allocates and destroys the
     contents of the associated memory region.
     So, 'shm_unlink' is not a deletion, until the DDI job is finished!
     However, the "file" disappears from the memory subsystem, so it
     is hard for users to see that POSIX memory is in use.
*/

  iret = shm_unlink(namebuf);

  if (iret != 0) {
    fprintf(stdout,"Cannot unlink POSIX shared memory segment '%s'\n",namebuf);
    fprintf(stdout,"error value=%i\n",iret);
    fflush(stdout);
    Fatal_error(911);
    return; /* never returns */
  }
  return iret;
}

#endif
     /*   the line above closes the POSIX shared memory implementation   */

/*
             SEMAPHORES.
     There is only a SYSV semaphore implementation here,
     so this must be compiled for POSIX, too
*/

# if defined USE_SHMEM_SYSV || defined USE_SHMEM_POSIX

/* --------------------------- *\
   Wrapper function for semget
\* --------------------------- */
   int Semget(key_t key,int nsems,int semflg) {
      int ret;

      if((ret = semget(key,nsems,semflg)) == -1) {
         fprintf(stdout,"%s: semget return an error.\n",DDI_Id());
         switch(errno) {
           case EACCES: fprintf(stdout," semget errno=EACCES.\n"); break;
           case EINVAL: fprintf(stdout," semget errno=EINVAL.\n"); break;
           case ENOENT: fprintf(stdout," semget errno=ENOENT.\n"); break;
           case ENOSPC: fprintf(stdout," semget errno=ENOSPC -- check system limit for sysv semaphores.\n"); break;
           case ENOMEM: fprintf(stdout," semget errno=ENOMEM.\n"); break;
           case EEXIST: fprintf(stdout," semget errno=EEXIST.\n"); break;
           default:
             fprintf(stdout," unusual semget errno=%d\n",errno); break;
         }
         Fatal_error(911);
      }

      return ret;
   }


/* -------------------------- *\
   Wrapper function for semop
\* -------------------------- */
   int Semop(int semid,struct sembuf *opers,size_t nops) {
      int ret;

/*
    Note that it is apparently possible for a run to receive an
    interupt signal in the middle of checking a semaphore.  This
    will interupt prior to learning the state of the semaphore.
    The "fix" is to just retry (recursive call) to learn what the
    state of the semaphore really is.
*/
      if((ret = semop(semid,opers,nops)) == -1) {
         fprintf(stdout,"%s: semop return an error performing %i operation(s) on semid %i.\n",DDI_Id(),(int) nops,semid);
         switch(errno) {
           case EFBIG:  fprintf(stdout," semop errno=EFBIG.\n"); break;
           case E2BIG:  fprintf(stdout," semop errno=E2BIG.\n"); break;
           case EINTR:  fprintf(stdout," semop errno=EINTR.\n");
                        return Semop(semid,opers,nops);
           case EINVAL: fprintf(stdout," semop errno=EINVAL.\n"); break;
           case EACCES: fprintf(stdout," semop errno=EACCES.\n"); break;
           case EAGAIN: fprintf(stdout," semop errno=EAGAIN.\n"); break;
           case ENOSPC: fprintf(stdout," semop errno=ENOSPC.\n"); break;
           case ERANGE: fprintf(stdout," semop errno=ERANGE.\n"); break;
           case EFAULT: fprintf(stdout," semop errno=EFAULT.\n"); break;
           default:
              fprintf(stdout," unusual semop errno=%d\n",errno); break;
         }
      }

      return ret;
   }


/* ------------ *\
   DDI_Sem_oper
\* ------------ */
   int DDI_Sem_oper(int semid,int semnum,int semop) {
      struct sembuf op;

      op.sem_op  = semop;
      op.sem_num = semnum;
      op.sem_flg = 0;

      return Semop(semid,&op,1);
   }


/* ------------------------------------------------- *\
   DDI_Sem_acquire
   Acquire a user defined access level to a resource
\* ------------------------------------------------- */
   void DDI_Sem_acquire(int semid,int semnum,int access) {
      struct sembuf op;

      op.sem_op = -access;
      op.sem_num = semnum;
      op.sem_flg = 0;

      Semop(semid,&op,1);
   }


/* --------------------------------------------- *\
   DDI_Sem_release
   Release an acquire access to a given resource
\* --------------------------------------------- */
   void DDI_Sem_release(int semid,int semnum,int access) {
      struct sembuf op;

      op.sem_op  = access;
      op.sem_num = semnum;
      op.sem_flg = 0;

      Semop(semid,&op,1);
   }


/* ------------------------------------------- *\
   DDI_Sem_remove
   Remove a System V semaphore from the system
\* ------------------------------------------- */
   void DDI_Sem_remove(int semid) {
    # if defined MACOSX
      union semun arg;
      semctl(semid,0,IPC_RMID,arg);
    # else
      semctl(semid,0,IPC_RMID);
    # endif
   }

# endif
     /*   the line above closes the SYSV semaphore implementation   */

     /*   below takes care of not using either kind of shared memory  */
# else
   void DDI_Sysv_dummy(void) { return;  }
# endif
