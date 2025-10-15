#define CHECK_CUDECOMP_EXIT(f) if (f /= CUDECOMP_RESULT_SUCCESS) call exit(1)

program main
use cudafor
use cudecomp
use cufft
use mpi
use velocity 
use phase
use particles
use param
use mpivar
use cudecompvar

implicit none
! timer for scaling test
real :: t_start, t_end, elapsed
! grid dimensions
integer :: comm_backend
integer :: pr, pc
! cudecomp
! cuFFT
integer :: planXf, planXb
integer :: planY, planZ
integer :: batchsize
integer :: status
integer :: i,j,k,il,jl,kl,ig,jg,kg,t
integer :: im,ip,jm,jp,km,kp,last
integer :: inY,enY,inZ,enZ
!beginDEBUG
integer, parameter :: Mx = 1, My = 2, Mz = 1
!endDEBUG
real(8), device, allocatable :: kx_d(:)
! working arrays
complex(4), allocatable :: psi(:)
real(4), allocatable :: ua(:,:,:)
real(4), allocatable :: uaa(:,:,:)
real(4), allocatable :: psi_real(:)
complex(4), device, allocatable :: psi_d(:)
complex(8), pointer, device, contiguous :: work_d(:), work_halo_d(:)
complex(4), pointer, device, contiguous :: work_d_d2z(:)
character(len=40) :: namefile
character(len=4) :: itcount
! Code variables

real(8)::err,maxErr

complex(4), device, pointer :: phi3d(:,:,:)
real(8) :: k2
!integer :: il, jl, ig, jg
integer :: offsets(3), xoff, yoff
integer :: np(3)

! Enable or disable phase field (acceleration eneabled by default)
#define phiflag 0
! Enable or disable particle Lagrangian tracking (tracers)
#define partflag 0

!########################################################################################################################################
! 1. INITIALIZATION OF MPI AND cuDECOMP AUTOTUNING : START
!########################################################################################################################################
! MPI initialization, put in rank the local MPI rank number and ranks total number
! Same procedura defined in the cuDecomp documentation
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
call mpi_comm_size(MPI_COMM_WORLD, ranks, ierr)

call mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, localComm, ierr)
call mpi_comm_rank(localComm, localRank, ierr)
ierr = cudaSetDevice(localRank) !assign GPU to MPI rank

! Define grid and decomposition
call readinput

! number of process along each direction (r means y and c means z)
pr = 0
pc = 0
halo_ext=1
! comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
comm_backend = 0 ! Enable full autotuning

CHECK_CUDECOMP_EXIT(cudecompInit(handle, MPI_COMM_WORLD))

! config is a struct and pr and pc are the number of pencils along the two directions
! gdims is the global grid
! create an uninitialized configuration struct and initialize it to defaults using cudecompGridDescConfigSetDefaults. 
! Initializing to default values is required to ensure no entries are left uninitialized.
CHECK_CUDECOMP_EXIT(cudecompGridDescConfigSetDefaults(config))
pdims = [pr, pc] !pr and pc are the number of pencil along the different directions
config%pdims = pdims
! gdims = [nx, ny, nz]
! config%gdims = gdims
halo = [0, halo_ext, halo_ext] ! no halo along x neeed because is periodic and in physical space i have x-pencil
! for transpositions
config%transpose_comm_backend = comm_backend
config%transpose_axis_contiguous = .true.
! for halo exchanges
config%halo_comm_backend = CUDECOMP_HALO_COMM_MPI
! Setting for periodic halos in all directions (non required to be in config)
halo_periods = [.true., .true., .true.]

! create spectral grid descriptor first to select pdims for optimal transposes
gdims = [nx/2+1, ny, nz]
config%gdims = gdims

! Set up autotuning options for spectral grid (transpose related settings)
CHECK_CUDECOMP_EXIT(cudecompGridDescAutotuneOptionsSetDefaults(options))
options%dtype = CUDECOMP_FLOAT_COMPLEX
if (comm_backend == 0) then
   options%autotune_transpose_backend = .true.
   options%autotune_halo_backend = .false.
endif
options%transpose_use_inplace_buffers = .true.
options%transpose_input_halo_extents(:, 1) = halo
options%transpose_output_halo_extents(:, 4) = halo

CHECK_CUDECOMP_EXIT(cudecompGridDescCreate(handle, grid_descD2Z, config, options))

! create physical grid descriptor
! take previous config and modify the global grid (nx instead of nx/2+1)
! reset transpose_comm_backend to default value to avoid picking up possible nvshmem
! transpose backend selection (this impacts how workspaces are allocated)
gdims = [nx, ny, nz]
config%gdims = gdims
config%transpose_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P

! Set up autotuning options for physical grid (halo related settings)
CHECK_CUDECOMP_EXIT(cudecompGridDescAutotuneOptionsSetDefaults(options))
options%dtype = CUDECOMP_DOUBLE_COMPLEX
if (comm_backend == 0) then
   options%autotune_halo_backend = .true.
endif
options%halo_extents(:) = halo
options%halo_periods(:) = halo_periods
options%halo_axis = 1
CHECK_CUDECOMP_EXIT(cudecompGridDescCreate(handle, grid_desc, config, options))

!Print information on configuration
! issue with NVHPC > 25.X, simply avoid this call (not critical)
!if (rank == 0) then
!   write(*,"(' Running on ', i0, ' x ', i0, ' process grid ...')") config%pdims(1), config%pdims(2)
!   write(*,"(' Using ', a, ' transpose backend ...')") &
!         cudecompTransposeCommBackendToString(config%transpose_comm_backend)
!   write(*,"(' Using ', a, ' halo backend ...')") &
!         cudecompHaloCommBackendToString(config%halo_comm_backend)
!endif

! Get pencil info for the grid descriptor in the physical space
! This function returns a pencil struct (piX, piY or piZ) that contains the shape, global lower and upper index bounds (lo and hi), 
! size of the pencil, and an order array to indicate the memory layout that will be used (to handle permuted, axis-contiguous layouts).
! Additionally, there is a halo_extents data member that indicates the depth of halos for the pencil, by axis.
! If no halo regions are necessary, a NULL pointer can be provided in place of this array (or omitted)
! Pencil info in x-configuration present in PiX (shape,lo,hi,halo_extents,size)
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piX, 1, halo))
nElemX = piX%size !<- number of total elments in x-configuratiion (including halo)
! Pencil info in Y-configuration present in PiY
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piY, 2))
nElemY = piY%size
! Pencil info in Z-configuration present in PiZ
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piZ, 3))
nElemZ = piZ%size

! Get workspace sizes for transpose (1st row, not used) and halo (2nd row, used)
CHECK_CUDECOMP_EXIT(cudecompGetTransposeWorkspaceSize(handle, grid_desc, nElemWork))
CHECK_CUDECOMP_EXIT(cudecompGetHaloWorkspaceSize(handle, grid_desc, 1, halo, nElemWork_halo))

CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_descD2Z, piX_d2z, 1,halo))
nElemX_d2z = piX_d2z%size !<- number of total elments in x-configuratiion (include halo)
! Pencil info in Y-configuration present in PiY
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_descD2Z, piY_d2z, 2))
nElemY_d2z = piY_d2z%size
! Pencil info in Z-configuration present in PiZ
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_descD2Z, piZ_d2z, 3))
nElemZ_d2z = piZ_d2z%size
! Get workspace sizes for transpose (1st row,used) and halo (2nd row, not used)
CHECK_CUDECOMP_EXIT(cudecompGetTransposeWorkspaceSize(handle, grid_descD2Z, nElemWork_d2z))
CHECK_CUDECOMP_EXIT(cudecompGetHaloWorkspaceSize(handle, grid_descD2Z, 1, halo, nElemWork_halo_d2z))

! Get the global rank of neighboring processes in PiX config (required only for particles)
CHECK_CUDECOMP_EXIT(cudecompGetShiftedRank(handle, grid_desc, 1, 2, 1  , .true. , nidp1y))
CHECK_CUDECOMP_EXIT(cudecompGetShiftedRank(handle, grid_desc, 1, 2, -1 , .true. , nidm1y))
CHECK_CUDECOMP_EXIT(cudecompGetShiftedRank(handle, grid_desc, 1, 3, 1  , .true. , nidp1z))
CHECK_CUDECOMP_EXIT(cudecompGetShiftedRank(handle, grid_desc, 1, 3, -1 , .true. , nidm1z))


! CUFFT initialization -- Create Plans
! Forward 1D FFT in X: D2Z
batchSize = piX_d2z%shape(2)*piX_d2z%shape(3) !<- number of FFT (from x-pencil dimension)
status = cufftPlan1D(planXf, nx, CUFFT_R2C, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating X plan Forward'

! Backward 1D FFT in X: Z2D
batchSize = piX_d2z%shape(2)*piX_d2z%shape(3) !<- number of FFT (from x-pencil dimension)
status = cufftPlan1D(planXb, nx, CUFFT_C2R, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating X plan Backward'

! it's always 2 and 3 because y-pencil have coordinates y,z,x
batchSize = piY_d2z%shape(2)*piY_d2z%shape(3)
status = cufftPlan1D(planY, ny, CUFFT_C2C, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating Y plan Forward & Backward'

! it's always 2 and 3 because y-pencil have coordinates z,y,x
batchSize = piZ_d2z%shape(2)*piZ_d2z%shape(3)
status = cufftPlan1D(planZ, nz, CUFFT_C2C, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating Z plan Forward & Backward'



! define grid
allocate(x(nx),x_ext(nx+1))
x(1)= 0
do i = 2, nx
   x(i) = x(i-1) + dx
enddo

x_ext(1:nx) = x(1:nx)
x_ext(nx+1) = lx

! Offsets in the X-pencil decomposition
pix_yoff = pix%lo(2)-1
pix_zoff = pix%lo(3)-1

! Boundaries of each PiX pencil
yinf = x_ext( pix%lo(2) )
ysup = x_ext( pix%hi(2) + 1 )
zinf = x_ext( pix%lo(3) )
zsup = x_ext( pix%hi(3) + 1 )   

! Physical size of each PiX pencil
lyloc = ysup - yinf
lzloc = zsup - zinf

! define wavenumbers
allocate(kx(nx))
do i = 1, nx/2
   kx(i) = (i-1)*(twoPi/lx)
enddo
do i = nx/2+1, nx
   kx(i) = (i-1-nx)*(twoPi/LX)
enddo
! allocate k_d on the device (later on remove and use OpenACC + managed memory?)
allocate(kx_d, source=kx)

allocate(mysin(nx), mycos(nx))
do i=1,nx
   ! compute here the sin to avoid multiple computations of sin
   mysin(i)=sin(k0*x(i)+dx/2)
   ! compute here the cos to avoid multiple computations of cos
   mycos(i)=cos(k0*x(i)+dx/2)
enddo

! Initial distribution of particles among the processes
nploc = npart/ranks
nplocmax = nploc*2

!########################################################################################################################################
! 1. INITIALIZATION AND cuDECOMP AUTOTUNING : END
!########################################################################################################################################



!########################################################################################################################################
! START STEP 2: ALLOCATE ARRAYS
!########################################################################################################################################

! allocate arrays
allocate(psi(max(nElemX, nElemY, nElemZ))) !largest among the pencil
allocate(psi_real(max(nElemX, nElemY, nElemZ))) !largest among the pencil
allocate(psi_d(max(nElemX_d2z, nElemY_d2z, nElemZ_d2z))) ! phi on device
allocate(ua(nx, piX%shape(2), piX%shape(3)))

! Pressure variable
allocate(rhsp(piX%shape(1), piX%shape(2), piX%shape(3))) 
allocate(p(piX%shape(1), piX%shape(2), piX%shape(3))) 

!allocate variables
!NS variables
allocate(u(piX%shape(1),piX%shape(2),piX%shape(3)),v(piX%shape(1),piX%shape(2),piX%shape(3)),w(piX%shape(1),piX%shape(2),piX%shape(3))) !velocity vector
! allocate(ustar(piX%shape(1),piX%shape(2),piX%shape(3)),vstar(piX%shape(1),piX%shape(2),piX%shape(3)),wstar(piX%shape(1),piX%shape(2),piX%shape(3))) ! provisional velocity field
allocate(rhsu(piX%shape(1),piX%shape(2),piX%shape(3)),rhsv(piX%shape(1),piX%shape(2),piX%shape(3)),rhsw(piX%shape(1),piX%shape(2),piX%shape(3))) ! right hand side u,v,w
allocate(rhsu_o(piX%shape(1),piX%shape(2),piX%shape(3)),rhsv_o(piX%shape(1),piX%shape(2),piX%shape(3)),rhsw_o(piX%shape(1),piX%shape(2),piX%shape(3))) ! right hand side u,v,w
allocate(div(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(rhspp(piX%shape(1),piX%shape(2),piX%shape(3))) ! right hand side pressure in FP32
allocate(pp(piX%shape(1),piX%shape(2),piX%shape(3)))     ! pressure in FP32

!PFM variables
#if phiflag == 1
   allocate(phi(piX%shape(1),piX%shape(2),piX%shape(3)),rhsphi(piX%shape(1),piX%shape(2),piX%shape(3)),rhsphi_o(piX%shape(1),piX%shape(2),piX%shape(3)))
   allocate(psidi(piX%shape(1),piX%shape(2),piX%shape(3)))
   allocate(tanh_psi(piX%shape(1),piX%shape(2),piX%shape(3)))
   allocate(normx(piX%shape(1),piX%shape(2),piX%shape(3)),normy(piX%shape(1),piX%shape(2),piX%shape(3)),normz(piX%shape(1),piX%shape(2),piX%shape(3)))
   allocate(fxst(piX%shape(1),piX%shape(2),piX%shape(3)),fyst(piX%shape(1),piX%shape(2),piX%shape(3)),fzst(piX%shape(1),piX%shape(2),piX%shape(3))) ! surface tension forces
#endif

! allocate arrays for transpositions and halo exchanges 
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_desc, work_d, nElemWork))
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_desc, work_halo_d, nElemWork_halo))
! allocate arrays for transpositions
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_descD2Z, work_d_d2z, nElemWork_d2z))
!CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_descD2Z, work_halo_d_d2z, nElemWork_halo_d2z))

#if partflag == 1
! Particle variables
allocate(part(1:nplocmax,1:ninfop))
allocate(partbuff(1:nplocmax,1:ninfop))
allocate(vec_p(1:nplocmax))
allocate(order_p(1:nplocmax))
allocate(buffvar1(1:ninfop,1:nploc))
allocate(buffvar2(1:ninfop,1:nploc))
#endif
!########################################################################################################################################
! END STEP2: ALLOCATE ARRAYS
!########################################################################################################################################



!########################################################################################################################################
! START STEP 3: FLOW, PHASE FIELD AND PARTICLES INIT
!########################################################################################################################################
! 3.1 Read/initialize from data without halo grid points (avoid out-of-bound if reading usin MPI I/O)
! 3.2 Call halo exchnages along Y and Z for u,v,w and phi
if (restart .eq. 0) then !fresh start Taylor Green or read from file in init folder
   if (rank.eq.0) write(*,*) "Initialize velocity field (fresh start)"
   if (inflow .eq. 0) then
      if (rank.eq.0) write(*,*) "Initialize Taylor-green"
      do k = 1+halo_ext, piX%shape(3)-halo_ext
         kg = piX%lo(3) + k - 1 
         do j = 1+halo_ext, piX%shape(2)-halo_ext
            jg = piX%lo(2) + j - 1 
            do i = 1, piX%shape(1)
               u(i,j,k) =   sin(x(i)-dx/2)*cos(x(jg))*cos(x(kg))
               v(i,j,k) =  -cos(x(i))*sin(x(jg)-dx/2)*cos(x(kg))
               w(i,j,k) =  0.d0
            enddo
         enddo
      enddo
   endif
   if (inflow .eq. 1) then
      if (rank.eq.0)  write(*,*) "Initialize from data"
      call readfield(1)
      call readfield(2)
      call readfield(3)
   endif
endif
if (restart .eq. 1) then !restart, ignore inflow and read the tstart field 
   if (rank.eq.0)  write(*,*) "Initialize velocity field (from output folder), iteration:", tstart
   call readfield_restart(tstart,1)
   call readfield_restart(tstart,2)
   call readfield_restart(tstart,3)
endif

! update halo cells along y and z directions (enough only if pr and pc are non-unitary)
!$acc host_data use_device(u)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 
!$acc host_data use_device(v)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 
!$acc host_data use_device(w)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 

! initialize phase-field
#if phiflag == 1
   if (restart .eq. 0) then
      if (rank.eq.0) write(*,*) 'Initialize phase field (fresh start)'
      if (inphi .eq. 0) then
         if (rank.eq.0) write(*,*) 'Spherical drop'
         do k = 1+halo_ext, piX%shape(3)-halo_ext
            kg = piX%lo(3) + k - 1 
            do j = 1+halo_ext, piX%shape(2)-halo_ext
               jg = piX%lo(2) + j - 1 
               do i = 1, piX%shape(1)
                  pos=(x(i)-lx/2)**2d0 +  (x(jg)-lx/2)**2d0 + (x(kg)-lx/2)**2d0
                  phi(i,j,k) = 0.5d0*(1.d0-tanh((sqrt(pos)-radius)/2/eps))
               enddo
            enddo
         enddo
      endif
      if (inphi .eq. 1) then
         if (rank.eq.0)  write(*,*) "Initialize phase-field from data"
         call readfield(5)
      endif
   endif
   if (restart .eq. 1) then
      write(*,*) "Initialize phase-field (restart, from output folder), iteration:", tstart
      call readfield_restart(tstart,5)
   endif
   ! update halo
   !$acc host_data use_device(phi)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
#endif

#if partflag == 1
  if (restart .eq. 0) then
     if (rank.eq.0) write(*,*) 'Initialize particles (fresh start)'
     if (inpart .eq. 1) then
        if (rank.eq.0) write(*,*) 'Random Position whole Domain'
        call particlegenerator(inpart)
     endif
      if (inpart .eq. 2) then
         if (rank.eq.0)  write(*,*) 'Random Position in Drops'
         call particlegenerator(inpart)
      endif
  endif
 !   if (restart .eq. 1) then
 !      write(*,*) "Initialize phase-field (restart, from output folder), iteration:", tstart
 !      call readfield_restart(tstart,5)
 !   endif
#endif

!Save initial fields (only if a fresh start)
if (restart .eq. 0) then
   if (rank.eq.0) write(*,*) "Save initial fields"
   call writefield(tstart,1)
   call writefield(tstart,2)
   call writefield(tstart,3)
   call writefield(tstart,4)
   #if phiflag == 1
      call writefield(tstart,5)
   #endif
endif

!########################################################################################################################################
! END STEP 3: FLOW, PHASE FIELD AND PARTICLES INIT
!########################################################################################################################################



!########################################################################################################################################
! START TEMPORAL LOOP: STEP 4 to 9 REPEATED AT EVERY TIME STEP
!########################################################################################################################################

! First step use Euler
alpha=1.0d0
beta=0.0d0
gumax=1.d0
tstart=tstart+1
gamma=1.d0*gumax
!$acc data copyin(piX)
!$acc data create(rhsu_o, rhsv_o, rhsw_o)
!$acc data copyin(mysin, mycos)
#if partflag == 1
!$acc data copy(part,partbuff)
!$acc data create(vec_p, order_p)
!$acc data create(buffvar1,buffvar2)
#endif
call cpu_time(t_start)

! Start temporal loop
do t=tstart,tfin
   ! Create custom label for each marker (profiling)
    write(itcount,'(i4)') t
   ! Range with custom  color (uncomment for profiling)
   ! call nvtxStartRange("Iteration "//itcount,t)

    if (rank.eq.0) write(*,*) "Time step",t,"of",tfin
    call cpu_time(times)


   !########################################################################################################################################
   ! START STEP 4: PARTICLES (TRACERS)
   !########################################################################################################################################
   #if partflag == 1
      ! Operations:
      ! 4.1 Perform Interpolation (consider passing to trilinear: more accurate, but way more expensive)
      ! 4.2 Integrate with Adams-Bashforth
      ! 4.3 Order and transfer in y
      ! 4.4 Order and transfer in z 
      ! 4.5 Check Leakage of Particles

      call LinearInterpolation()
      ! Particle Tracker Integration
      ! Two-Step Adams-Bashfort (Euler for first step)
      !$acc parallel loop collapse(2) default(present)
      do j = 0, 2
         do i = 1, nploc
            part(i,Ixp+j)=part(i,Ixp+j)+&
                     dt*(alpha*part(i,Iup+j)-beta*part(i,Iup1+j))
         enddo
      enddo

      ! Transfer in y
      call SortPartY()
      call CountPartTransfY()
      call SendPartUP(2)
      call SendPartDOWN(2)
 
      ! Transfer in z
      call SortPartZ()
      call CountPartTransfZ()
      call SendPartUP(3)
      call SendPartDOWN(3)

      ! Check Particles esacping the domain (leakage)
      call ParticlesLeakage()
      ! Shift data for next step
      !$acc parallel loop collapse(2) default(present)
      do j = 0, 2
         do i = 1, nploc
            part(i,Iup1+j)=part(i,Iup+j)
         enddo
      enddo
      write(*,*) 'rank',rank, 'nploc',nploc
   #endif
   !########################################################################################################################################
   ! END STEP 4: PARTICLES
   !########################################################################################################################################

   
   ! (uncomment for profiling)
   ! call nvtxStartRange("Phase-field")
   !########################################################################################################################################
   ! START STEP 5: PHASE-FIELD SOLVER (EXPLICIT)
   !########################################################################################################################################
   #if phiflag == 1
      !$acc kernels
      do k=1, piX%shape(3)
         do j=1, piX%shape(2)
            do i=1,nx
               ! compute distance function psi (used to compute normals)
               val = min(phi(i,j,k),1.0d0) ! avoid machine precision overshoots in phi that leads to problem with log
               psidi(i,j,k) = eps*log((val+enum)/(1.d0-val+enum))
               ! compute here the tanh of distance function psi (used in the sharpening term) to avoid multiple computations of tanh
               tanh_psi(i,j,k) = tanh(0.5d0*psidi(i,j,k)*epsi)
            enddo
         enddo
      enddo
      !$acc end kernels

      gamma=1.d0*gumax
      !$acc parallel loop tile(16,4,2)
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               ! 4.1 RHS computation
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               if (ip .gt. nx) ip=1
               if (im .lt. 1) im=nx
               ! convective (first three lines) and diffusive (last three lines)
               rhsphi(i,j,k) =   &
                     - (u(ip,j,k)*0.5d0*(phi(ip,j,k)+phi(i,j,k)) - u(i,j,k)*0.5d0*(phi(i,j,k)+phi(im,j,k)))*dxi  &  
                     - (v(i,jp,k)*0.5d0*(phi(i,jp,k)+phi(i,j,k)) - v(i,j,k)*0.5d0*(phi(i,j,k)+phi(i,jm,k)))*dxi  &  
                     - (w(i,j,kp)*0.5d0*(phi(i,j,kp)+phi(i,j,k)) - w(i,j,k)*0.5d0*(phi(i,j,k)+phi(i,j,km)))*dxi  &  
                           + gamma*(eps*(phi(ip,j,k)-2.d0*phi(i,j,k)+phi(im,j,k))*ddxi + &                   
                                    eps*(phi(i,jp,k)-2.d0*phi(i,j,k)+phi(i,jm,k))*ddxi + &                   
                                    eps*(phi(i,j,kp)-2.d0*phi(i,j,k)+phi(i,j,km))*ddxi)                      
               ! 4.1.3. Compute normals for sharpening term (gradient)
               normx(i,j,k) = (psidi(ip,j,k) - psidi(im,j,k))
               normy(i,j,k) = (psidi(i,jp,k) - psidi(i,jm,k))
               normz(i,j,k) = (psidi(i,j,kp) - psidi(i,j,km))
            enddo
         enddo
      enddo

      ! Update normx,normy and normz halos, required to then compute normal derivative
      !$acc host_data use_device(normx)
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normx, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normx, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      !$acc end host_data 
      !$acc host_data use_device(normy)
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normy, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normy, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      !$acc end host_data 
      !$acc host_data use_device(normz)
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normz, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normz, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      !$acc end host_data 

      ! 4.1.3. Compute Sharpening term (gradient)
      ! Substep 2: Compute normals (1.e-16 is a numerical tollerance to avoid 0/0)
      !$acc kernels
      do k=1, piX%shape(3)
         do j=1, piX%shape(2)
            do i=1,nx
               normag = 1.d0/(sqrt(normx(i,j,k)*normx(i,j,k) + normy(i,j,k)*normy(i,j,k) + normz(i,j,k)*normz(i,j,k)) + enum)
               normx(i,j,k) = normx(i,j,k)*normag
               normy(i,j,k) = normy(i,j,k)*normag
               normz(i,j,k) = normz(i,j,k)*normag
            enddo
         enddo
      enddo
      !$acc end kernels

      ! Compute sharpening term
      !$acc kernels
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               if (ip .gt. nx) ip=1
               if (im .lt. 1) im=nx
               ! ACDI with pre-computed tanh                                
               rhsphi(i,j,k)=rhsphi(i,j,k)-gamma*((0.25d0*(1.d0-tanh_psi(ip,j,k)*tanh_psi(ip,j,k))*normx(ip,j,k) - &
                                                   0.25d0*(1.d0-tanh_psi(im,j,k)*tanh_psi(im,j,k))*normx(im,j,k))*0.5d0*dxi + &
                                                  (0.25d0*(1.d0-tanh_psi(i,jp,k)*tanh_psi(i,jp,k))*normy(i,jp,k) - &
                                                   0.25d0*(1.d0-tanh_psi(i,jm,k)*tanh_psi(i,jm,k))*normy(i,jm,k))*0.5d0*dxi + &
                                                  (0.25d0*(1.d0-tanh_psi(i,j,kp)*tanh_psi(i,j,kp))*normz(i,j,kp) - &
                                                   0.25d0*(1.d0-tanh_psi(i,j,km)*tanh_psi(i,j,km))*normz(i,j,km))*0.5d0*dxi)
            enddo
         enddo
      enddo
      !$acc end kernels

      ! 4.2 Get phi at n+1 using AB2
      !$acc kernels
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               phi(i,j,k) = phi(i,j,k) + dt*(alpha*rhsphi(i,j,k)-beta*rhsphi_o(i,j,k))
               rhsphi_o(i,j,k)=rhsphi(i,j,k)
            enddo
         enddo
      enddo
      !$acc end kernels

      ! 4.3 Call halo exchnages along Y and Z for phi 
      !$acc host_data use_device(phi)
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      !$acc end host_data 
   #endif

   ! (uncomment for profiling)
   ! call nvtxEndRange

   !########################################################################################################################################
   ! END STEP 5: PHASE-FIELD SOLVER 
   !########################################################################################################################################


   !########################################################################################################################################
   ! START STEP 6: USTAR COMPUTATION (PROJECTION STEP)
   !########################################################################################################################################
   ! 6.1 compute rhs 
   ! 6.2 obtain ustar and store old rhs in rhs_o
   ! 6.3 Call halo exchnages along Y and Z for u,v,w

   ! (uncomment for profiling)
   ! call nvtxStartRange("Projection")
   ! 6.1a Convective and diffusiver terms NS
   ! Loop on inner nodes
   !$acc parallel loop tile(16,4,2) 
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            ip=i+1
            jp=j+1
            kp=k+1
            im=i-1
            jm=j-1
            km=k-1
            ! Manual periodicity ony along x (x-pencil), along y and z directions use halos
            if (ip .gt. nx) ip=1  
            if (im .lt. 1) im=nx
            ! compute the products (conservative form)
            h11 = (u(ip,j,k)+u(i,j,k))*(u(ip,j,k)+u(i,j,k))     - (u(i,j,k)+u(im,j,k))*(u(i,j,k)+u(im,j,k))
            h12 = (u(i,jp,k)+u(i,j,k))*(v(i,jp,k)+v(im,jp,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k))
            h13 = (u(i,j,kp)+u(i,j,k))*(w(i,j,kp)+w(im,j,kp))   - (u(i,j,k)+u(i,j,km))*(w(i,j,k)+w(im,j,k))
            h21 = (u(ip,j,k)+u(ip,jm,k))*(v(ip,j,k)+v(i,j,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k))
            h22 = (v(i,jp,k)+v(i,j,k))*(v(i,jp,k)+v(i,j,k))     - (v(i,j,k)+v(i,jm,k))*(v(i,j,k)+v(i,jm,k))
            h23 = (w(i,j,kp)+w(i,jm,kp))*(v(i,j,kp)+v(i,j,k))   - (w(i,j,k)+w(i,jm,k))*(v(i,j,k)+v(i,j,km))
            h31 = (w(ip,j,k)+w(i,j,k))*(u(ip,j,k)+u(ip,j,km))   - (w(i,j,k)+w(im,j,k))*(u(i,j,k)+u(i,j,km))
            h32 = (v(i,jp,k)+v(i,jp,km))*(w(i,jp,k)+w(i,j,k))   - (v(i,j,k)+v(i,j,km))*(w(i,j,k)+w(i,jm,k))
            h33 = (w(i,j,kp)+w(i,j,k))*(w(i,j,kp)+w(i,j,k))     - (w(i,j,k)+w(i,j,km))*(w(i,j,k)+w(i,j,km))
            ! compute the derivative
            h11=h11*0.25d0*dxi
            h12=h12*0.25d0*dxi
            h13=h13*0.25d0*dxi
            h21=h21*0.25d0*dxi
            h22=h22*0.25d0*dxi
            h23=h23*0.25d0*dxi
            h31=h31*0.25d0*dxi
            h32=h32*0.25d0*dxi
            h33=h33*0.25d0*dxi
            ! add to the rhs
            rhsu(i,j,k)=-(h11+h12+h13)
            rhsv(i,j,k)=-(h21+h22+h23)
            rhsw(i,j,k)=-(h31+h32+h33)
            ! viscos term
            h11 = mu*(u(ip,j,k)-2.d0*u(i,j,k)+u(im,j,k))*ddxi
            h12 = mu*(u(i,jp,k)-2.d0*u(i,j,k)+u(i,jm,k))*ddxi
            h13 = mu*(u(i,j,kp)-2.d0*u(i,j,k)+u(i,j,km))*ddxi
            h21 = mu*(v(ip,j,k)-2.d0*v(i,j,k)+v(im,j,k))*ddxi
            h22 = mu*(v(i,jp,k)-2.d0*v(i,j,k)+v(i,jm,k))*ddxi
            h23 = mu*(v(i,j,kp)-2.d0*v(i,j,k)+v(i,j,km))*ddxi
            h31 = mu*(w(ip,j,k)-2.d0*w(i,j,k)+w(im,j,k))*ddxi
            h32 = mu*(w(i,jp,k)-2.d0*w(i,j,k)+w(i,jm,k))*ddxi
            h33 = mu*(w(i,j,kp)-2.d0*w(i,j,k)+w(i,j,km))*ddxi
            rhsu(i,j,k)=rhsu(i,j,k)+(h11+h12+h13)*rhoi
            rhsv(i,j,k)=rhsv(i,j,k)+(h21+h22+h23)*rhoi
            rhsw(i,j,k)=rhsw(i,j,k)+(h31+h32+h33)*rhoi
            ! NS forcing
            kg = piX%lo(3) + k - 1 
            jg = piX%lo(2) + j - 1
            ! ABC forcing
            rhsu(i,j,k)= rhsu(i,j,k) + f3*mysin(kg)+f2*mycos(jg)
            rhsv(i,j,k)= rhsv(i,j,k) + f1*mysin(i)+f3*mycos(kg)
            rhsw(i,j,k)= rhsw(i,j,k) + f2*mysin(jg)+f1*mycos(i)
         enddo
      enddo
   enddo

   ! Surface tension forces
   #if phiflag == 1
      !$acc kernels
      !Obtain surface tension forces evaluated at the center of the cell (where phi is located)
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               if (ip .gt. nx) ip=1
               if (im .lt. 1) im=nx
               ! chempot
               chempot=phi(i,j,k)*(1.d0-phi(i,j,k))*(1.d0-2.d0*phi(i,j,k))*epsi-eps*(phi(ip,j,k)+phi(im,j,k)+phi(i,jp,k)+phi(i,jm,k)+phi(i,j,kp)+phi(i,j,km)- 6.d0*phi(i,j,k))*ddxi
               ! chempot*gradphi
               fxst(i,j,k)=6.d0*sigma*chempot*0.5d0*(phi(ip,j,k)-phi(im,j,k))*dxi
               fyst(i,j,k)=6.d0*sigma*chempot*0.5d0*(phi(i,jp,k)-phi(i,jm,k))*dxi
               fzst(i,j,k)=6.d0*sigma*chempot*0.5d0*(phi(i,j,kp)-phi(i,j,km))*dxi
            enddo
         enddo
      enddo
      !$acc end kernels

      ! Update halo of fxst, fyst and fzst (required then to interpolate at velocity points)
      !$acc host_data use_device(fxst)
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fxst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fxst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      !$acc end host_data 
      !$acc host_data use_device(fyst)
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fyst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fyst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      !$acc end host_data 
      !$acc host_data use_device(fzst)
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fzst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fzst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      !$acc end host_data 
      
      ! Interpolate force at velocity points
      !$acc kernels
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               im=i-1
               jm=j-1
               km=k-1
               if (im .lt. 1) im=nx
               rhsu(i,j,k)=rhsu(i,j,k) + 0.5d0*(fxst(im,j,k)+fxst(i,j,k))*rhoi
               rhsv(i,j,k)=rhsv(i,j,k) + 0.5d0*(fyst(i,jm,k)+fyst(i,j,k))*rhoi
               rhsw(i,j,k)=rhsw(i,j,k) + 0.5d0*(fzst(i,j,km)+fzst(i,j,k))*rhoi
               u(i,j,k) = u(i,j,k) + dt*(alpha*rhsu(i,j,k)-beta*rhsu_o(i,j,k))
               v(i,j,k) = v(i,j,k) + dt*(alpha*rhsv(i,j,k)-beta*rhsv_o(i,j,k))
               w(i,j,k) = w(i,j,k) + dt*(alpha*rhsw(i,j,k)-beta*rhsw_o(i,j,k))
               rhsu_o(i,j,k)=rhsu(i,j,k)
               rhsv_o(i,j,k)=rhsv(i,j,k)
               rhsw_o(i,j,k)=rhsw(i,j,k)
            enddo
         enddo
      enddo
      !$acc end kernels
   #else
      ! 6.2 find u, v and w star only in the inner nodes 
      !$acc kernels
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               u(i,j,k) = u(i,j,k) + dt*(alpha*rhsu(i,j,k)-beta*rhsu_o(i,j,k))
               v(i,j,k) = v(i,j,k) + dt*(alpha*rhsv(i,j,k)-beta*rhsv_o(i,j,k))
               w(i,j,k) = w(i,j,k) + dt*(alpha*rhsw(i,j,k)-beta*rhsw_o(i,j,k))
               rhsu_o(i,j,k)=rhsu(i,j,k)
               rhsv_o(i,j,k)=rhsv(i,j,k)
               rhsw_o(i,j,k)=rhsw(i,j,k)
            enddo
         enddo
      enddo
      !$acc end kernels
   #endif

   ! store rhs* in rhs*_o 
   ! First step is done with Euler explicit and then move to AB2 
   alpha=1.5d0
   beta= 0.5d0

   ! 5.3 update halos (y and z directions), required to then compute the RHS of Poisson equation because of staggered grid
   !$acc host_data use_device(u)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(v)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(w)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   ! (uncomment for profiling)
   ! call nvtxEndRange
   !########################################################################################################################################
   ! END STEP 6: USTAR COMPUTATION 
   !########################################################################################################################################


   
   !########################################################################################################################################
   ! START STEP 7: POISSON SOLVER FOR PRESSURE
   !########################################################################################################################################
   ! initialize rhs and analytical solution
   ! 7.1 Compute rhs of Poisson equation div*ustar: divergence at the cell center 
   ! (uncomment for profiling)
   ! call nvtxStartRange("Poisson")
   ! call nvtxStartRange("compute RHS")
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
               ip=i+1
               jp=j+1
               kp=k+1
               if (ip > nx) ip=1
               rhsp(i,j,k) =               (rho*dxi/dt)*(u(ip,j,k)-u(i,j,k))
               rhsp(i,j,k) = rhsp(i,j,k) + (rho*dxi/dt)*(v(i,jp,k)-v(i,j,k))
               rhsp(i,j,k) = rhsp(i,j,k) + (rho*dxi/dt)*(w(i,j,kp)-w(i,j,k))
               rhspp(i,j,k)=real(rhsp(i,j,k),kind=4)
         enddo
      enddo
   enddo
   !$acc end kernels
   ! call nvtxEndRange

   ! call nvtxStartRange("FFT forward w/ transpositions")
   !$acc host_data use_device(rhsp)
   status = cufftExecR2C(planXf, rhspp, psi_d)
   if (status /= CUFFT_SUCCESS) write(*,*) 'X forward error: ', status
   !$acc end host_data

   ! psi(kx,y,z) -> psi(y,z,kx)
   CHECK_CUDECOMP_EXIT(cudecompTransposeXToY(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_FLOAT_COMPLEX,piX_d2z%halo_extents, [0,0,0]))
   ! psi(y,z,kx) -> psi(ky,z,kx)
   status = cufftExecC2C(planY, psi_d, psi_d, CUFFT_FORWARD)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Y forward error: ', status
   ! psi(ky,z,kx) -> psi(z,kx,ky)
   CHECK_CUDECOMP_EXIT(cudecompTransposeYToZ(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_FLOAT_COMPLEX)) 
   ! psi(z,kx,ky) -> psi(kz,kx,ky)
   status = cufftExecC2C(planZ, psi_d, psi_d, CUFFT_FORWARD)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Z forward error: ', status
   ! END of FFT3D forward

   ! call nvtxEndRange
   np(piZ_d2z%order(1)) = piZ_d2z%shape(1)
   np(piZ_d2z%order(2)) = piZ_d2z%shape(2)
   np(piZ_d2z%order(3)) = piZ_d2z%shape(3)
   call c_f_pointer(c_devloc(psi_d), phi3d, piZ_d2z%shape)
   ! divide by -K**2, and normalize
   offsets(piZ_d2z%order(1)) = piZ_d2z%lo(1) - 1
   offsets(piZ_d2z%order(2)) = piZ_d2z%lo(2) - 1
   offsets(piZ_d2z%order(3)) = piZ_d2z%lo(3) - 1
   xoff = offsets(1)
   yoff = offsets(2)
   npx = np(1)
   npy = np(2)
   ! call nvtxStartRange("Solution")
   !$acc kernels
   do jl = 1, npy
      jg = yoff + jl
      do il = 1, npx
         ig = xoff + il
         do k = 1, nz
            k2 = kx_d(ig)**2 + kx_d(jg)**2 + kx_d(k)**2    
            phi3d(k,il,jl) = -phi3d(k,il,jl)/k2/(int(nx,8)*int(ny,8)*int(nz,8))
         enddo
      enddo
   enddo
   ! specify mean (corrects division by zero wavenumber above)
   if (xoff == 0 .and. yoff == 0) phi3d(1,1,1) = 0.0
   !$acc end kernels
   ! call nvtxEndRange
   ! call nvtxStartRange("FFT backwards w/ transpositions")

   ! psi(kz,kx,ky) -> psi(z,kx,ky)
   status = cufftExecC2C(planZ, psi_d, psi_d, CUFFT_INVERSE)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Z inverse error: ', status
   ! psi(z,kx,ky) -> psi(ky,z,kx)
   CHECK_CUDECOMP_EXIT(cudecompTransposeZToY(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_FLOAT_COMPLEX))
   ! psi(ky,z,kx) -> psi(y,z,kx)
   status = cufftExecC2C(planY, psi_d, psi_d, CUFFT_INVERSE)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Y inverse error: ', status
   ! psi(y,z,kx) -> psi(kx,y,z)
   CHECK_CUDECOMP_EXIT(cudecompTransposeYToX(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_FLOAT_COMPLEX,[0,0,0], piX_d2z%halo_extents))
   !$acc host_data use_device(p)
   ! psi(kx,y,z) -> psi(x,y,z)
   status = cufftExecC2R(planXb, psi_d, pp)
   if (status /= CUFFT_SUCCESS) write(*,*) 'X inverse error: ', status
   !$acc end host_data
   ! call nvtxEndRange


   !$acc host_data use_device(p)
   ! update halo nodes with pressure (needed for the pressure correction step), using device variable no need to use host-data
   ! Update X-pencil halos in Y and Z direction
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, pp, work_halo_d, CUDECOMP_FLOAT, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, pp, work_halo_d, CUDECOMP_FLOAT, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   ! call nvtxEndRange
   !########################################################################################################################################
   ! END STEP 7: POISSON SOLVER FOR PRESSURE
   !########################################################################################################################################

 

   ! (uncomment for profiling

   !########################################################################################################################################
   ! START STEP 8: VELOCITY CORRECTION
   ! ########################################################################################################################################
   ! 8.1 Correct velocity 
   ! 8.2 Remove mean velocity if using ABC forcing
   ! 8.3 Call halo exchnages along Y and Z for u,v,w
   ! Correct velocity, pressure has also the halo
   ! call nvtxStartRange("Correction")

   !$acc kernels 
   umean=0.d0
   vmean=0.d0
   wmean=0.d0
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i = 1, piX%shape(1) ! equal to nx (no halo on x)
            im=i-1
            jm=j-1
            km=k-1
            if (im < 1) im=nx
            u(i,j,k)=u(i,j,k) - dt/rho*(real(p(i,j,k),kind=8)-real(p(im,j,k),kind=8))*dxi
            v(i,j,k)=v(i,j,k) - dt/rho*(real(p(i,j,k),kind=8)-real(p(i,jm,k),kind=8))*dxi
            w(i,j,k)=w(i,j,k) - dt/rho*(real(p(i,j,k),kind=8)-real(p(i,j,km),kind=8))*dxi
            umean=umean + u(i,j,k)
            vmean=vmean + v(i,j,k)
            wmean=wmean + w(i,j,k)
          enddo
      enddo
   enddo
   !$acc end kernels 

   ! Remove mean velocity (get local mean of the rank)

   ! Divide by total number of points in the pencil
   umean=umean/nx/(piX%shape(2)-2*halo_ext)/(piX%shape(3)-2*halo_ext)
   vmean=vmean/nx/(piX%shape(2)-2*halo_ext)/(piX%shape(3)-2*halo_ext)
   wmean=wmean/nx/(piX%shape(2)-2*halo_ext)/(piX%shape(3)-2*halo_ext)

   ! Find global mean (MPI_SUM and then divide by number of ranks)
   call MPI_Allreduce(umean,gumean,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
   call MPI_Allreduce(vmean,gvmean,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
   call MPI_Allreduce(wmean,gwmean,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)

   ! remove mean value
   !$acc kernels 
   u=u-(gumean/ranks)
   v=v-(gvmean/ranks)
   w=w-(gwmean/ranks)
   !$acc end kernels 

   ! 8.3 update halos (y and z directions), required to then compute the RHS of Poisson equation because of staggered grid
   !$acc host_data use_device(u)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(v)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(w)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 

   ! find local maximum velocity
   uc=maxval(u)
   vc=maxval(v)
   wc=maxval(w)
   umax=max(wc,max(uc,vc))
   call MPI_Allreduce(umax,gumax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD, ierr)
   cou=gumax*dt*dxi
   if (rank.eq.0) then
      write(*,*) "CFL (max among tasks)", cou
      if (cou .gt. 7) stop
   endif

   call cpu_time(timef)
   if (rank.eq.0) print '(" Time elapsed = ",f8.1," ms")',1000*(timef-times)

   ! (uncomment for profiling)
   ! call nvtxEndRange

   !########################################################################################################################################
   ! END STEP 8: VELOCITY CORRECTION  
   !########################################################################################################################################


   !########################################################################################################################################
   ! START STEP 9: OUTPUT FIELDS 
   ! ########################################################################################################################################
   if (mod(t,dump) .eq. 0) then
      if (rank .eq. 0) write(*,*) "Saving output files"
      ! write velocity and pressure fiels (1-4)
      call writefield(t,1)
      call writefield(t,2)
      call writefield(t,3)
      call writefield(t,4)
      #if phiflag == 1
         ! write phase-field (5)
         call writefield(t,5)
      #endif
   endif
   !########################################################################################################################################
   ! END STEP 9: OUTPUT FIELDS N  
   !########################################################################################################################################

! (uncomment for profiling)
! call nvtxEndRange

enddo
call cpu_time(t_end)
elapsed = t_end-t_start
if (rank .eq. 0) write(*,*)  'Elapsed time (seconds):', elapsed
#if partflag == 1
!$acc end data
!$acc end data
!$acc end data
#endif
!$acc end data
!$acc end data
!$acc end data

#if partflag == 1
! Particle variables
deallocate(part)
deallocate(partbuff)
deallocate(vec_p)
deallocate(order_p)
deallocate(buffvar1)
deallocate(buffvar2)
#endif

! Remove allocated variables (add new)
deallocate(x_ext)
deallocate(u,v,w)
deallocate(tanh_psi, mysin, mycos)
deallocate(rhsu,rhsv,rhsw)
deallocate(rhsu_o,rhsv_o,rhsw_o)
deallocate(phi,rhsphi,rhsphi_o,normx,normy,normz)

call mpi_finalize(ierr)

end program main

