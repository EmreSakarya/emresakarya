!=======================================================================
!  C5G7  --  2-D, 7-group neutron TRANSPORT benchmark (OECD/NEA)
!            Discrete-ordinates (S_N) eigenvalue solver
!
!  Geometry  : 3x3 assemblies (2x2 UO2/MOX fuel block + L-shaped water
!              reflector), 17x17 pins/assembly, pin pitch 1.26 cm,
!              cylinder radius 0.54 cm.  Domain 64.26 x 64.26 cm.
!  Boundary  : West & North  = reflective (core symmetry planes)
!              East & South  = vacuum     (outer reflector edge)
!  Physics   : 7-group transport-corrected (P0) cross sections,
!              up-scattering in thermal groups, fission.
!  Method    : product angular quadrature (uniform azimuth x Gauss-
!              Legendre polar), diamond-difference spatial sweep with
!              negative-flux fix-up, power (eigenvalue) iteration,
!              Gauss-Seidel scattering source, OpenMP over directions.
!
!  Reference k_eff (MCNP, NEA/NSC/DOC(2003)16) = 1.18655
!
!  Build : gfortran -O3 -fopenmp c5g7_sn.f90 -o c5g7
!  Run   : ./c5g7 [cells_per_pin] [n_azimuthal] [n_polar]
!          e.g.  ./c5g7 8 16 3
!=======================================================================
module prec
  implicit none
  integer, parameter :: dp = kind(1.0d0)
end module prec

module data
  use prec
  implicit none
  integer, parameter :: ng = 7          ! energy groups
  integer, parameter :: nmat = 7        ! materials
  integer, parameter :: npa = 17        ! pins per assembly side
  integer, parameter :: nac = 3         ! assemblies per core side
  real(dp), parameter :: pitch = 1.26_dp, rfuel = 0.54_dp
  real(dp), parameter :: pi = 3.141592653589793_dp
  real(dp), parameter :: kref = 1.18655_dp

  ! cross sections
  real(dp) :: st(ng,nmat), nsf(ng,nmat), chi(ng,nmat), ss(ng,ng,nmat)
  ! assembly pin-material maps and core layout
  integer  :: amap_uo2(npa,npa), amap_mox(npa,npa)
  integer  :: core(nac,nac)             ! 0=reflector 1=UO2 2=MOX  (core(ai,aj))

  ! mesh
  integer  :: nx, ny, cpp
  real(dp) :: h
  integer,  allocatable :: matid(:,:)   ! material of fine cell (i,j)
  real(dp), allocatable :: phi(:,:,:)   ! scalar flux (i,j,g)
  real(dp), allocatable :: fd(:,:)      ! fission density
  real(dp), allocatable :: qsrc(:,:)    ! angular source for current group
  integer  :: gnow                      ! current group being swept

  ! quadrature
  integer  :: nazi, npol, ndir
  real(dp), allocatable :: omx(:), omy(:), wt(:)
  integer,  allocatable :: mirx(:), miry(:)

  ! reflective boundary angular fluxes, per group (lagged, double-buffered)
  real(dp), allocatable :: bWE_o(:,:,:), bWE_n(:,:,:)   ! west face  (j,m,g)
  real(dp), allocatable :: bNO_o(:,:,:), bNO_n(:,:,:)   ! north face (i,m,g)
end module data

!-----------------------------------------------------------------------
program c5g7_sn
  use prec; use data
  !$ use omp_lib
  implicit none
  integer  :: i,j,g,gp,m, outit, init, narg, ios
  integer  :: nin                       ! inner scattering passes / outer
  integer  :: maxout, npass
  real(dp) :: keff, knew, psum_o, psum_n
  real(dp) :: res, fnorm, t0, t1, tol
  character(len=32) :: arg
  real(dp), allocatable :: buf(:,:), fdold(:,:)

  ! ---- parameters (defaults, overridable on command line) ----
  ! usage: ./c5g7 [cells_per_pin] [n_azimuthal] [n_polar] [inner_passes]
  cpp = 8;  nazi = 16;  npol = 3;  nin = 3;  maxout = 800;  tol = 1.0e-5_dp
  narg = command_argument_count()
  if (narg>=1) then; call get_command_argument(1,arg); read(arg,*,iostat=ios) cpp;  endif
  if (narg>=2) then; call get_command_argument(2,arg); read(arg,*,iostat=ios) nazi; endif
  if (narg>=3) then; call get_command_argument(3,arg); read(arg,*,iostat=ios) npol; endif
  if (narg>=4) then; call get_command_argument(4,arg); read(arg,*,iostat=ios) nin;  endif
  if (mod(nazi,2)/=0) nazi = nazi+1     ! need even azimuth for reflective mirrors

  nx = nac*npa*cpp;  ny = nx;  h = pitch/real(cpp,dp)

  call set_xsec()
  call set_maps()
  call build_quad()

  write(*,'(a)') '======================================================'
  write(*,'(a)') '  C5G7  2-D 7-group S_N transport eigenvalue solver'
  write(*,'(a)') '======================================================'
  write(*,'(a,i5,a,i5,a,f7.4,a)') '  mesh        : ',nx,' x ',ny,'  (h =',h,' cm)'
  write(*,'(a,i9)')               '  fine cells  : ',nx*ny
  write(*,'(a,i3,a,i3,a,i5,a)')   '  quadrature  : ',nazi,' azim x ',npol, &
                                  ' polar = ',ndir,' directions/hemisphere'
  write(*,'(a,i3,a,i4)')          '  inner/outer : ',nin,'   max outer : ',maxout
  !$ write(*,'(a,i3)')            '  OpenMP thr. : ',omp_get_max_threads()
  write(*,'(a)') '------------------------------------------------------'

  ! ---- allocate fields ----
  allocate(matid(nx,ny), phi(nx,ny,ng), fd(nx,ny), qsrc(nx,ny))
  allocate(buf(nx,ny), fdold(nx,ny))
  allocate(bWE_o(ny,ndir,ng), bWE_n(ny,ndir,ng), bNO_o(nx,ndir,ng), bNO_n(nx,ndir,ng))
  call assign_materials()

  phi = 1.0_dp
  fd  = 0.0_dp
  do j=1,ny; do i=1,nx
     do g=1,ng; fd(i,j) = fd(i,j) + nsf(g,matid(i,j))*phi(i,j,g); enddo
  enddo; enddo
  fd = fd/sum(fd)                       ! normalise fission source (integral = 1)
  keff = 1.15_dp                        ! warm start near the expected eigenvalue
  bWE_o = 0.0_dp; bNO_o = 0.0_dp

  t0 = 0.0_dp
  !$ t0 = omp_get_wtime()
  ! ================= POWER (EIGENVALUE) ITERATION =================
  do outit = 1, maxout
     psum_o = sum(fd)
     fdold  = fd
     ! ---- inner: relax scattering (incl. up-scatter) for fixed fission source ----
     ! extra passes on the first outer settle the flux shape (kills startup spike)
     npass = nin
     if (outit==1) npass = max(nin, 25)
     do init = 1, npass
        do g = 1, ng
           gnow = g
           ! isotropic source  q = (1/4pi)[ sum_g' Ss(g<-g') phi_g' + chi_g Fd/keff ]
           !$omp parallel do collapse(2) private(i,j,gp) schedule(static)
           do j=1,ny
              do i=1,nx
                 qsrc(i,j) = chi(g,matid(i,j))*fd(i,j)/keff
                 do gp=1,ng
                    qsrc(i,j) = qsrc(i,j) + ss(g,gp,matid(i,j))*phi(i,j,gp)
                 enddo
                 qsrc(i,j) = qsrc(i,j)/(4.0_dp*pi)
              enddo
           enddo
           !$omp end parallel do

           ! transport sweep over all directions  (phi_g = sum_m w_m psi_m)
           bWE_n(:,:,g) = 0.0_dp;  bNO_n(:,:,g) = 0.0_dp
           phi(:,:,g) = 0.0_dp
           buf = 0.0_dp
           !$omp parallel private(m) firstprivate(buf)
           !$omp do schedule(dynamic)
           do m = 1, ndir
              call sweep_dir(m, buf)
           enddo
           !$omp end do
           !$omp critical
           phi(:,:,g) = phi(:,:,g) + buf
           !$omp end critical
           !$omp end parallel
           bWE_o(:,:,g) = bWE_n(:,:,g);  bNO_o(:,:,g) = bNO_n(:,:,g)   ! lagged update
        enddo
     enddo

     ! ---- eigenvalue + fission-source update ----
     !$omp parallel do collapse(2) private(i,j,g) schedule(static)
     do j=1,ny; do i=1,nx
        fd(i,j) = 0.0_dp
        do g=1,ng; fd(i,j) = fd(i,j) + nsf(g,matid(i,j))*phi(i,j,g); enddo
     enddo; enddo
     !$omp end parallel do
     psum_n = sum(fd)
     knew = keff * psum_n/psum_o

     ! renormalise (keep integral fission source = 1) and measure residual
     fd  = fd  / psum_n
     phi = phi / psum_n
     bWE_o = bWE_o / psum_n;  bNO_o = bNO_o / psum_n
     res = sum(abs(fd-fdold))
     fnorm = abs(knew-keff)
     keff = knew

     if (mod(outit,10)==0 .or. outit<=5) &
        write(*,'(a,i4,a,f10.6,a,es10.3,a,es10.3)') &
        '  outer ',outit,'   k = ',knew,'   |dk| = ',fnorm,'   |dFd| = ',res
     if (res < tol .and. fnorm < tol .and. outit>15) exit
  enddo
  t1 = 0.0_dp
  !$ t1 = omp_get_wtime()

  write(*,'(a)') '------------------------------------------------------'
  write(*,'(a,f10.6)') '  C5G7  k_eff  (this solver) = ', keff
  write(*,'(a,f10.6)') '  C5G7  k_eff  (NEA reference)= ', kref
  write(*,'(a,f9.1,a)')'  difference  = ', (keff-kref)*1.0e5_dp, ' pcm'
  write(*,'(a,f9.1,a)')'  wall time   = ', t1-t0, ' s'
  write(*,'(a)') '======================================================'

  call write_flux_csv()
  call write_pin_powers()
end program c5g7_sn

!-----------------------------------------------------------------------
! One discrete-ordinates sweep for direction m; accumulate w*psi into buf
! and record outgoing reflective-boundary angular fluxes.
subroutine sweep_dir(m, buf)
  use prec; use data
  implicit none
  integer, intent(in) :: m
  real(dp), intent(inout) :: buf(nx,ny)
  integer  :: i,j, i0,i1,di, j0,j1,dj
  real(dp) :: ox,oy,w, cx,cy, xin,yin,psi,xout,yout, stot, denom
  real(dp) :: yedge(nx)

  ox = omx(m);  oy = omy(m);  w = wt(m)
  cx = 2.0_dp*abs(ox)/h;  cy = 2.0_dp*abs(oy)/h

  if (ox > 0.0_dp) then; i0=1; i1=nx; di=1; else; i0=nx; i1=1; di=-1; endif
  if (oy > 0.0_dp) then; j0=1; j1=ny; dj=1; else; j0=ny; j1=1; dj=-1; endif

  ! incoming y-face for the first row
  if (oy > 0.0_dp) then
     yedge = 0.0_dp                      ! south = vacuum
  else
     yedge = bNO_o(:,miry(m),gnow)       ! north = reflective
  endif

  j = j0
  do
     ! incoming x for this row
     if (ox > 0.0_dp) then
        xin = bWE_o(j, mirx(m), gnow)    ! west = reflective
     else
        xin = 0.0_dp                     ! east = vacuum
     endif

     i = i0
     do
        yin = yedge(i)
        stot = st(gnow, matid(i,j))
        denom = cx + cy + stot
        psi = (qsrc(i,j) + cx*xin + cy*yin) / denom
        if (psi < 0.0_dp) psi = 0.0_dp
        xout = 2.0_dp*psi - xin;  if (xout < 0.0_dp) xout = 0.0_dp
        yout = 2.0_dp*psi - yin;  if (yout < 0.0_dp) yout = 0.0_dp
        buf(i,j) = buf(i,j) + w*psi
        ! store reflective outgoing
        if (ox < 0.0_dp .and. i==1)  bWE_n(j,m,gnow) = xout   ! west outgoing
        if (oy > 0.0_dp .and. j==ny) bNO_n(i,m,gnow) = yout   ! north outgoing
        xin = xout
        yedge(i) = yout
        if (i==i1) exit
        i = i + di
     enddo
     if (j==j1) exit
     j = j + dj
  enddo
end subroutine sweep_dir

!-----------------------------------------------------------------------
subroutine set_xsec()
  use prec; use data
  implicit none
  ! === C5G7 macroscopic cross sections (transport-corrected, P0) ===
  ! materials: 1 UO2  2 MOX4.3  3 MOX7.0  4 MOX8.7  5 FissChamber  6 GuideTube  7 Moderator
  ss = 0.0_dp
  st(:,1)  = (/ 1.779490E-01_dp, 3.298050E-01_dp, 4.803880E-01_dp, 5.543670E-01_dp, &
                 3.118010E-01_dp, 3.951680E-01_dp, 5.644060E-01_dp /)
  nsf(:,1) = (/ 2.005998E-02_dp, 2.027303E-03_dp, 1.570599E-02_dp, 4.518301E-02_dp, &
                 4.334208E-02_dp, 2.020901E-01_dp, 5.257105E-01_dp /)
  chi(:,1) = (/ 5.878190E-01_dp, 4.117600E-01_dp, 3.390600E-04_dp, 1.176100E-07_dp, &
                 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp /)
  ss(1,1,1) = 1.275370E-01_dp
  ss(2,1,1) = 4.237800E-02_dp
  ss(2,2,1) = 3.244560E-01_dp
  ss(3,1,1) = 9.437400E-06_dp
  ss(3,2,1) = 1.631400E-03_dp
  ss(3,3,1) = 4.509400E-01_dp
  ss(4,1,1) = 5.516300E-09_dp
  ss(4,2,1) = 3.142700E-09_dp
  ss(4,3,1) = 2.679200E-03_dp
  ss(4,4,1) = 4.525650E-01_dp
  ss(4,5,1) = 1.252500E-04_dp
  ss(5,4,1) = 5.566400E-03_dp
  ss(5,5,1) = 2.714010E-01_dp
  ss(5,6,1) = 1.296800E-03_dp
  ss(6,5,1) = 1.025500E-02_dp
  ss(6,6,1) = 2.658020E-01_dp
  ss(6,7,1) = 8.545800E-03_dp
  ss(7,5,1) = 1.002100E-08_dp
  ss(7,6,1) = 1.680900E-02_dp
  ss(7,7,1) = 2.730800E-01_dp
  st(:,2)  = (/ 1.787310E-01_dp, 3.308490E-01_dp, 4.837720E-01_dp, 5.669220E-01_dp, &
                 4.262270E-01_dp, 6.789970E-01_dp, 6.828520E-01_dp /)
  nsf(:,2) = (/ 2.175300E-02_dp, 2.535103E-03_dp, 1.626799E-02_dp, 6.547410E-02_dp, &
                 3.072409E-02_dp, 6.666510E-01_dp, 7.139904E-01_dp /)
  chi(:,2) = (/ 5.878190E-01_dp, 4.117600E-01_dp, 3.390600E-04_dp, 1.176100E-07_dp, &
                 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp /)
  ss(1,1,2) = 1.288760E-01_dp
  ss(2,1,2) = 4.141300E-02_dp
  ss(2,2,2) = 3.254520E-01_dp
  ss(3,1,2) = 8.229000E-06_dp
  ss(3,2,2) = 1.639500E-03_dp
  ss(3,3,2) = 4.531880E-01_dp
  ss(4,1,2) = 5.040500E-09_dp
  ss(4,2,2) = 1.598200E-09_dp
  ss(4,3,2) = 2.614200E-03_dp
  ss(4,4,2) = 4.571730E-01_dp
  ss(4,5,2) = 1.604600E-04_dp
  ss(5,4,2) = 5.539400E-03_dp
  ss(5,5,2) = 2.768140E-01_dp
  ss(5,6,2) = 2.005100E-03_dp
  ss(6,5,2) = 9.312700E-03_dp
  ss(6,6,2) = 2.529620E-01_dp
  ss(6,7,2) = 8.494800E-03_dp
  ss(7,5,2) = 9.165600E-09_dp
  ss(7,6,2) = 1.485000E-02_dp
  ss(7,7,2) = 2.650070E-01_dp
  st(:,3)  = (/ 1.813230E-01_dp, 3.343680E-01_dp, 4.937850E-01_dp, 5.912160E-01_dp, &
                 4.741980E-01_dp, 8.336010E-01_dp, 8.536030E-01_dp /)
  nsf(:,3) = (/ 2.381395E-02_dp, 3.858689E-03_dp, 2.413400E-02_dp, 9.436622E-02_dp, &
                 4.576988E-02_dp, 9.281814E-01_dp, 1.043200E+00_dp /)
  chi(:,3) = (/ 5.878190E-01_dp, 4.117600E-01_dp, 3.390600E-04_dp, 1.176100E-07_dp, &
                 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp /)
  ss(1,1,3) = 1.304570E-01_dp
  ss(2,1,3) = 4.179200E-02_dp
  ss(2,2,3) = 3.284280E-01_dp
  ss(3,1,3) = 8.510500E-06_dp
  ss(3,2,3) = 1.643600E-03_dp
  ss(3,3,3) = 4.583710E-01_dp
  ss(4,1,3) = 5.132900E-09_dp
  ss(4,2,3) = 2.201700E-09_dp
  ss(4,3,3) = 2.533100E-03_dp
  ss(4,4,3) = 4.637090E-01_dp
  ss(4,5,3) = 1.761900E-04_dp
  ss(5,4,3) = 5.476600E-03_dp
  ss(5,5,3) = 2.823130E-01_dp
  ss(5,6,3) = 2.276000E-03_dp
  ss(6,5,3) = 8.728900E-03_dp
  ss(6,6,3) = 2.497510E-01_dp
  ss(6,7,3) = 8.864500E-03_dp
  ss(7,5,3) = 9.001600E-09_dp
  ss(7,6,3) = 1.311400E-02_dp
  ss(7,7,3) = 2.595290E-01_dp
  st(:,4)  = (/ 1.830450E-01_dp, 3.367050E-01_dp, 5.005070E-01_dp, 6.061740E-01_dp, &
                 5.027540E-01_dp, 9.210280E-01_dp, 9.552310E-01_dp /)
  nsf(:,4) = (/ 2.518600E-02_dp, 4.739509E-03_dp, 2.947805E-02_dp, 1.122500E-01_dp, &
                 5.530301E-02_dp, 1.074999E+00_dp, 1.239298E+00_dp /)
  chi(:,4) = (/ 5.878190E-01_dp, 4.117600E-01_dp, 3.390600E-04_dp, 1.176100E-07_dp, &
                 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp /)
  ss(1,1,4) = 1.315040E-01_dp
  ss(2,1,4) = 4.204600E-02_dp
  ss(2,2,4) = 3.304030E-01_dp
  ss(3,1,4) = 8.697200E-06_dp
  ss(3,2,4) = 1.646300E-03_dp
  ss(3,3,4) = 4.617920E-01_dp
  ss(4,1,4) = 5.193800E-09_dp
  ss(4,2,4) = 2.600600E-09_dp
  ss(4,3,4) = 2.474900E-03_dp
  ss(4,4,4) = 4.680210E-01_dp
  ss(4,5,4) = 1.859700E-04_dp
  ss(5,4,4) = 5.433000E-03_dp
  ss(5,5,4) = 2.857710E-01_dp
  ss(5,6,4) = 2.391600E-03_dp
  ss(6,5,4) = 8.397300E-03_dp
  ss(6,6,4) = 2.476140E-01_dp
  ss(6,7,4) = 8.968100E-03_dp
  ss(7,5,4) = 8.928000E-09_dp
  ss(7,6,4) = 1.232200E-02_dp
  ss(7,7,4) = 2.560930E-01_dp
  st(:,5)  = (/ 1.260320E-01_dp, 2.931600E-01_dp, 2.842500E-01_dp, 2.810200E-01_dp, &
                 3.344600E-01_dp, 5.656400E-01_dp, 1.172140E+00_dp /)
  nsf(:,5) = (/ 1.323401E-08_dp, 1.434500E-08_dp, 1.128599E-06_dp, 1.276299E-05_dp, &
                 3.538502E-07_dp, 1.740099E-06_dp, 5.063302E-06_dp /)
  chi(:,5) = (/ 5.878190E-01_dp, 4.117600E-01_dp, 3.390600E-04_dp, 1.176100E-07_dp, &
                 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp /)
  ss(1,1,5) = 6.616590E-02_dp
  ss(2,1,5) = 5.907000E-02_dp
  ss(2,2,5) = 2.403770E-01_dp
  ss(3,1,5) = 2.833400E-04_dp
  ss(3,2,5) = 5.243500E-02_dp
  ss(3,3,5) = 1.834250E-01_dp
  ss(4,1,5) = 1.462200E-06_dp
  ss(4,2,5) = 2.499000E-04_dp
  ss(4,3,5) = 9.228800E-02_dp
  ss(4,4,5) = 7.907690E-02_dp
  ss(4,5,5) = 3.734000E-05_dp
  ss(5,1,5) = 2.064200E-08_dp
  ss(5,2,5) = 1.923900E-05_dp
  ss(5,3,5) = 6.936500E-03_dp
  ss(5,4,5) = 1.699900E-01_dp
  ss(5,5,5) = 9.975700E-02_dp
  ss(5,6,5) = 9.174200E-04_dp
  ss(6,2,5) = 2.987500E-06_dp
  ss(6,3,5) = 1.079000E-03_dp
  ss(6,4,5) = 2.586000E-02_dp
  ss(6,5,5) = 2.067900E-01_dp
  ss(6,6,5) = 3.167740E-01_dp
  ss(6,7,5) = 4.979300E-02_dp
  ss(7,2,5) = 4.214000E-07_dp
  ss(7,3,5) = 2.054300E-04_dp
  ss(7,4,5) = 4.925600E-03_dp
  ss(7,5,5) = 2.447800E-02_dp
  ss(7,6,5) = 2.387600E-01_dp
  ss(7,7,5) = 1.099100E+00_dp
  st(:,6)  = (/ 1.260320E-01_dp, 2.931600E-01_dp, 2.842400E-01_dp, 2.809600E-01_dp, &
                 3.344400E-01_dp, 5.656400E-01_dp, 1.172150E+00_dp /)
  nsf(:,6) = (/ 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp, &
                 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp /)
  chi(:,6) = (/ 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp, &
                 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp /)
  ss(1,1,6) = 6.616590E-02_dp
  ss(2,1,6) = 5.907000E-02_dp
  ss(2,2,6) = 2.403770E-01_dp
  ss(3,1,6) = 2.833400E-04_dp
  ss(3,2,6) = 5.243500E-02_dp
  ss(3,3,6) = 1.832970E-01_dp
  ss(4,1,6) = 1.462200E-06_dp
  ss(4,2,6) = 2.499000E-04_dp
  ss(4,3,6) = 9.239700E-02_dp
  ss(4,4,6) = 7.885110E-02_dp
  ss(4,5,6) = 3.733300E-05_dp
  ss(5,1,6) = 2.064200E-08_dp
  ss(5,2,6) = 1.923900E-05_dp
  ss(5,3,6) = 6.944600E-03_dp
  ss(5,4,6) = 1.701400E-01_dp
  ss(5,5,6) = 9.973720E-02_dp
  ss(5,6,6) = 9.172600E-04_dp
  ss(6,2,6) = 2.987500E-06_dp
  ss(6,3,6) = 1.079000E-03_dp
  ss(6,4,6) = 2.586000E-02_dp
  ss(6,5,6) = 2.067900E-01_dp
  ss(6,6,6) = 3.167740E-01_dp
  ss(6,7,6) = 4.979300E-02_dp
  ss(7,2,6) = 4.214000E-07_dp
  ss(7,3,6) = 2.054300E-04_dp
  ss(7,4,6) = 4.925600E-03_dp
  ss(7,5,6) = 2.447800E-02_dp
  ss(7,6,6) = 2.387600E-01_dp
  ss(7,7,6) = 1.099100E+00_dp
  st(:,7)  = (/ 1.592060E-01_dp, 4.129700E-01_dp, 5.903100E-01_dp, 5.843500E-01_dp, &
                 7.180000E-01_dp, 1.254450E+00_dp, 2.650380E+00_dp /)
  nsf(:,7) = (/ 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp, &
                 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp /)
  chi(:,7) = (/ 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp, &
                 0.000000E+00_dp, 0.000000E+00_dp, 0.000000E+00_dp /)
  ss(1,1,7) = 4.447770E-02_dp
  ss(2,1,7) = 1.134000E-01_dp
  ss(2,2,7) = 2.823340E-01_dp
  ss(3,1,7) = 7.234700E-04_dp
  ss(3,2,7) = 1.299400E-01_dp
  ss(3,3,7) = 3.452560E-01_dp
  ss(4,1,7) = 3.749900E-06_dp
  ss(4,2,7) = 6.234000E-04_dp
  ss(4,3,7) = 2.245700E-01_dp
  ss(4,4,7) = 9.102840E-02_dp
  ss(4,5,7) = 7.143700E-05_dp
  ss(5,1,7) = 5.318400E-08_dp
  ss(5,2,7) = 4.800200E-05_dp
  ss(5,3,7) = 1.699900E-02_dp
  ss(5,4,7) = 4.155100E-01_dp
  ss(5,5,7) = 1.391380E-01_dp
  ss(5,6,7) = 2.215700E-03_dp
  ss(6,2,7) = 7.448600E-06_dp
  ss(6,3,7) = 2.644300E-03_dp
  ss(6,4,7) = 6.373200E-02_dp
  ss(6,5,7) = 5.118200E-01_dp
  ss(6,6,7) = 6.999130E-01_dp
  ss(6,7,7) = 1.324400E-01_dp
  ss(7,2,7) = 1.045500E-06_dp
  ss(7,3,7) = 5.034400E-04_dp
  ss(7,4,7) = 1.213900E-02_dp
  ss(7,5,7) = 6.122900E-02_dp
  ss(7,6,7) = 5.373200E-01_dp
  ss(7,7,7) = 2.480700E+00_dp
end subroutine set_xsec

!-----------------------------------------------------------------------
subroutine set_maps()
  use prec; use data
  implicit none
  ! 17x17 assembly pin-material maps (row 1 = first text row; pattern is symmetric)
  ! amap_uo2
  amap_uo2(1,:) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  amap_uo2(2,:) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  amap_uo2(3,:) = (/ 1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1 /)
  amap_uo2(4,:) = (/ 1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1 /)
  amap_uo2(5,:) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  amap_uo2(6,:) = (/ 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1 /)
  amap_uo2(7,:) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  amap_uo2(8,:) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  amap_uo2(9,:) = (/ 1, 1, 6, 1, 1, 6, 1, 1, 5, 1, 1, 6, 1, 1, 6, 1, 1 /)
  amap_uo2(10,:) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  amap_uo2(11,:) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  amap_uo2(12,:) = (/ 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1 /)
  amap_uo2(13,:) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  amap_uo2(14,:) = (/ 1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1 /)
  amap_uo2(15,:) = (/ 1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1 /)
  amap_uo2(16,:) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  amap_uo2(17,:) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  ! amap_mox
  amap_mox(1,:) = (/ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 /)
  amap_mox(2,:) = (/ 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2 /)
  amap_mox(3,:) = (/ 2, 3, 3, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 3, 3, 2 /)
  amap_mox(4,:) = (/ 2, 3, 3, 6, 3, 4, 4, 4, 4, 4, 4, 4, 3, 6, 3, 3, 2 /)
  amap_mox(5,:) = (/ 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2 /)
  amap_mox(6,:) = (/ 2, 3, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 3, 2 /)
  amap_mox(7,:) = (/ 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2 /)
  amap_mox(8,:) = (/ 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2 /)
  amap_mox(9,:) = (/ 2, 3, 6, 4, 4, 6, 4, 4, 5, 4, 4, 6, 4, 4, 6, 3, 2 /)
  amap_mox(10,:) = (/ 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2 /)
  amap_mox(11,:) = (/ 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2 /)
  amap_mox(12,:) = (/ 2, 3, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 3, 2 /)
  amap_mox(13,:) = (/ 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2 /)
  amap_mox(14,:) = (/ 2, 3, 3, 6, 3, 4, 4, 4, 4, 4, 4, 4, 3, 6, 3, 3, 2 /)
  amap_mox(15,:) = (/ 2, 3, 3, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 3, 3, 2 /)
  amap_mox(16,:) = (/ 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2 /)
  amap_mox(17,:) = (/ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 /)
  ! core layout (core(ai,aj): ai=1 West..3 East ; aj=1 South..3 North)
  ! OpenMOC root:  top  row = UO2 MOX REFL ;  mid = MOX UO2 REFL ; bot = REFL*3
  core = 0
  core(1,3)=1; core(2,3)=2; core(3,3)=0     ! North row
  core(1,2)=2; core(2,2)=1; core(3,2)=0     ! Middle row
  core(1,1)=0; core(2,1)=0; core(3,1)=0     ! South row (all reflector)
end subroutine set_maps

!-----------------------------------------------------------------------
subroutine build_mesh()
  use prec; use data
  implicit none
  ! nothing (mesh sizes set in main); placeholder for clarity
end subroutine build_mesh

!-----------------------------------------------------------------------
! Assign a material to every fine cell.  Volume-preserving "digital disk":
! within each pin cell the N_target cells closest to the pin centre are
! made fuel (N_target = round(pi r^2 / h^2)), the rest moderator.  This
! keeps the fuel volume exact (and the pin symmetric) at any mesh size.
subroutine assign_materials()
  use prec; use data
  implicit none
  integer  :: pc,pr, gi,gj, ntar
  real(dp) :: vfuel, vtot
  ntar = nint(pi*rfuel*rfuel/(h*h))                 ! fuel cells per pin
  !$omp parallel do collapse(2) schedule(static)
  do pr = 1, nac*npa
     do pc = 1, nac*npa
        call assign_pin(pc, pr, ntar)
     enddo
  enddo
  !$omp end parallel do
  ! report realised fuel volume fraction
  vfuel = 0.0_dp; vtot = real(nx,dp)*real(ny,dp)
  do gj=1,ny; do gi=1,nx
     if (matid(gi,gj)>=1 .and. matid(gi,gj)<=4) vfuel = vfuel + 1.0_dp
  enddo; enddo
  write(*,'(a,f7.4,a,i5,a)') '  fuel vol fr.: ',vfuel/vtot,'   (',ntar,' fuel cells/pin)'
end subroutine assign_materials

! Assign all fine cells of one pin (pc,pr). Local automatic arrays make this
! thread-safe inside the OpenMP loop.
subroutine assign_pin(pc, pr, ntar)
  use prec; use data
  implicit none
  integer, intent(in) :: pc, pr, ntar
  integer  :: ai,aj, lc,lr, mfuel, atype, li,lj, gi,gj, k
  real(dp) :: pcx,pcy, xc,yc, thr, d2(cpp*cpp)

  ai = (pc-1)/npa + 1;   aj = (pr-1)/npa + 1
  lc = pc - npa*(ai-1);  lr = pr - npa*(aj-1)
  atype = core(ai,aj)
  if (atype==0) then
     mfuel = -1
  elseif (atype==1) then
     mfuel = amap_uo2(lr,lc)
  else
     mfuel = amap_mox(lr,lc)
  endif
  if (mfuel < 0) then                                ! reflector pin = all water
     do lj=1,cpp; do li=1,cpp
        matid((pc-1)*cpp+li,(pr-1)*cpp+lj) = 7
     enddo; enddo
     return
  endif
  pcx = (real(pc,dp)-0.5_dp)*pitch
  pcy = (real(pr,dp)-0.5_dp)*pitch
  k = 0
  do lj=1,cpp; do li=1,cpp
     gi = (pc-1)*cpp+li;  gj = (pr-1)*cpp+lj
     xc = (real(gi,dp)-0.5_dp)*h;  yc = (real(gj,dp)-0.5_dp)*h
     k = k+1;  d2(k) = (xc-pcx)**2 + (yc-pcy)**2
  enddo; enddo
  call kth_smallest(d2, cpp*cpp, max(1,min(ntar,cpp*cpp)), thr)
  k = 0
  do lj=1,cpp; do li=1,cpp
     gi = (pc-1)*cpp+li;  gj = (pr-1)*cpp+lj
     k = k+1
     if (d2(k) <= thr) then
        matid(gi,gj) = mfuel
     else
        matid(gi,gj) = 7
     endif
  enddo; enddo
end subroutine assign_pin

! threshold = the k-th smallest value of a(1:n) (copy, partial selection sort)
subroutine kth_smallest(a, n, kk, thr)
  use prec
  implicit none
  integer, intent(in) :: n, kk
  real(dp), intent(in) :: a(n)
  real(dp), intent(out) :: thr
  real(dp) :: b(n), t
  integer :: i,j,mn
  b = a
  do i = 1, kk
     mn = i
     do j = i+1, n
        if (b(j) < b(mn)) mn = j
     enddo
     t = b(i); b(i) = b(mn); b(mn) = t
  enddo
  thr = b(kk)
end subroutine kth_smallest

!-----------------------------------------------------------------------
! Product angular quadrature for the xi>0 hemisphere (doubled weight for
! the xi<0 mirror).  Uniform azimuth, Gauss-Legendre polar in xi=cos(theta).
! Weights sum to 4*pi.  Builds x- and y-reflection mirror indices.
subroutine build_quad()
  use prec; use data
  implicit none
  integer  :: a,p,m, mm
  real(dp) :: dom, om, xi, sth, wpol(8), xipol(8), wazi
  real(dp) :: tox,toy, best
  integer  :: bi

  call gl01(npol, xipol, wpol)          ! GL nodes/weights on (0,1), sum(w)=1
  ndir = nazi*npol
  allocate(omx(ndir), omy(ndir), wt(ndir), mirx(ndir), miry(ndir))
  dom = 2.0_dp*pi/real(nazi,dp)
  wazi = dom
  m = 0
  do p = 1, npol
     xi  = xipol(p)
     sth = sqrt(max(0.0_dp,1.0_dp-xi*xi))
     do a = 1, nazi
        om = (real(a,dp)-0.5_dp)*dom
        m = m + 1
        omx(m) = sth*cos(om)
        omy(m) = sth*sin(om)
        wt(m)  = 2.0_dp*wpol(p)*wazi     ! factor 2 = xi<0 mirror
     enddo
  enddo

  ! mirror indices (same direction set): x-mirror flips omx, y-mirror flips omy
  do m = 1, ndir
     tox = -omx(m); toy =  omy(m); call finddir(tox,toy,bi); mirx(m)=bi
     tox =  omx(m); toy = -omy(m); call finddir(tox,toy,bi); miry(m)=bi
  enddo
end subroutine build_quad

subroutine finddir(tx,ty,bi)
  use prec; use data
  implicit none
  real(dp), intent(in) :: tx,ty
  integer, intent(out) :: bi
  integer :: m
  real(dp) :: d, best
  best = 1.0e30_dp; bi = 1
  do m=1,ndir
     d = (omx(m)-tx)**2 + (omy(m)-ty)**2
     if (d < best) then; best = d; bi = m; endif
  enddo
end subroutine finddir

! Gauss-Legendre nodes/weights mapped to (0,1); sum(w)=1, n=1..4
subroutine gl01(n, x, w)
  use prec
  implicit none
  integer, intent(in) :: n
  real(dp), intent(out) :: x(*), w(*)
  select case(n)
  case(1)
     x(1)=0.5_dp;                 w(1)=1.0_dp
  case(2)
     x(1)=0.211324865_dp; x(2)=0.788675135_dp
     w(1)=0.5_dp;         w(2)=0.5_dp
  case(3)
     x(1)=0.112701665_dp; x(2)=0.5_dp;          x(3)=0.887298335_dp
     w(1)=0.277777778_dp; w(2)=0.444444444_dp;  w(3)=0.277777778_dp
  case default      ! n>=4 -> use 4
     x(1)=0.069431844_dp; x(2)=0.330009478_dp;  x(3)=0.669990522_dp; x(4)=0.930568156_dp
     w(1)=0.173927423_dp; w(2)=0.326072577_dp;  w(3)=0.326072577_dp; w(4)=0.173927423_dp
  end select
end subroutine gl01

!-----------------------------------------------------------------------
! Compute and write normalised pin power distribution + benchmark metrics.
!   Assembly map (core(ai,aj)):
!     A1 = UO2-inner  core(1,3)   pc=1..17,  pr=35..51
!     A2 = MOX-A      core(2,3)   pc=18..34, pr=35..51
!     A3 = MOX-B      core(1,2)   pc=1..17,  pr=18..34
!     A4 = UO2-outer  core(2,2)   pc=18..34, pr=18..34
!
!   Normalisation: average power over all fuel-producing pins = 1.0
!   (same convention as NEA/NSC/DOC(2003)16 Table A-6)
subroutine write_pin_powers()
  use prec; use data
  implicit none
  integer  :: pc, pr, li, lj, gi, gj, g, ai, aj, aidx
  integer  :: nfuel, lr, lc
  real(dp) :: ppin(nac*npa, nac*npa)
  real(dp) :: ptot, pavg, pnorm(nac*npa, nac*npa)
  real(dp) :: apow(4), atot, pmax, pmin
  integer  :: assm_map(3,3)

  assm_map      = 0
  assm_map(1,3) = 1    ! UO2-inner  (near reflective corner)
  assm_map(2,3) = 2    ! MOX-A
  assm_map(1,2) = 3    ! MOX-B
  assm_map(2,2) = 4    ! UO2-outer

  ! raw pin power = sum_{cells in pin} sum_g nsf_g * phi_g
  ppin = 0.0_dp
  do pr = 1, nac*npa
     do pc = 1, nac*npa
        do lj = 1, cpp; do li = 1, cpp
           gi = (pc-1)*cpp + li
           gj = (pr-1)*cpp + lj
           do g = 1, ng
              ppin(pc,pr) = ppin(pc,pr) + nsf(g,matid(gi,gj))*phi(gi,gj,g)
           enddo
        enddo; enddo
     enddo
  enddo

  ! normalise: average over fuel-producing pins = 1.0
  ptot = 0.0_dp;  nfuel = 0
  do pr = 1, nac*npa; do pc = 1, nac*npa
     if (ppin(pc,pr) > 0.0_dp) then
        ptot = ptot + ppin(pc,pr);  nfuel = nfuel + 1
     endif
  enddo; enddo
  pavg  = ptot / max(1.0_dp, real(nfuel,dp))
  pnorm = ppin / max(pavg, 1.0e-30_dp)

  ! assembly power fractions
  apow = 0.0_dp
  do pr = 1, nac*npa
     aj = (pr-1)/npa + 1
     do pc = 1, nac*npa
        ai   = (pc-1)/npa + 1
        aidx = assm_map(ai,aj)
        if (aidx > 0) apow(aidx) = apow(aidx) + ppin(pc,pr)
     enddo
  enddo
  atot = sum(apow)
  if (atot > 0.0_dp) apow = apow / atot * 100.0_dp

  pmax = 0.0_dp;  pmin = 1.0e30_dp
  do pr = 1, nac*npa; do pc = 1, nac*npa
     if (pnorm(pc,pr) > 0.0_dp) then
        if (pnorm(pc,pr) > pmax) pmax = pnorm(pc,pr)
        if (pnorm(pc,pr) < pmin) pmin = pnorm(pc,pr)
     endif
  enddo; enddo

  write(*,'(a)') '======================================================'
  write(*,'(a,i5,a,f7.4,a,f7.4)') &
     '  Pin power:  fuel pins=', nfuel,'  max=',pmax,'  min=',pmin
  write(*,'(a)') '  Assembly power fractions (% of total fuel power):'
  write(*,'(a,f7.3,a)') '    A1 UO2-inner  : ', apow(1), ' %'
  write(*,'(a,f7.3,a)') '    A2 MOX-A      : ', apow(2), ' %'
  write(*,'(a,f7.3,a)') '    A3 MOX-B      : ', apow(3), ' %'
  write(*,'(a,f7.3,a)') '    A4 UO2-outer  : ', apow(4), ' %'
  write(*,'(a)') '======================================================'

  open(21, file='c5g7_pinpower.csv', status='replace')
  write(21,'(a)') 'pc,pr,assembly,local_row,local_col,norm_power'
  do pr = 1, nac*npa
     aj = (pr-1)/npa + 1
     lr = mod(pr-1, npa) + 1
     do pc = 1, nac*npa
        ai   = (pc-1)/npa + 1
        lc   = mod(pc-1, npa) + 1
        aidx = assm_map(ai,aj)
        if (aidx > 0) &
           write(21,'(i3,a,i3,a,i1,a,i2,a,i2,a,f10.6)') &
              pc,',',pr,',',aidx,',',lr,',',lc,',',pnorm(pc,pr)
     enddo
  enddo
  close(21)
  write(*,'(a)') '  wrote c5g7_pinpower.csv'
end subroutine write_pin_powers

!-----------------------------------------------------------------------
! Dump downsampled fast (g1) and thermal (g7) scalar flux for plotting.
subroutine write_flux_csv()
  use prec; use data
  implicit none
  integer :: i,j, ds, ii,jj
  real(dp) :: f1,f7
  ds = max(1, cpp)                       ! ~one value per pin cell
  open(20,file='c5g7_flux.csv',status='replace')
  write(20,'(a)') 'x_cm,y_cm,flux_g1,flux_g7'
  do j=1,ny,ds
     do i=1,nx,ds
        f1=0; f7=0
        do jj=j,min(j+ds-1,ny); do ii=i,min(i+ds-1,nx)
           f1=f1+phi(ii,jj,1); f7=f7+phi(ii,jj,7)
        enddo; enddo
        write(20,'(f8.3,a,f8.3,a,es13.6,a,es13.6)') &
           (real(i,dp)-0.5_dp)*h,',',(real(j,dp)-0.5_dp)*h,',',f1,',',f7
     enddo
  enddo
  close(20)
  write(*,'(a)') '  wrote c5g7_flux.csv  (fast & thermal flux map)'
end subroutine write_flux_csv
