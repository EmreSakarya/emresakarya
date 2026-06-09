program iaea3d_pwr
    !==========================================================================
    ! IAEA 3D PWR Benchmark (ANL-7416 Suppl. 2, Benchmark Problem 11)
    ! Two-group neutron diffusion, quarter core, mesh-centered finite
    ! differences. Reference: keff = 1.02903.
    !
    ! Geometry (verified against the original specification and the FeeNoX
    ! iaea-3dpwr-quarter CAD geometry):
    !   - Quarter core, x,y in [0,170] cm, z in [0,380] cm.
    !   - Assembly pitch 20 cm; first row/column is the HALF central
    !     assembly (10 cm). Grid lines at 0,10,30,50,...,170.
    !   - Stepped external boundary: the domain is NOT the full 170x170
    !     square; outside the stepped outline there is vacuum, not water.
    !   - Axial: bottom reflector 0-20, active core 20-360, top refl 360-380.
    !   - 4 fully inserted rods at assemblies (1,1),(5,1),(1,5),(5,5)
    !     (material 3 over the whole core height, material 5 above them).
    !   - 1 PARTIALLY inserted rod at assembly (3,3): material 3 only for
    !     z in [280,360] (80 cm from the top of the core), material 2 below,
    !     material 5 above. Full core total: 13 rodded assemblies.
    !
    ! Numerics:
    !   - Harmonic-mean diffusion coefficients on cell faces (correct
    !     treatment of fuel/reflector interfaces).
    !   - Mirror (zero current) on the x=0 and y=0 symmetry planes.
    !   - Vacuum Robin BC  -D dphi/dn = 0.4692*phi  on every other external
    !     face, built into the cell balance (no flux overwriting).
    !   - Power iteration + red-black SOR inner sweeps (OpenMP-safe,
    !     deterministic).
    !
    ! Build:  ifx -O3 -qopenmp iaea3d_pwr.f90 -o iaea3d
    !    or:  gfortran -O3 -fopenmp iaea3d_pwr.f90 -o iaea3d
    !
    ! Run:    ./iaea3d            benchmark config (4 full rods + 80 cm partial)
    !         ./iaea3d aro        all rods out
    !         ./iaea3d ari        all rods fully inserted
    !==========================================================================
    implicit none
    integer, parameter :: dp = selected_real_kind(15)

    !---------------- user knobs --------------------------------------------
    integer, parameter :: nref = 5          ! cells per 10 cm  (5 -> dx=2 cm)
    real(dp), parameter :: omega   = 1.9_dp ! SOR relaxation factor
    integer, parameter :: nsweep   = 10     ! inner sweeps per group per outer
    integer, parameter :: max_outer = 5000
    real(dp), parameter :: tol_k   = 1.0e-7_dp
    real(dp), parameter :: tol_src = 1.0e-6_dp
    !-------------------------------------------------------------------------

    real(dp), parameter :: h  = 10.0_dp/real(nref,dp)      ! dx=dy=dz [cm]
    integer, parameter :: NXc = 17*nref                    ! 170 cm
    integer, parameter :: NYc = 17*nref
    integer, parameter :: NZc = 38*nref                    ! 380 cm
    real(dp), parameter :: gam = 0.4692_dp                 ! Robin constant

    ! two-group constants, indexed by material 1..5
    ! 1 fuel1 (outer fuel), 2 fuel2 (inner fuel), 3 fuel2+rod,
    ! 4 reflector, 5 reflector+rod
    real(dp), parameter :: xD1(5)  = [1.5_dp,1.5_dp,1.5_dp,2.0_dp,2.0_dp]
    real(dp), parameter :: xD2(5)  = [0.4_dp,0.4_dp,0.4_dp,0.3_dp,0.3_dp]
    real(dp), parameter :: xS12(5) = [0.02_dp,0.02_dp,0.02_dp,0.04_dp,0.04_dp]
    real(dp), parameter :: xSa1(5) = [0.01_dp,0.01_dp,0.01_dp,0.0_dp,0.0_dp]
    real(dp), parameter :: xSa2(5) = [0.080_dp,0.085_dp,0.130_dp,0.010_dp,0.055_dp]
    real(dp), parameter :: xNSf2(5)= [0.135_dp,0.135_dp,0.135_dp,0.0_dp,0.0_dp]

    ! quarter-core assembly map, row = y band, col = x band, 0 = outside
    ! (stepped external boundary). Rod positions: (1,1),(5,1),(1,5),(5,5)
    ! fully inserted, (3,3) partially inserted.
    integer :: map9(9,9)

    integer,  allocatable :: mat(:,:,:)
    real(dp), allocatable :: phi1(:,:,:), phi2(:,:,:)
    real(dp), allocatable :: fx1(:,:,:), fy1(:,:,:), fz1(:,:,:), dg1(:,:,:)
    real(dp), allocatable :: fx2(:,:,:), fy2(:,:,:), fz2(:,:,:), dg2(:,:,:)
    real(dp), allocatable :: q(:,:,:), fis(:,:,:), fisold(:,:,:)

    real(dp) :: keff, keff_old, fsum, fsum_new, err_src, scale, fmax
    real(dp) :: zc, pw(9,9), pnorm, z_rod
    integer  :: npw(9,9)
    integer  :: i, j, k, it, m, row, col, nfuel
    logical  :: rodpos(9,9)
    character(len=32) :: rodcase

    ! assembly band boundaries along x and y [cm]
    real(dp), parameter :: bnd(0:9) = &
        [0.0_dp,10.0_dp,30.0_dp,50.0_dp,70.0_dp,90.0_dp,110.0_dp,130.0_dp,150.0_dp,170.0_dp]

    map9(1,:) = [3,2,2,2,3,2,2,1,4]
    map9(2,:) = [2,2,2,2,2,2,2,1,4]
    map9(3,:) = [2,2,2,2,2,2,1,1,4]
    map9(4,:) = [2,2,2,2,2,2,1,4,4]
    map9(5,:) = [3,2,2,2,3,1,1,4,0]
    map9(6,:) = [2,2,2,2,1,1,4,4,0]
    map9(7,:) = [2,2,1,1,1,4,4,0,0]
    map9(8,:) = [1,1,1,4,4,4,0,0,0]
    map9(9,:) = [4,4,4,4,0,0,0,0,0]

    rodpos = .false.
    rodpos(1,1) = .true.; rodpos(1,5) = .true.
    rodpos(5,1) = .true.; rodpos(5,5) = .true.
    rodpos(3,3) = .true.            ! partially inserted rod

    ! rod configuration from command line:
    !   'benchmark' (default) = 4 full rods + partial rod over z = 280-360
    !   'aro' = all rods out, 'ari' = all rods fully inserted
    call get_command_argument(1, rodcase)
    if (len_trim(rodcase) == 0) rodcase = 'benchmark'
    select case (trim(rodcase))
    case ('benchmark')
        z_rod = 280.0_dp              ! (3,3) rodded only in the top 80 cm
    case ('ari')
        z_rod = 20.0_dp               ! (3,3) rodded over the whole core height
    case ('aro')
        z_rod = 1.0e9_dp              ! (3,3) never rodded
        where (map9 == 3) map9 = 2    ! remove the 4 fully inserted rods too
        rodpos = .false.              ! plain reflector everywhere on top
    case default
        print '(3a)', ' unknown rod case "', trim(rodcase), &
                      '" (use: benchmark | aro | ari)'
        stop 1
    end select

    print '(a)', '====================================================='
    print '(a)', ' IAEA 3D PWR Benchmark - two-group finite differences'
    print '(a,f6.3,a,i4,a,i4,a,i4)', ' h =', h, ' cm,  grid ', NXc, ' x', NYc, ' x', NZc
    print '(2a)', ' rod configuration: ', trim(rodcase)
    print '(a)', '====================================================='

    !--------------------------------------------------------------------
    ! material map (with one layer of halo cells, mat=0 outside)
    !--------------------------------------------------------------------
    allocate(mat(0:NXc+1,0:NYc+1,0:NZc+1)); mat = 0

    do k = 1, NZc
        zc = (real(k,dp)-0.5_dp)*h
        do j = 1, NYc
            row = band_of((real(j,dp)-0.5_dp)*h)
            do i = 1, NXc
                col = band_of((real(i,dp)-0.5_dp)*h)
                m = map9(row,col)
                if (m == 0) cycle                       ! outside stepped boundary
                if (zc < 20.0_dp) then
                    m = 4                               ! bottom reflector
                else if (zc > 360.0_dp) then
                    if (rodpos(row,col)) then; m = 5    ! top reflector (+rod)
                    else;                      m = 4
                    end if
                else
                    ! rod assembly (3,3): rodded above z_rod (case dependent)
                    if (row == 3 .and. col == 3 .and. zc > z_rod) m = 3
                end if
                mat(i,j,k) = m
            end do
        end do
    end do

    !--------------------------------------------------------------------
    ! face coupling coefficients (harmonic mean D) and diagonals
    !--------------------------------------------------------------------
    allocate(fx1(0:NXc,NYc,NZc), fy1(NXc,0:NYc,NZc), fz1(NXc,NYc,0:NZc))
    allocate(fx2(0:NXc,NYc,NZc), fy2(NXc,0:NYc,NZc), fz2(NXc,NYc,0:NZc))
    allocate(dg1(NXc,NYc,NZc), dg2(NXc,NYc,NZc))
    fx1 = 0; fy1 = 0; fz1 = 0; fx2 = 0; fy2 = 0; fz2 = 0

    do k = 1, NZc; do j = 1, NYc; do i = 0, NXc
        if (mat(i,j,k) > 0 .and. mat(i+1,j,k) > 0) then
            fx1(i,j,k) = harm(xD1(mat(i,j,k)), xD1(mat(i+1,j,k)))/(h*h)
            fx2(i,j,k) = harm(xD2(mat(i,j,k)), xD2(mat(i+1,j,k)))/(h*h)
        end if
    end do; end do; end do
    do k = 1, NZc; do j = 0, NYc; do i = 1, NXc
        if (mat(i,j,k) > 0 .and. mat(i,j+1,k) > 0) then
            fy1(i,j,k) = harm(xD1(mat(i,j,k)), xD1(mat(i,j+1,k)))/(h*h)
            fy2(i,j,k) = harm(xD2(mat(i,j,k)), xD2(mat(i,j+1,k)))/(h*h)
        end if
    end do; end do; end do
    do k = 0, NZc; do j = 1, NYc; do i = 1, NXc
        if (mat(i,j,k) > 0 .and. mat(i,j,k+1) > 0) then
            fz1(i,j,k) = harm(xD1(mat(i,j,k)), xD1(mat(i,j,k+1)))/(h*h)
            fz2(i,j,k) = harm(xD2(mat(i,j,k)), xD2(mat(i,j,k+1)))/(h*h)
        end if
    end do; end do; end do

    dg1 = 1.0_dp; dg2 = 1.0_dp
    do k = 1, NZc; do j = 1, NYc; do i = 1, NXc
        m = mat(i,j,k)
        if (m == 0) cycle
        dg1(i,j,k) = xSa1(m) + xS12(m) &
            + fx1(i-1,j,k) + fx1(i,j,k) + fy1(i,j-1,k) + fy1(i,j,k) &
            + fz1(i,j,k-1) + fz1(i,j,k)
        dg2(i,j,k) = xSa2(m) &
            + fx2(i-1,j,k) + fx2(i,j,k) + fy2(i,j-1,k) + fy2(i,j,k) &
            + fz2(i,j,k-1) + fz2(i,j,k)
        ! vacuum (Robin) leakage on faces against the outside.
        ! x=0 / y=0 faces are mirror planes -> no term there.
        if (mat(i+1,j,k) == 0) then
            dg1(i,j,k) = dg1(i,j,k) + vleak(xD1(m))
            dg2(i,j,k) = dg2(i,j,k) + vleak(xD2(m))
        end if
        if (mat(i-1,j,k) == 0 .and. i > 1) then
            dg1(i,j,k) = dg1(i,j,k) + vleak(xD1(m))
            dg2(i,j,k) = dg2(i,j,k) + vleak(xD2(m))
        end if
        if (mat(i,j+1,k) == 0) then
            dg1(i,j,k) = dg1(i,j,k) + vleak(xD1(m))
            dg2(i,j,k) = dg2(i,j,k) + vleak(xD2(m))
        end if
        if (mat(i,j-1,k) == 0 .and. j > 1) then
            dg1(i,j,k) = dg1(i,j,k) + vleak(xD1(m))
            dg2(i,j,k) = dg2(i,j,k) + vleak(xD2(m))
        end if
        if (mat(i,j,k+1) == 0) then
            dg1(i,j,k) = dg1(i,j,k) + vleak(xD1(m))
            dg2(i,j,k) = dg2(i,j,k) + vleak(xD2(m))
        end if
        if (mat(i,j,k-1) == 0) then
            dg1(i,j,k) = dg1(i,j,k) + vleak(xD1(m))
            dg2(i,j,k) = dg2(i,j,k) + vleak(xD2(m))
        end if
    end do; end do; end do

    !--------------------------------------------------------------------
    ! power iteration
    !--------------------------------------------------------------------
    allocate(phi1(0:NXc+1,0:NYc+1,0:NZc+1), phi2(0:NXc+1,0:NYc+1,0:NZc+1))
    allocate(q(NXc,NYc,NZc), fis(NXc,NYc,NZc), fisold(NXc,NYc,NZc))
    phi1 = 0; phi2 = 0
    where (mat(1:NXc,1:NYc,1:NZc) > 0)
        phi1(1:NXc,1:NYc,1:NZc) = 1.0_dp
        phi2(1:NXc,1:NYc,1:NZc) = 1.0_dp
    end where

    nfuel = 0
    fisold = 0
    do k = 1, NZc; do j = 1, NYc; do i = 1, NXc
        m = mat(i,j,k)
        if (m > 0) then
            fisold(i,j,k) = xNSf2(m)*phi2(i,j,k)
            if (xNSf2(m) > 0.0_dp) nfuel = nfuel + 1
        end if
    end do; end do; end do
    keff = 1.0_dp

    do it = 1, max_outer
        fsum = sum(fisold)

        ! ---- group 1:  removal*phi1 - div(D1 grad phi1) = fis/keff -------
        !$omp parallel do collapse(3)
        do k = 1, NZc; do j = 1, NYc; do i = 1, NXc
            q(i,j,k) = fisold(i,j,k)/keff
        end do; end do; end do
        !$omp end parallel do
        call sor_group(phi1, fx1, fy1, fz1, dg1, q)

        ! ---- group 2:  Sa2*phi2 - div(D2 grad phi2) = S12*phi1 -----------
        !$omp parallel do collapse(3) private(m)
        do k = 1, NZc; do j = 1, NYc; do i = 1, NXc
            m = mat(i,j,k)
            if (m > 0) then; q(i,j,k) = xS12(m)*phi1(i,j,k)
            else;            q(i,j,k) = 0.0_dp
            end if
        end do; end do; end do
        !$omp end parallel do
        call sor_group(phi2, fx2, fy2, fz2, dg2, q)

        ! ---- new fission source, keff, convergence -----------------------
        fsum_new = 0.0_dp
        !$omp parallel do collapse(3) private(m) reduction(+:fsum_new)
        do k = 1, NZc; do j = 1, NYc; do i = 1, NXc
            m = mat(i,j,k)
            if (m > 0) then; fis(i,j,k) = xNSf2(m)*phi2(i,j,k)
            else;            fis(i,j,k) = 0.0_dp
            end if
            fsum_new = fsum_new + fis(i,j,k)
        end do; end do; end do
        !$omp end parallel do

        keff_old = keff
        keff = keff_old*fsum_new/fsum

        fmax = maxval(fis)
        err_src = maxval(abs(fis - fisold*(fsum_new/fsum)))/fmax

        ! normalize: mean fission source over fuel cells = 1
        scale = real(nfuel,dp)/fsum_new
        !$omp parallel do collapse(3)
        do k = 1, NZc; do j = 1, NYc; do i = 1, NXc
            phi1(i,j,k) = phi1(i,j,k)*scale
            phi2(i,j,k) = phi2(i,j,k)*scale
            fisold(i,j,k) = fis(i,j,k)*scale
        end do; end do; end do
        !$omp end parallel do

        if (mod(it,10) == 0 .or. it == 1) &
            print '(a,i5,a,f10.6,a,es10.3,a,es10.3)', &
                ' iter', it, '  keff =', keff, '  dk =', abs(keff-keff_old), &
                '  src err =', err_src

        if (abs(keff-keff_old) < tol_k .and. err_src < tol_src) exit
    end do

    print '(a)', '====================================================='
    print '(4a)', ' FINAL (', trim(rodcase), '):'
    print '(a,f10.6,a,i5,a)', '   keff =', keff, '   (', it, ' outer iterations)'
    if (trim(rodcase) == 'benchmark') &
        print '(a,f10.6)', '   reference keff = 1.029030  ->  diff [pcm] =', &
                           (keff-1.02903_dp)*1.0e5_dp
    print '(a)', '====================================================='

    ! assembly-averaged radial power map (axially integrated, normalized)
    pw = 0; npw = 0
    do k = 1, NZc; do j = 1, NYc; do i = 1, NXc
        m = mat(i,j,k)
        if (m > 0) then
            if (xNSf2(m) > 0.0_dp) then
                row = band_of((real(j,dp)-0.5_dp)*h)
                col = band_of((real(i,dp)-0.5_dp)*h)
                pw(row,col) = pw(row,col) + xNSf2(m)*phi2(i,j,k)
                npw(row,col) = npw(row,col) + 1
            end if
        end if
    end do; end do; end do
    pnorm = sum(pw)/real(sum(npw),dp)
    print '(a)', ' normalized assembly powers (row 1 = core center line):'
    do row = 1, 9
        print '(9f8.4)', (merge(pw(row,col)/(real(max(npw(row,col),1),dp)*pnorm), &
                                0.0_dp, npw(row,col) > 0), col = 1, 9)
    end do

    open(unit=20, file='iaea3d_midplane.dat', access='stream', status='replace')
    write(20) phi1(1:NXc,1:NYc,NZc/2)
    write(20) phi2(1:NXc,1:NYc,NZc/2)
    close(20)

contains

    pure integer function band_of(x) result(b)
        real(dp), intent(in) :: x
        integer :: n
        b = 9
        do n = 1, 9
            if (x < bnd(n)) then; b = n; return; end if
        end do
    end function band_of

    pure real(dp) function harm(da, db)
        real(dp), intent(in) :: da, db
        harm = 2.0_dp*da*db/(da + db)
    end function harm

    ! per-volume leakage coefficient of a vacuum (Robin) boundary face:
    ! J = gam*phi_s with phi_s extrapolated from the cell center
    pure real(dp) function vleak(d)
        real(dp), intent(in) :: d
        vleak = (gam*d/(d + 0.5_dp*gam*h))/h
    end function vleak

    subroutine sor_group(phi, fx, fy, fz, dg, src)
        real(dp), intent(inout) :: phi(0:,0:,0:)
        real(dp), intent(in)    :: fx(0:,:,:), fy(:,0:,:), fz(:,:,0:)
        real(dp), intent(in)    :: dg(:,:,:), src(:,:,:)
        real(dp) :: pnew
        integer  :: s, color, i0, ii, jj, kk

        do s = 1, nsweep
            do color = 0, 1
                !$omp parallel do collapse(2) private(ii,jj,kk,i0,pnew)
                do kk = 1, NZc
                    do jj = 1, NYc
                        i0 = 1 + mod(jj + kk + color, 2)
                        do ii = i0, NXc, 2
                            if (mat(ii,jj,kk) > 0) then
                                pnew = (src(ii,jj,kk) &
                                    + fx(ii-1,jj,kk)*phi(ii-1,jj,kk) &
                                    + fx(ii,  jj,kk)*phi(ii+1,jj,kk) &
                                    + fy(ii,jj-1,kk)*phi(ii,jj-1,kk) &
                                    + fy(ii,jj,  kk)*phi(ii,jj+1,kk) &
                                    + fz(ii,jj,kk-1)*phi(ii,jj,kk-1) &
                                    + fz(ii,jj,kk  )*phi(ii,jj,kk+1)) / dg(ii,jj,kk)
                                phi(ii,jj,kk) = phi(ii,jj,kk) + omega*(pnew - phi(ii,jj,kk))
                            end if
                        end do
                    end do
                end do
                !$omp end parallel do
            end do
            ! mirror ghost cells (zero current across x=0 and y=0): not
            ! needed explicitly because the face couplings fx(0,..), fy(..,0)
            ! are zero there; ghost values are never used with weight > 0.
        end do
    end subroutine sor_group

end program iaea3d_pwr
