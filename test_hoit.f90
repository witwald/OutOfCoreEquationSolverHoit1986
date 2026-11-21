! Test program for the Hoit (1986) blocked assembly routines and the
! modified Wilson-Dovey solver SORED.
!
! Test A: the paper's own example - the two-bay, two-storey plane frame of
!         Figures 1-3 (18 equations, 10 elements, block size 25). Checks
!         every intermediate array against the values printed in the paper:
!         profile pointers (Fig. 1 profile), KA table and diagonal
!         locations (Fig. 1), element-to-block pointers NELP/NBKP
!         (Figs. 2 and 3). Each assembled block is compared term by term
!         with a dense assembly, and the system is solved out-of-core
!         (5 blocks, as in the figure) and in-core.
!
! Test B: a larger frame (3 bays x 15 storeys, 180 equations, 7 load cases)
!         plus spring elements: a long-range spring (wide profile, varied
!         "lowest connected block") and a spring whose two DOFs are the same
!         equation (exercises ASEMBL's "ELEM HAS 2 DOF SAME" line). Element
!         loads come through ASEMLD, nodal loads through MLOAD (with load
!         generation, a skipped line, and jumps between load blocks).
!         Solved for several memory budgets: in-core, out-of-core with a few
!         or many blocks, and out-of-core forced by the record-size limit.
!         Displacements and the assembled mass vector are checked against a
!         dense double-precision reference.
!
! Test C: the same model solved with separate calls KKK=3, 4, 5.
!
program test_hoit
  implicit none
  ! ---- SIMPAL common blocks (must match the Fortran 77 routines) ----------
  integer :: ISTAG,NUMNP,NEQ,NEQS,LOAD,MAXD,MAXNA,NUMEL,NEL(20),LENEL(20), &
             NKBLK,NLBLK,KSIZ,KMEM,KYDYN,KYSUB,KYLIN
  common /CLT/ ISTAG,NUMNP,NEQ,NEQS,LOAD,MAXD,MAXNA,NUMEL,NEL,LENEL, &
               NKBLK,NLBLK,KSIZ,KMEM,KYDYN,KYSUB,KYLIN
  character(len=1) :: NAM(12), ELNAM(4,20)
  common /PARC/ NAM, ELNAM
  integer :: NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
  common /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
  integer :: MAXCDA
  common /DAPAR/ MAXCDA
  integer :: NDAR,NDAW,NCRLD
  common /DASTAT/ NDAR,NDAW,NCRLD
  integer, external :: LENELR

  integer, parameter :: MXE = 400, MXQ = 400, MXL = 8, MXN = 200
  ! element data kept in memory for the dense reference
  integer :: nel_t, elND(MXE), elLM(6,MXE), elGRP(MXE)
  double precision :: elSK(6,6,MXE), elR(6,MXL,MXE), elEM(6,MXE)
  integer :: ID(MXN,6)
  double precision :: xn(MXN), yn(MXN)
  ! reference
  double precision :: Kd(MXQ,MXQ), Fd(MXQ,MXL), Md(MXQ), Xd(MXQ,MXL)
  double precision :: X(MXQ,MXL), TMAS(MXQ)
  integer :: IST(12)
  ! load lines for MLOAD: N1 N2 N3 LOD + 6 values
  integer :: nll, lli(4,20)
  double precision :: llr(6,20)
  integer :: nfail = 0
  integer(8) :: iseed = 12345

  NTM = 30; NTR = 0; NIN = 32; NOT = 31; NT1 = 33; NFL = 0; NT2 = 0; NT3 = 34
  open (NTM, file='hoit.log', status='replace')
  open (NOT, file='hoit.out', status='replace')
  ISTAG = 0; NEQS = 0; KYDYN = 0; KYSUB = 0; KYLIN = 0

  call test_a()
  call test_b()

  print *
  if (nfail == 0) then
    print *, 'ALL TESTS PASSED'
  else
    print '(a,i0,a)', ' ', nfail, ' CHECK(S) FAILED'
    stop 1
  end if

contains

  ! ======================================================================
  subroutine test_a()
    integer :: NA(MXQ), LM(6), KA(4,50), KSIZB(50), NELP(MXE), NBKP(2000)
    integer :: MA(MXQ), KAX(4), NALAP, KSZKM, ib, j, i, k, nw, ia, ih
    double precision :: STIF(200), ESTF(21), R(6*MXL), EMAS(6), err
    integer, parameter :: na_exp(18) = (/ 1,3,6,10,15,21,25, 30,36,46, 57,69, &
                                          73,78,84,94, 105,117 /)
    integer, parameter :: ka_exp(4,5) = reshape((/ 1,1,1,7, 2,1,8,10, &
                      3,1,11,12, 4,1,13,16, 5,1,17,18 /), (/4,5/))
    integer, parameter :: ksb_exp(5) = (/ 25,21,23,25,23 /)
    integer, parameter :: ma_exp(18) = (/ 1,3,6,10,15,21,25, 5,11,21, 11,23, &
                                          4,9,15,25, 11,23 /)
    integer, parameter :: nelp_exp(10) = (/ 2,4,6,9,10,11,15,18,20,23 /)
    integer, parameter :: nbkp_exp(23) = (/ 4,5, 4,5, 1,2, 1,2,3, 1, 1, &
                                            1,2,4,5, 2,3,4, 1,2, 1,2,3 /)
    logical :: ok

    print '(/a)', ' Test A: Hoit (1986) Figures 1-3 example frame'
    call build_fig1_frame()
    call write_element_files()

    ! --- profile (Figure 1) ---
    call PROFIL(NA, LM, MAXD, NEQ)
    ok = all(NA(1:18) == na_exp) .and. KSIZ == 117
    call report('PROFIL diagonal pointers = Fig. 1 profile (117 terms)', ok)

    ! --- blocking with a 25-term block (Figure 1 KA array) ---
    call BLKSET(KA, KSIZB, NA, KSZKM, 25, NKBLK, NEQ, MAXNA)
    ok = NKBLK == 5 .and. all(KA(:,1:5) == ka_exp) .and. &
         all(KSIZB(1:5) == ksb_exp) .and. KSZKM == 25 .and. MAXNA == 7
    call report('BLKSET KA table, block sizes = Fig. 1', ok)

    ! --- element-to-block pointer (Figures 2 and 3) ---
    call BLKPNT(LM, KA, NELP, NBKP, NALAP, 0, nel_t)
    ok = NALAP == 23 .and. all(NELP(1:10) == nelp_exp) .and. &
         all(NBKP(1:23) == nbkp_exp)
    call report('BLKPNT NELP/NBKP = Figs. 2 and 3', ok)

    ! --- blocked assembly, then compare each block with dense K ---
    KSIZ = 25
    MAXCDA = MAXNA
    call SOPEN(NT3, (/'K',' ',' ',' '/), -(KSIZ*8 + (MAXCDA+5)*4))
    call ASEMBL(STIF, NA, ESTF, R, EMAS, KSIZB, LM, KA, NELP, NBKP, MA, &
                MAXD, NEQ, i, LOAD, KSIZ)
    call dense_reference()
    ok = .true.; err = 0.d0
    do ib = 1, NKBLK
      call DABLK(NT3, ib, STIF, nw, 200, MA, MAXCDA, KAX, 1)
      if (any(KAX /= KA(:,ib)) .or. nw /= KSIZB(ib)) ok = .false.
      do j = KAX(3), KAX(4)
        k = j - KAX(3) + 1
        if (MA(k) /= ma_exp(j)) ok = .false.
        ia = MA(k)
        ih = MA(k) - merge(0, MA(max(k-1,1)), k == 1)       ! column height
        do i = j, j-ih+1, -1
          err = max(err, abs(STIF(ia) - Kd(i,j)))
          ia = ia - 1
        end do
        if (j-ih >= 1) err = max(err, maxval(abs(Kd(1:j-ih, j))))  ! outside profile
      end do
    end do
    call FCLOSE(NT3)
    call report('ASEMBL block diagonal locations M(j) = Fig. 1', ok)
    call report('ASEMBL blocks equal dense K term by term', err < 1.d-9*maxval(abs(Kd(1:NEQ,1:NEQ))))

    ! --- full assemble + solve, out-of-core 5 blocks (KMEM=50 -> KSIZ=25) ---
    call run_case('A1 out-of-core, KMEM=50 (Fig. 1 blocking)', 50, 8000, 0, 1, 5)
    call run_case('A2 in-core, KMEM=10000', 10000, 8000, 0, 1, 1)
  end subroutine test_a

  ! ======================================================================
  subroutine test_b()
    integer :: NA(MXQ), LM(6), S
    print '(/a)', ' Test B: 3-bay x 15-storey frame + springs, 7 load cases'
    call build_grid_frame(3, 15)
    call write_element_files()
    call set_load_lines()
    call PROFIL(NA, LM, MAXD, NEQ)
    S = KSIZ
    print '(a,i0,a,i0)', '   NEQ = ', NEQ, ', profile terms = ', S
    call dense_reference()
    !                label                              KMEM        MAXREC NODLD ISOLV  expect-blocks(0=any)
    call run_case('B1 in-core, all loads fit       ', S+NEQ*LOAD+5, 8000, 1, 1, 1)
    call run_case('B2 in-core, one load vector fits', S+NEQ,        8000, 1, 1, 1)
    call run_case('B3 out-of-core, KMEM = profile  ', S,            8000, 1, 1, 0)
    call run_case('B4 out-of-core, KMEM = 2*NEQ    ', 2*NEQ,        8000, 1, 1, 0)
    call run_case('B5 out-of-core by MAXREC limit  ', 100000,       S/3,  1, 1, 0)
    print '(/a)', ' Test C: same model, solver called as KKK=3, then 4, then 5'
    call run_case('C1 out-of-core, KMEM = profile  ', S,            8000, 1, 2, 0)
    call run_case('C2 in-core                      ', S+NEQ*LOAD+5, 8000, 1, 2, 1)
  end subroutine test_b

  ! ======================================================================
  subroutine run_case(label, kmemin, maxrec, nodld, isolv, nblk_expect)
    character(len=*), intent(in) :: label
    integer, intent(in) :: kmemin, maxrec, nodld, isolv, nblk_expect
    double precision :: ex, em, Xv(MXQ*MXL)
    integer :: IDc(NUMNP,6)          ! ID(NUMNP,6) as declared in MLOAD
    logical :: ok
    if (nodld == 1) call write_load_file()
    Xv = 0.d0
    IDc = ID(1:NUMNP,:)
    call ADDDRV(kmemin, maxrec, IDc, nodld, isolv, Xv, TMAS, IST)
    X(1:NEQ,1:LOAD) = reshape(Xv(1:NEQ*LOAD), (/ NEQ, LOAD /))   ! X(NEQ,LOAD) in ADDDRV
    ex = maxval(abs(X(1:NEQ,1:LOAD) - Xd(1:NEQ,1:LOAD))) / maxval(abs(Xd(1:NEQ,1:LOAD)))
    em = maxval(abs(TMAS(1:NEQ) - Md(1:NEQ))) / maxval(abs(Md(1:NEQ)))
    ok = ex < 1.d-10 .and. em < 1.d-12
    if (nblk_expect > 0) ok = ok .and. IST(1) == nblk_expect
    ! in-core, single KKK=1 call: the coefficient matrix is never re-read.
    ! (separate KKK=4 and KKK=5 calls each read the single block once)
    if (IST(3) == 1 .and. isolv == 1) ok = ok .and. IST(6) == 0
    if (IST(3) == 1 .and. isolv == 2) ok = ok .and. IST(6) == 2
    print '(3x,a,a,i2,a,i3,a,i2,a,i4,a,es8.1,a,es8.1,2x,a)', label, &
      '  in-core=', IST(3), '  K blocks=', IST(1), '  load blocks=', IST(2), &
      '  K re-reads=', IST(6), '  err x=', ex, '  mass=', em, merge('ok  ','FAIL',ok)
    if (.not. ok) nfail = nfail + 1
  end subroutine run_case

  subroutine report(label, ok)
    character(len=*), intent(in) :: label
    logical, intent(in) :: ok
    print '(3x,a,1x,a)', merge('ok  ','FAIL',ok), label
    if (.not. ok) nfail = nfail + 1
  end subroutine report

  ! ======================================================================
  ! Models
  ! ======================================================================
  subroutine build_fig1_frame()
    integer :: n
    ! nodes: 1-3 supports, 4 mid-left, 5 mid-mid, 6 mid-right,
    !        7 top-left, 8 top-mid, 9 top-right. DOF numbers as in Fig. 1.
    NUMNP = 9; NEQ = 18; LOAD = 2; MAXD = 6
    xn(1:9) = (/ 0.d0, 6.d0, 12.d0, 0.d0, 6.d0, 12.d0, 0.d0, 6.d0, 12.d0 /)
    yn(1:9) = (/ 0.d0, 0.d0, 0.d0, 4.d0, 4.d0, 4.d0, 8.d0, 8.d0, 8.d0 /)
    ID = 0
    call set_id(4, 16); call set_id(5, 7);  call set_id(6, 4)
    call set_id(7, 13); call set_id(8, 10); call set_id(9, 1)
    nel_t = 0
    NUMEL = 2
    call set_name(1, 'COLS'); call set_name(2, 'BEAM')
    NEL(1) = 6; NEL(2) = 4
    ! columns (elements 1-6), then beams (7-10) - numbering of Fig. 1
    call add_frame(1,4,1); call add_frame(4,7,1); call add_frame(2,5,1)
    call add_frame(5,8,1); call add_frame(3,6,1); call add_frame(6,9,1)
    call add_frame(4,5,2); call add_frame(7,8,2); call add_frame(5,6,2)
    call add_frame(8,9,2)
    n = nel_t
  end subroutine build_fig1_frame

  subroutine build_grid_frame(nbay, nsto)
    integer, intent(in) :: nbay, nsto
    integer :: i, j, nd, eq, e1, e2
    NUMNP = (nbay+1)*(nsto+1); LOAD = 7; MAXD = 6
    ID = 0; eq = 0
    do j = 0, nsto
      do i = 0, nbay
        nd = j*(nbay+1) + i + 1
        xn(nd) = 5.d0*i + 0.3d0*j          ! slightly leaning columns
        yn(nd) = 3.5d0*j
        if (j > 0) then
          call set_id(nd, eq+1); eq = eq + 3
        end if
      end do
    end do
    NEQ = eq
    nel_t = 0; NUMEL = 3
    call set_name(1, 'COLS'); call set_name(2, 'BEAM'); call set_name(3, 'SPRG')
    do j = 1, nsto
      do i = 0, nbay
        call add_frame((j-1)*(nbay+1)+i+1, j*(nbay+1)+i+1, 1)
      end do
    end do
    NEL(1) = nel_t
    do j = 1, nsto
      do i = 1, nbay
        call add_frame(j*(nbay+1)+i, j*(nbay+1)+i+1, 2)
      end do
    end do
    NEL(2) = nel_t - NEL(1)
    ! springs: long-range, a "2 DOF same" element, and a short one
    e1 = ID(nsto*(nbay+1)+1, 1)            ! top-left u
    e2 = ID(1*(nbay+1)+nbay+1, 1)          ! first-storey right u
    call add_spring(min(e1,e2), max(e1,e2), 4.d2, -1.5d2, 3.d2)
    e1 = ID(8*(nbay+1)+2, 2)
    call add_spring(e1, e1, 7.d2, -2.d2, 5.d2)
    e1 = ID(5*(nbay+1)+3, 6); e2 = ID(6*(nbay+1)+3, 6)
    call add_spring(e1, e2, 9.d1, -3.d1, 6.d1)
    NEL(3) = nel_t - NEL(1) - NEL(2)
  end subroutine build_grid_frame

  subroutine set_id(nd, e1)
    integer, intent(in) :: nd, e1
    ID(nd,1) = e1; ID(nd,2) = e1+1; ID(nd,6) = e1+2   ! u, v, rotation z
  end subroutine set_id

  subroutine set_name(g, s)
    integer, intent(in) :: g
    character(len=4), intent(in) :: s
    integer :: k
    do k = 1, 4
      ELNAM(k,g) = s(k:k)
    end do
  end subroutine set_name

  subroutine add_frame(n1, n2, grp)
    integer, intent(in) :: n1, n2, grp
    double precision :: E, AR, ZI
    integer :: k, i
    nel_t = nel_t + 1
    if (grp == 1) then
      E = 2.0d5; AR = 1.2d-2; ZI = 2.0d-4
    else
      E = 2.0d5; AR = 8.0d-3; ZI = 3.5d-4
    end if
    call frame_k(xn(n1), yn(n1), xn(n2), yn(n2), E, AR, ZI, elSK(:,:,nel_t))
    elND(nel_t) = 6
    elLM(:,nel_t) = (/ ID(n1,1), ID(n1,2), ID(n1,6), ID(n2,1), ID(n2,2), ID(n2,6) /)
    elGRP(nel_t) = grp
    do k = 1, 6
      elEM(k,nel_t) = 1.d0 + rnd()
    end do
    do k = 1, LOAD
      do i = 1, 6
        elR(i,k,nel_t) = 10.d0*(rnd()-0.5d0)
      end do
    end do
  end subroutine add_frame

  subroutine add_spring(e1, e2, k11, k12, k22)
    integer, intent(in) :: e1, e2
    double precision, intent(in) :: k11, k12, k22
    integer :: k
    nel_t = nel_t + 1
    elND(nel_t) = 2
    elLM(:,nel_t) = 0
    elLM(1:2,nel_t) = (/ e1, e2 /)
    elSK(:,:,nel_t) = 0.d0
    elSK(1,1,nel_t) = k11; elSK(1,2,nel_t) = k12
    elSK(2,1,nel_t) = k12; elSK(2,2,nel_t) = k22
    elGRP(nel_t) = 3
    elEM(:,nel_t) = 0.d0
    elEM(1:2,nel_t) = (/ 0.5d0, 0.25d0 /)
    elR(:,:,nel_t) = 0.d0
    do k = 1, LOAD
      elR(1:2,k,nel_t) = (/ rnd(), rnd() /)
    end do
  end subroutine add_spring

  subroutine frame_k(x1, y1, x2, y2, E, AR, ZI, SK)
    double precision, intent(in) :: x1, y1, x2, y2, E, AR, ZI
    double precision, intent(out) :: SK(6,6)
    double precision :: L, c, s, KL(6,6), T(6,6)
    L = sqrt((x2-x1)**2 + (y2-y1)**2); c = (x2-x1)/L; s = (y2-y1)/L
    KL = 0.d0
    KL(1,1) = E*AR/L;  KL(4,4) = KL(1,1); KL(1,4) = -KL(1,1); KL(4,1) = KL(1,4)
    KL(2,2) = 12*E*ZI/L**3; KL(5,5) = KL(2,2); KL(2,5) = -KL(2,2); KL(5,2) = KL(2,5)
    KL(2,3) = 6*E*ZI/L**2; KL(3,2) = KL(2,3); KL(2,6) = KL(2,3); KL(6,2) = KL(2,3)
    KL(3,5) = -KL(2,3); KL(5,3) = KL(3,5); KL(5,6) = -KL(2,3); KL(6,5) = KL(5,6)
    KL(3,3) = 4*E*ZI/L; KL(6,6) = KL(3,3); KL(3,6) = 2*E*ZI/L; KL(6,3) = KL(3,6)
    T = 0.d0
    T(1,1) = c; T(1,2) = s; T(2,1) = -s; T(2,2) = c; T(3,3) = 1.d0
    T(4:6,4:6) = T(1:3,1:3)
    SK = matmul(transpose(T), matmul(KL, T))
  end subroutine frame_k

  ! ---- write each element group to its direct access file ----------------
  subroutine write_element_files()
    integer :: g, j, k, l, iel, n, first
    integer :: LM(6)
    double precision :: ESTF(21), R(6*MXL), EMAS(6)
    first = 0
    do g = 1, NUMEL
      LENEL(g) = LENELR(MAXD, LOAD)
      call SOPEN(NT1, ELNAM(1,g), -LENEL(g))
      do j = 1, NEL(g)
        n = first + j
        LM = elLM(:,n)
        iel = 0
        ESTF = 0.d0
        do k = 1, elND(n)                 ! upper triangle by columns
          do l = 1, k
            iel = iel + 1
            ESTF(iel) = elSK(l,k,n)
          end do
        end do
        R = 0.d0
        do k = 1, LOAD
          R((k-1)*MAXD+1:(k-1)*MAXD+6) = elR(:,k,n)
        end do
        EMAS = elEM(:,n)
        call ELWRIT(NT1, j, LM, ESTF, R, EMAS, elND(n), LOAD, MAXD)
      end do
      call FCLOSE(NT1)
      first = first + NEL(g)
    end do
  end subroutine write_element_files

  ! ---- nodal load lines for MLOAD --------------------------------------------
  subroutine set_load_lines()
    integer :: k
    nll = 0
    call add_ll(7, 7, 0, 99, (/ 1.d0,1.d0,0.d0,0.d0,0.d0,1.d0 /))  ! LOD>LOAD: skipped
    call add_ll(5, 60, 5, 1, (/ 12.d0,-4.d0,0.d0,0.d0,0.d0,2.5d0 /)) ! generation
    call add_ll(30, 0, 0, 5, (/ 0.d0,-20.d0,0.d0,0.d0,0.d0,0.d0 /))
    call add_ll(9, 0, 0, 0, (/ 3.d0,0.d0,0.d0,0.d0,0.d0,-1.d0 /))   ! LOD 0 -> 1
    call add_ll(1, 4, 1, 2, (/ 50.d0,50.d0,0.d0,0.d0,0.d0,0.d0 /)) ! support nodes: no DOF
    do k = 1, 7                                                      ! hop between blocks
      call add_ll(60-k, 0, 0, 8-k, (/ dble(k),-dble(k),0.d0,0.d0,0.d0,0.1d0*k /))
    end do
    call add_ll(33, 41, 4, 3, (/ -6.d0,0.d0,9.d0,9.d0,9.d0,1.5d0 /))
    call add_ll(64, 64, 1, 7, (/ 2.d0,2.d0,0.d0,0.d0,0.d0,2.d0 /))
  end subroutine set_load_lines

  subroutine add_ll(n1, n2, n3, lod, rt)
    integer, intent(in) :: n1, n2, n3, lod
    double precision, intent(in) :: rt(6)
    nll = nll + 1
    lli(:,nll) = (/ n1, n2, n3, lod /)
    llr(:,nll) = rt
  end subroutine add_ll

  subroutine write_load_file()
    integer :: k
    close (NIN)
    open (NIN, file='loads.inp', status='replace')
    do k = 1, nll
      write (NIN, '(4i6,6es16.8)') lli(:,k), llr(:,k)
    end do
    write (NIN, '(a)') '0 /'
    rewind NIN
  end subroutine write_load_file

  ! ---- dense reference: K, F (element + nodal loads), mass, solution -------
  subroutine dense_reference()
    integer :: n, a, b, ia, ib, k, j, lod, n2, n3
    Kd = 0.d0; Fd = 0.d0; Md = 0.d0
    do n = 1, nel_t
      do a = 1, elND(n)
        ia = elLM(a,n)
        if (ia <= 0 .or. ia > NEQ) cycle
        Md(ia) = Md(ia) + elEM(a,n)
        Fd(ia,1:LOAD) = Fd(ia,1:LOAD) + elR(a,1:LOAD,n)
        do b = 1, elND(n)
          ib = elLM(b,n)
          if (ib <= 0 .or. ib > NEQ) cycle
          Kd(ia,ib) = Kd(ia,ib) + elSK(a,b,n)
        end do
      end do
    end do
    if (NUMEL == 3) then                     ! test B model has nodal loads
      do k = 1, nll
        if (lli(1,k) <= 0) exit
        lod = lli(4,k); if (lod <= 0) lod = 1
        if (lod > LOAD) cycle
        n2 = lli(2,k); if (n2 <= lli(1,k)) n2 = lli(1,k)
        n3 = lli(3,k); if (n3 <= 0) n3 = 1
        do j = lli(1,k), n2, n3
          do a = 1, 6
            ia = ID(j,a)
            if (ia > 0 .and. ia <= NEQ) Fd(ia,lod) = Fd(ia,lod) + llr(a,k)
          end do
        end do
      end do
    end if
    call gauss(Kd, Fd, Xd, NEQ, LOAD)
  end subroutine dense_reference

  subroutine gauss(A, F, Xs, n, nr)               ! partial pivoting
    integer, intent(in) :: n, nr
    double precision, intent(in) :: A(MXQ,MXQ), F(MXQ,MXL)
    double precision, intent(out) :: Xs(MXQ,MXL)
    double precision :: W(n,n), Y(n,nr), f_, rowW(n), rowY(nr)
    integer :: i, k, p
    W = A(1:n,1:n); Y = F(1:n,1:nr)
    do k = 1, n-1
      p = k - 1 + maxloc(abs(W(k:n,k)), 1)
      if (p /= k) then
        rowW = W(k,:); W(k,:) = W(p,:); W(p,:) = rowW
        rowY = Y(k,:); Y(k,:) = Y(p,:); Y(p,:) = rowY
      end if
      do i = k+1, n
        f_ = W(i,k)/W(k,k)
        W(i,k:n) = W(i,k:n) - f_*W(k,k:n)
        Y(i,:) = Y(i,:) - f_*Y(k,:)
      end do
    end do
    do i = n, 1, -1
      Xs(i,1:nr) = (Y(i,:) - matmul(W(i,i+1:n), Xs(i+1:n,1:nr))) / W(i,i)
    end do
  end subroutine gauss

  double precision function rnd()
    iseed = mod(iseed*16807_8, 2147483647_8)
    rnd = dble(iseed)/2147483647.d0
  end function rnd

end program test_hoit
