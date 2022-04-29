!===============================================================================
! Description
! -----------
!   Different from Fortran program, to call subroutine in Julia we have to use
!   explicit-shape array. Therefore, many subroutine interfaces are changed
!   compared to the original fortran program.
!
!   To hide unnecessary details, we pack table-generating subroutines inside a
!   single subroutine initialize().
!
!   Version:  0.0
!   Created:  2021-04-29 17:01
!    Author:  Zhiyuan Yao, zhiyuan.yao@icloud.com
! Institute:  Insitute for Advanced Study, Tsinghua University
!===============================================================================
MODULE symmetryED
    USE, INTRINSIC:: ISO_FORTRAN_ENV
    USE, INTRINSIC:: ISO_C_BINDING
    IMPLICIT NONE
    ! -------------------------------------------------------------------------
    ! System Parameters
    ! -------------------------------------------------------------------------
    integer, parameter  :: L = 20
    double precision, parameter :: mass = 0.d0
    integer, parameter  :: momentum  =  0
    integer, parameter  :: inversion =  1
    logical, parameter  :: hasPBC = .True.     ! limitted to PBC for translational symmetry

    integer(KIND=INT64) :: stateN              ! Hilbert space dimension
    integer :: classN !number of distinct translational classes
    integer :: basisN !basis number in k=0+ (translation+parity) symmetry sector

    ! -------------------------------------------------------------------------
    ! k = 0 or π, I = even or odd
    ! -------------------------------------------------------------------------
    INTEGER, parameter :: I_EVEN = 1, I_ODD = -1
    INTEGER :: k_ZERO, k_PI
    INTEGER, parameter :: NOTHING = -1

    !---------------------------------------------------------------------------
    ! Hierarchical tables and arrays for efficient storage and fast access
    !---------------------------------------------------------------------------
    integer(kind=INT64), allocatable, dimension(:), TARGET :: stateTable
    integer, allocatable, dimension(:), TARGET :: istate_state
    integer, allocatable, dimension(:), TARGET :: iclassTable, shiftTable
    integer, allocatable, dimension(:), TARGET :: istate_class, period_class
    ! state_all_class(0, iclass) store the number of members in [a]
    integer(kind=INT64), allocatable, TARGET :: state_all_class(:, :)
    integer, allocatable, dimension(:), TARGET :: shift_basis
    ! iclass == class_parity_class(iclass) means Pa = b ∈ [a]
    integer, allocatable, dimension(:), TARGET :: iclass_parity_class
    integer, allocatable, dimension(:), TARGET :: ibasis_class, iclass_basis
    logical, allocatable, dimension(:), TARGET :: has_parity_basis
    double precision, allocatable, dimension(:), TARGET :: norm_basis

    ! -------------------------------------------------------------------------
    ! eigvals, eigvects and the translated wavefunction
    !---------------------------------------------------------------------------
    double precision, allocatable, TARGET:: Ham(:, :), eigens(:)
    double precision, allocatable, TARGET :: psi(:)

contains
    !=========================================================================!
    ! read in parameters and allocate necessary array
    !=========================================================================!
    SUBROUTINE initialize()
        implicit none

        k_ZERO = 0; k_PI  = L/2
        if (momentum /= k_ZERO .AND. momentum /= k_PI) then
            call raiseError("initialize(): momentum limitted to 0 and π")
        end if

        if (.NOT. allocated(stateTable)) then
            print *, "L = ", L
            call createStateTable(stateN, stateTable)
            allocate (iclassTable(stateN)); allocate (shiftTable(stateN));
        end if

        call createTranslationClass(stateN, stateTable)
        call createZeroPiMomentumBasis()

        call constructHam(Ham, basisN)

        allocate(eigens(basisN))
        call getDSYEigenSystem(Ham, eigens)
        print '(G20.8)', eigens(1:min(10, basisN))

        allocate(psi(stateN))

        call get_szList()
    END SUBROUTINE initialize


    !=========================================================================!
    ! filter the physical states by problem-specific rules.
    !=========================================================================!
    SUBROUTINE createStateTable(stateN, stateTable)
        integer(KIND=INT64), intent(out) :: stateN
        integer(kind=INT64), allocatable, intent(out) :: stateTable(:)
        integer(kind=INT64) :: state
        integer, parameter  :: NOTHING = -1
        integer :: k, istate
        integer :: FibonacciTable(40) = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, &
        144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, &
        46368, 75025, 121393, 196418, 317811, 514229, 832040, 1346269, 2178309, &
        3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986, 102334155]

        allocate (stateTable(FibonacciTable(L - 1) + FibonacciTable(L + 1)))
        stateTable = NOTHING

        allocate(istate_state(0:2_INT64**L-1)); istate_state = NOTHING

        istate = 0
        stateloop: do state = 0, 2_INT64**L - 1
            if (ibits(state, 0, 1)*ibits(state, L - 1, 1) == 1) cycle stateloop
            do k = 0, L - 2
                if (ibits(state, k, 1)*ibits(state, k + 1, 1) == 1) cycle stateloop
            end do
            istate = istate + 1
            stateTable(istate) = state
            istate_state(state) = istate
        end do stateloop
        stateN = istate
        print "(A, I10)", " stateN = ", stateN
    END SUBROUTINE createStateTable

    !=========================================================================!
    ! get the period of a given state
    !=========================================================================!
    INTEGER FUNCTION get_period(state) result(res)
        integer(kind=INT64), intent(in) :: state
        integer(kind=INT64) :: state_shift
        integer :: n

        state_shift = state
        res = 0
        do n = 1, L
            state_shift = ishftc(state_shift, -1, SIZE=L)
            res = res + 1
            if (state_shift == state) return
        end do
    END FUNCTION get_period

    !=========================================================================!
    ! create the translational class with given stateTable(1:stateN)
    !
    ! dependent variables:
    ! ---------------------
    ! shiftTable, iclassTable, classN, istate_class, period_class
    !=========================================================================!
    SUBROUTINE createTranslationClass(stateN, stateTable)
        integer(KIND=INT64), intent(in) :: stateN
        integer(kind=INT64), dimension(1:stateN), intent(in) :: stateTable
        integer :: n, iclass, istate, istate_shift, period
        integer(kind=INT64) :: state, state_shift
        logical, dimension(stateN) :: state_marked
        integer, dimension(stateN) :: period_class_tmp
        integer, dimension(stateN) :: istate_class_tmp
        integer, parameter  :: NOTHING  = -1

        iclassTable = NOTHING;  shiftTable = NOTHING
        state_marked = .False.
        istate_class_tmp = NOTHING
        period_class_tmp = NOTHING
        !----------------------------------------------------------------------
        ! deal with translation symmetry only
        !   * iclass: index of the translation class, [a]
        !----------------------------------------------------------------------
        iclass = 0
        do istate = 1, int(stateN, KIND=4)
            if (state_marked(istate)) cycle
            state = stateTable(istate)
            period = get_period(state)
            if (mod(momentum, L/period) /= 0) cycle
            iclass = iclass + 1
            istate_class_tmp(iclass) = istate
            iclassTable(istate) = iclass
            ! the starting state a* of a translation class [a]
            ! now mark all state related by translation
            state_shift = state
            shiftTable(istate) = 0
            Tloop: do n = 1, L
                state_shift = ishftc(state_shift, -1, SIZE=L)
                if (state_shift == state) then
                    ! loop invariant: has shifted n times
                    period_class_tmp(iclass) = n
                    if (mod(L, n) /= 0 .OR. period /= n) call raiseError("creatTable(): period failed")
                    exit Tloop
                else
                    ! istate_shift = istate_state(state_shift)
                    istate_shift = get_istate(state_shift, stateTable, stateN)
                    ! loop invariant: has shifted n times
                    shiftTable(istate_shift) = n
                    state_marked(istate_shift) = .True.
                    iclassTable(istate_shift) = iclass
                end if
            end do Tloop
            ! if (mod(iclass, 10000) == 0) print *, 'iclass = ', iclass
        end do
        classN = iclass

        print "(A, I10)", " classN = ", classN

        if (allocated(period_class)) deallocate(period_class); allocate (period_class(classN))
        period_class = period_class_tmp(1:classN)

        if (allocated(istate_class)) deallocate(istate_class); allocate (istate_class(classN))
        istate_class = istate_class_tmp(1:classN)
    END SUBROUTINE createTranslationClass


    !=========================================================================!
    ! create zero momentum basis with **POSITIVE** parity
    !
    ! dependent variables:
    ! ---------------------
    ! shiftTable, iclassTable, classN, istate_class, period_class
    !=========================================================================!
    SUBROUTINE createZeroPiMomentumBasis()
        integer :: istate_parity, shift
        integer(kind=INT64) :: state, state_shift, state_parity, weight(L)
        logical, dimension(classN) :: class_marked, has_parity_basis_tmp
        integer :: k, iclass, iclass_parity, ibasis
        integer, dimension(classN) :: iclass_basis_tmp
        integer, dimension(classN) :: shift_basis_tmp
        double precision, dimension(classN) :: norm_basis_tmp
        integer, parameter  :: NOTHING  = -1

        !----------------------------------------------------------------------
        ! deal with parity symmetry, get indice of |k_a*^p> where a* is the
        ! starting state of class [a], p = ± is the parity
        !----------------------------------------------------------------------
        if (allocated(ibasis_class)) deallocate(ibasis_class)
        allocate (ibasis_class(classN))

        if (allocated(iclass_parity_class)) deallocate(iclass_parity_class)
        allocate (iclass_parity_class(classN))

        if (allocated(state_all_class)) deallocate(state_all_class)
        allocate (state_all_class(0:L, classN))

        ibasis_class = NOTHING
        state_all_class = NOTHING
        iclass_parity_class = NOTHING

        class_marked = .False.
        shift_basis_tmp = NOTHING
        norm_basis_tmp = -1.d0
        iclass_basis_tmp = NOTHING
        has_parity_basis_tmp = .False.

        ! construct the weight table, note the sequence and convention, use to
        ! calculate the parity state
        weight = [(2_INT64**(L - k), k=1, L)]

        !----------------------------------------------------------------------
        ! get state_all_class(:)
        !----------------------------------------------------------------------
        do iclass = 1, classN
            state = stateTable(istate_class(iclass))
            state_all_class(0, iclass) = period_class(iclass)
            state_shift = state
            do k = 1, period_class(iclass)
                ! has shifted k-1 times, k-th element
                state_all_class(k, iclass) = state_shift
                state_shift = ishftc(state_shift, -1, SIZE=L)
            end do
        end do

        !----------------------------------------------------------------------
        ! get state_all_class(:)
        !----------------------------------------------------------------------
        ibasis = 0
        classloop: do iclass = 1, classN
            if (class_marked(iclass)) cycle classloop
            class_marked(iclass) = .True.
            state = stateTable(istate_class(iclass))

            state_parity = 0
            do k = 1, L
                state_parity = state_parity + ibits(state, k - 1, 1)*weight(k)
            end do
            ! istate_parity = istate_state(state_parity)
            istate_parity = get_istate(state_parity, stateTable, stateN)
            iclass_parity = iclassTable(istate_parity)
            iclass_parity_class(iclass) = iclass_parity
            iclass_parity_class(iclass_parity) = iclass

            shift = shiftTable(istate_parity)
            ! Filter out non-compatible class
            if (iclass_parity == iclass) then
                if (momentum == k_ZERO .AND. inversion  == I_ODD) then
                    cycle
                else if (momentum == k_PI) then
                    if ((inversion == I_EVEN) .AND. (mod(shift, 2) == 1)) cycle
                    if ((inversion == I_ODD)  .AND. (mod(shift, 2) == 0)) cycle
                end if
            end if

            ibasis = ibasis + 1
            iclass_basis_tmp(ibasis) = iclass
            if (iclass_parity == iclass) then
                ibasis_class(iclass) = ibasis
                norm_basis_tmp(ibasis) = sqrt(dble(L*L/period_class(iclass)))
            else
                has_parity_basis_tmp(ibasis) = .True.
                class_marked(iclass_parity) = .True.
                if (period_class(iclass) /= period_class(iclass_parity)) then
                    call raiseError("creatTable(): periodicity a ≠ Pa error")
                end if
                ibasis_class(iclass) = ibasis
                ibasis_class(iclass_parity) = ibasis
                norm_basis_tmp(ibasis) = sqrt(dble(2.d0*L*L/period_class(iclass)))
                shift_basis_tmp(ibasis) = shift
            end if
        end do classloop
        basisN = ibasis
        print "(A, I10)", " basisN = ", basisN

        if (allocated(shift_basis)) deallocate(shift_basis)
        allocate (shift_basis(basisN))
        shift_basis = shift_basis_tmp(1:basisN)

        if (allocated(iclass_basis)) deallocate(iclass_basis)
        allocate (iclass_basis(basisN))
        iclass_basis = iclass_basis_tmp(1:basisN)

        if (allocated(norm_basis)) deallocate(norm_basis)
        allocate (norm_basis(basisN))
        norm_basis = norm_basis_tmp(1:basisN)

        if (allocated(has_parity_basis)) deallocate(has_parity_basis)
        allocate (has_parity_basis(basisN))
        has_parity_basis = has_parity_basis_tmp(1:basisN)

        ! call testing()
    END SUBROUTINE createZeroPiMomentumBasis

    !=========================================================================!
    ! get state_parity of state
    !=========================================================================!
    INTEGER(kind=INT64) FUNCTION get_state_parity_state(state)
        integer(kind=INT64), intent(in) :: state
        integer :: k
        integer(kind=INT64), SAVE :: weight(L) = [(2_INT64**(L - k), k=1, L)]

        get_state_parity_state = 0
        do k = 1, L
            get_state_parity_state = get_state_parity_state + IBITS(state, k - 1, 1)*weight(k)
        end do
    END FUNCTION get_state_parity_state


    !=========================================================================!
    ! get the istate of state using state search in ordered stateTable
    ! Todo:
    !  * Pure table consume too much memeory 2**L, try out hash table later
    !  * modify interface and pack the mature version into mod_service.f90,
    !  name it as indexBinarySearch(value, table, iMax), so will work as
    !  long as table is ordered (check the order first)
    !=========================================================================!
    INTEGER FUNCTION get_istate(state, stateTable, stateN)
        integer(kind=INT64), intent(in) :: state
        integer(kind=INT64), intent(in) :: stateN
        integer(kind=INT64), intent(in), dimension(stateN) :: stateTable
        integer :: istateMin, istateMax, istate_new

        ! perform state exclusive search so that (imin + imax)/2 is always
        ! different from imin and imax
        if (state == stateTable(1)) then
            get_istate = 1
            return
        else
            istateMin = 1
        end if

        if (state == stateTable(stateN)) then
            get_istate = int(stateN, KIND=4)
            return
        else
            istateMax = int(stateN, KIND=4)
        end if

        do
            istate_new = (istateMin + istateMax)/2
            if (stateTable(istate_new) == state) then
                get_istate = istate_new
                return
            else if (stateTable(istate_new) < state) then
                istateMin = istate_new
            else
                istateMax = istate_new
            end if
        end do
    END FUNCTION get_istate


    !=========================================================================!
    ! constuct Hamiltonian:
    ! m=0 result agree with Papic/eigendecomposition/eigs_periodic_N12_k0_p0.h5
    !=========================================================================!
    SUBROUTINE constructHam(Ham, basisN)
        implicit none
        integer, intent(in) :: basisN
        double precision, allocatable, intent(out) :: Ham(:, :)
        ! i, j is short for ibasis, n is for H = ∑ H_n, sa_ = s_a', sb_ = s_b'
        integer :: i, j, n, sj, sa_, sb_
        integer(KIND=INT64) :: state, state_j
        integer :: iclass, iclass_j, iclass_j_parity, istate_j, istate_j_parity
        ! NOTE:
        ! For convenience, we have use the index from 0:L-1 instead of our
        ! common convention of 1:L. The reason is that mod(, L) takes 0:L-1
        integer :: sigma(0:L-1)
        integer(kind=INT64) :: weight(0:L-1)
        double precision :: factor
        ! for generic k values should be done separately
        integer, dimension(-L:L) :: expList

        weight = [(2_INT64**(L - n), n=1, L)]
        if (momentum == k_ZERO) then
            ! for momentum 0, exp[i k n] = 1.d0
            expList = 1
        else if (momentum == k_PI) then
            ! for momentum π, exp[i π n] = ±1
            expList = [(2*mod(n+2*L+1, 2) - 1, n=-L, L)]
        end if

        if (allocated(Ham)) deallocate(Ham); allocate (Ham(basisN, basisN))

        ! H_ji = <bra| H | ket> where j, i index bra and ket
        Ham = 0.d0
        !---------------------------------------------------------------------
        ! site :   0 1 2    ... L-2 L-1
        !          _________________
        !         |_|_|_|_| ... |_|_|
        !
        ! bit :  high weight     low weight
        !---------------------------------------------------------------------
        do i = 1, basisN
            iclass = iclass_basis(i)
            state = stateTable(istate_class(iclass))
            ! construc the sigma(0:L-1) array
            sigma = [(int(ibits(state, L-n-1, 1)), n=0, L - 1)]
            !     diagonal -m*σ_z term
            Ham(i, i) = Ham(i, i) - mass*(2*popcnt(state) - L)
            ! off-diagonal ∑ P_{n-1} X_n P_{n+1} term

            H_nLoop: do n = 0, L - 1
                if (sigma(mod(n - 1 + L, L)) /= 0 .or. sigma(mod(n + 1, L)) /= 0) cycle H_nLoop
                state_j = ieor(state, weight(n))
                istate_j = istate_state(state_j)
                ! istate_j = get_istate(state_j, stateTable, stateN)
                if (istate_j == NOTHING) then
                    ! Here physical means obey PXP restricted Hilbert space rule
                    call raiseError("constructHam(): istate_j = NOTHING for physical state_j")
                end if
                istate_j_parity = istate_state(get_state_parity_state(state_j))
                ! istate_j_parity = get_istate(get_state_parity_state(state_j), stateTable, stateN)

                iclass_j = iclassTable(istate_j)
                ! if istate_j is excluded due to momentum k requirement
                if (iclass_j == NOTHING) cycle

                iclass_j_parity = iclass_parity_class(iclass_j)
                if (iclass_j_parity == NOTHING) then
                    ! [Pa*] and [a*] have the same period ==> the same momentum k requirement
                    call raiseError("constructHam(): lonely parity class")
                end if

                j = ibasis_class(iclass_j)
                ! iclass_j is excluded due to parity I requirement
                if (j == NOTHING) cycle

                if (has_parity_basis(j)) then
                    if (iclass_j == iclass_j_parity) call raiseError("constructHam(): has_parity failed")
                    sj = shift_basis(j)
                    if (sj == NOTHING) call raiseError("constructHam(): sj == NOTHING")
                end if

                sa_ = shiftTable(istate_j)
                sb_ = shiftTable(istate_j_parity)
                if (sa_ == NOTHING .OR. sb_ == NOTHING) then
                    call raiseError("constructHam(): shift = nothing")
                end if

                ! common factor here
                factor = norm_basis(j)/norm_basis(i)
                if (iclass == iclass_parity_class(iclass)) then
                    !--------------------------------------------------------------
                    ! if [a] = [Pa]
                    !--------------------------------------------------------------
                    if (iclass_j == iclass_j_parity) then
                        ! if [a'] = [Pa']
                        Ham(j, i) = Ham(j, i) + factor*expList(shiftTable(istate_j))
                    else if (iclass_j < iclass_j_parity) then
                        ! a'* < (Pa')*
                        Ham(j, i) = Ham(j, i) + 0.5d0*factor*expList(sa_)
                    else
                        ! a'* > (Pa')*
                        Ham(j, i) = Ham(j, i) + inversion*0.5d0*factor*expList(sa_ - sj)
                    end if
                else
                    !--------------------------------------------------------------
                    ! [a] ≠ [Pa]
                    !--------------------------------------------------------------
                    if (iclass_j == iclass_j_parity) then
                        ! if [a'] = [Pa']
                        Ham(j, i) = Ham(j, i) + factor*(expList(sa_) + inversion*expList(sb_))
                    else if (iclass_j < iclass_j_parity) then
                        Ham(j, i) = Ham(j, i) + 0.5d0*factor*(expList(sa_) + expList(sb_ - sj))
                    else
                        Ham(j, i) = Ham(j, i) + inversion*0.5d0*factor*(expList(sb_) + expList(sa_ - sj))
                    end if
                end if
            end do H_nLoop
        end do
        print "(A)", "Subroutine constructHam() done."
    END SUBROUTINE constructHam


    !=========================================================================!
    ! translate wave function in momentum space to wave function in real-basis
    !  == Limitted to 0, pi momentum space since psi is taken to be real ==
    !=========================================================================!
    SUBROUTINE translatePsi(vect, psi)
        double precision, intent(in) :: vect(1:basisN)
        double precision, intent(out) :: psi(1:stateN)
        integer, dimension(0:L-1) :: expList
        integer :: j, n, istate, iclass, period, ibasis, shift
        double precision :: coeff

        !----------------------------------------------------------------------
        ! |k> = 1/norm Σ_{m=0}^{m=L-1} exp(-ikm) T^m |a>
        !     = 1/norm Σ_{m=0}^{m=t-1} exp(-ikm) T^m |a> * L/period
        !----------------------------------------------------------------------
        if (momentum == k_ZERO) then
            ! for momentum k=0, exp[-i k n] = 1.d0
            expList = 1
        else if (momentum == k_PI) then
            ! for momentum k=π, exp[-i k n] = ±1
            expList = [(2*mod(n+1, 2) - 1, n=0, L-1)]
        else
            call raiseError("translatePsi(): momentum k != 0, pi not done yet")
        end if

        psi = 0.d0
        do ibasis = 1, basisN
            iclass = iclass_basis(ibasis)
            period = period_class(iclass)
            coeff = 1.0/norm_basis(ibasis)*L/period*vect(ibasis)
            do j = 0, period-1
                istate = istate_state(state_all_class(j+1, iclass))
                if (istate == NOTHING) stop "nothing"
                psi(istate) = psi(istate) + expList(j)*coeff
            end do

            if (has_parity_basis(ibasis)) then
                iclass = iclass_parity_class(iclass)
                shift = shift_basis(ibasis)
                do j = 0, period-1
                    istate = istate_state(state_all_class(j+1, iclass))
                    if (istate == NOTHING) stop "nothing"
                    ! only for momentum 0 and π
                    psi(istate) = psi(istate) + inversion*expList(j)*expList(shift)*coeff
                end do
            end if
        end do
    END SUBROUTINE translatePsi

    !=========================================================================!
    ! get the order parameter of staggered magnetization
    !=========================================================================!
    DOUBLE PRECISION FUNCTION get_sz(psi) result(res)
        double precision, intent(in) :: psi(1:basisN)
        integer :: period, ibasis, iclass
        integer(kind=INT64) :: state
        ! s_z for single component |i> in an eigenstate average, sz_ave for
        ! s_z averaged over all component with coefficient c_i (real)
        double precision :: sz, sz_ave

        sz_ave = 0.d0
        do ibasis = 1, basisN
            iclass = iclass_basis(ibasis)
            period = int(state_all_class(0, iclass))
            state  = state_all_class(1, iclass)
            sz = (2.d0*POPCNT(state) - L)/dble(L)
            sz_ave = sz_ave + sz*psi(ibasis)**2
        end do
        res = sz_ave
    END FUNCTION get_sz

    SUBROUTINE get_szList()
        double precision :: szList(1:basisN)
        integer :: n
        do n = 1, basisN
            szList(n) = get_sz(Ham(:, n))
        end do

        call execute_command_line("mkdir -p ./data")
        open (unit=11, file="./data/L"//str(L)//"_E_sz.dat", status='replace')
            do n = 1, basisN
                write (11, *) eigens(n), szList(n)
            end do
        close (11)
    END SUBROUTINE get_szList

    function str(i) result(strout)
        integer, intent(in) :: i
        character(:), allocatable :: strout
        character(range(i)+2) :: strtmp
        write(strtmp,'(I0)') i
        strout = trim(strtmp)
    end function

    !=========================================================================!
    ! Get eigensystem of real symmetry matrix A and A is destroyed on exit
    !     ################################################################
    !     ###############  Note A is destroyed on exit  ##################
    !     ################################################################
    !=========================================================================!
    SUBROUTINE getDSYEigenSystem(A, eigens, errorBD)
        implicit none
        double precision, dimension(:, :), intent(inout):: A
        double precision, dimension(:), intent(out):: eigens
        double precision, optional, intent(out):: errorBD
        double precision, dimension(:), allocatable:: WORK
        integer:: N, LDA, LWORK, INFO, ierror

        N = size(A, 2)
        LDA = size(A, 1)
        if (N /= LDA) then
            call raiseError("getDSYEigenSystem(): A not square matrix")
        end if

        !----------------------------------------------------------------------
        ! http://physics.bu.edu/~py502/lectures4/examples/diatest.f90
        ! Note!!!: The parameter setting is not optimal and unverified.
        !----------------------------------------------------------------------
        ! LWORK = N*(3 + N/2)
        ! allocate (WORK(N*(3 + N/2)))

        !----------------------------------------------------------------------
        ! Query the optimal workspace and diagonalize, this step is very fast
        ! and the time consumed should not bother us.
        !----------------------------------------------------------------------
        LWORK = -1
        ! WORK has to be allocated to match argument type, size may be arbitray
        allocate(WORK(3*N - 1))
        CALL DSYEV('V', 'U', N, A, LDA, eigens, WORK, LWORK, INFO)
        if (INT(WORK(1)) < 3*N - 1) then
            print "(A)", "getDSYEigenSystem(): optimal LWORK < 3*N - 1 !!!"
        ! else
        !     print "(A, I10)", "optimal LWORK = ", INT(WORK(1))
        end if
        LWORK = MAX(3*N-1, INT(WORK(1)))
        deallocate(WORK)
        allocate (WORK(LWORK), stat=ierror)
        if (ierror /= 0) call raiseError("getDSYEigenSystem(): allocation of WORK(:) failed")

        ! ---------------------------------------------------------------------
        !www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442
        !c43fca5493590f8f26cf42fed4044.html#ga442c43fca5493590f8f26cf42fed4044
        ! ---------------------------------------------------------------------
        ! 'N/V'  : eigenvalue only/eigenvalues and eigenvectors.
        ! 'U/L'  : only upper/lower triangle of A is stored and used;
        ! A(:,:) : of dimension (LDA, N), on exit A(:, n) is n-th eigenvector
        ! eigens : If INFO = 0, contains the eigenvalues in ascending order;
        ! WORK   : On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
        ! LWORK  : size of the array WORK.
        ! INFO   : exit status, success if INFO = 0 ; i-th argument error if -i
        ! ---------------------------------------------------------------------
        call DSYEV('V', 'U', N, A, LDA, eigens, WORK, LWORK, INFO)
        if (INFO /= 0) call raiseError("getDSYEigenSystem(): Lapack DSYES failed")

        if (present(errorBD)) then
            errorBD = 1.11d-16*max(abs(eigens(1)), abs(eigens(N)))
        end if
    END SUBROUTINE getDSYEigenSystem


    !=========================================================================!
    ! Get eigenvalues of real symmetry matrix A and A is destroyed on exit
    !     ################################################################
    !     ###############  Note A is destroyed on exit  ##################
    !     ################################################################
    !=========================================================================!
    SUBROUTINE getDSYEigenValues(A, eigens, errorBD)
        implicit none
        double precision, dimension(:, :), intent(inout):: A
        double precision, dimension(:), intent(out):: eigens
        double precision, optional, intent(out):: errorBD
        double precision, dimension(:), allocatable:: WORK
        integer:: N, LDA, LWORK, INFO, ierror

        N = size(A, 2)
        LDA = size(A, 1)
        if (N /= LDA) then
            call raiseError("getDSYEigenValues(): A not square matrix")
        end if

        !----------------------------------------------------------------------
        ! http://physics.bu.edu/~py502/lectures4/examples/diatest.f90
        ! Note!!!: The parameter setting is not optimal and unverified.
        !----------------------------------------------------------------------
        ! LWORK = N*(3 + N/2)
        ! allocate (WORK(N*(3 + N/2)))

        !----------------------------------------------------------------------
        ! Query the optimal workspace and diagonalize, this step is very fast
        ! and the time consumed should not bother us.
        !----------------------------------------------------------------------
        LWORK = -1
        ! WORK has to be allocated to match argument type, size may be arbitray
        allocate(WORK(3*N - 1))
        call DSYEV('N', 'U', N, A, LDA, eigens, WORK, LWORK, INFO)
        if (INT(WORK(1)) < 3*N - 1) then
            print "(A)", "getDSYEigenValues(): optimal LWORK < 3*N - 1 !!!"
        ! else
        !     print "(A, I10)", "optimal LWORK = ", INT(WORK(1))
        end if
        LWORK = MAX(3*N-1, INT(WORK(1)))
        deallocate(WORK)
        allocate (WORK(LWORK), stat=ierror)
        if (ierror /= 0) call raiseError("getDSYEigenValues(): allocation of WORK(:) failed")

        ! ---------------------------------------------------------------------
        !www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga44
        !2c43fca5493590f8f26cf42fed4044.html#ga442c43fca5493590f8f26cf42fed4044
        ! ---------------------------------------------------------------------
        ! 'N/V'  : eigenvalue only/eigenvalues and eigenvectors.
        ! 'U/L'  : only upper/lower triangle of A is stored and used;
        ! A(:,:) : is of dimension (LDA, N), on exit A is destroyed
        ! eigens : If INFO = 0, contains the eigenvalues in ascending order;
        ! WORK   : On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
        ! LWORK  : size of the array WORK.
        ! INFO   : exit status, success if INFO = 0 ; i-th argument error if -i
        ! ---------------------------------------------------------------------
        call DSYEV('N', 'U', N, A, LDA, eigens, WORK, LWORK, INFO)
        if (INFO /= 0) call raiseError("getDSYEigenValues(): Lapack DSYEV failed")

        if (present(errorBD)) then
            errorBD = 1.11d-16*max(abs(eigens(1)), abs(eigens(N)))
        end if
    END SUBROUTINE getDSYEigenValues


    !=========================================================================!
    ! For fatal errors, create a file ERROR.txt and write a message to it,
    ! then stop the program.
    !=========================================================================!
    SUBROUTINE raiseError(message)
        character(len=*), OPTIONAL, intent(in) :: message
        integer :: ioflag

        ! 'truncate -s 0 ERROR.dat' might be unavailable
        call execute_command_line('> ERROR.txt')

        if(present(message)) then
            print *, message
            open (unit=110, file='ERROR.txt', status='old', IOSTAT=ioflag)
                if (ioflag /= 0) stop 'fail to create ERROR.txt file'
                write(110, '(A)') message
            close (110)
        end if

        stop
    END SUBROUTINE raiseError

    !=========================================================================!
    ! for constant parameter one has to use a C fortran to retrieve the value
    !=========================================================================!
    function get_L() bind(c)
        integer(C_INT) :: get_L
        get_L = int(L, C_INT)
    end function get_L

    function get_mass() bind(c)
        real(C_DOUBLE) :: get_mass
        get_mass = mass
    end function get_mass

    !--------------------------------------------------------------------------
    ! get various state-level tables
    !--------------------------------------------------------------------------
    function get_stateTable() bind(c)
        type(c_ptr) :: get_stateTable
        get_stateTable = c_loc(stateTable)
    end function get_stateTable

    function get_shiftTable() bind(c)
        type(c_ptr) :: get_shiftTable
        get_shiftTable = c_loc(shiftTable)
    end function get_shiftTable

    function get_iclassTable() bind(c)
        type(c_ptr) :: get_iclassTable
        get_iclassTable = c_loc(iclassTable)
    end function get_iclassTable

    ! function get_istate_state() bind(c)
    !     type(c_ptr) :: get_istate_state
    !     get_istate_state = c_loc(istate_state)
    ! end function get_istate_state

    !--------------------------------------------------------------------------
    ! get various class-level tables
    !--------------------------------------------------------------------------
    function get_istate_class() bind(c)
        type(c_ptr) :: get_istate_class
        get_istate_class = c_loc(istate_class)
    end function get_istate_class

    function get_period_class() bind(c)
        type(c_ptr) :: get_period_class
        get_period_class = c_loc(period_class)
    end function get_period_class

    function get_ibasis_class() bind(c)
        type(c_ptr) :: get_ibasis_class
        get_ibasis_class = c_loc(ibasis_class)
    end function get_ibasis_class

    function get_iclass_parity_class() bind(c)
        type(c_ptr) :: get_iclass_parity_class
        get_iclass_parity_class = c_loc(iclass_parity_class)
    end function get_iclass_parity_class

    !--------------------------------------------------------------------------
    ! get various basis-level tables
    !--------------------------------------------------------------------------
    function get_iclass_basis() bind(c)
        type(c_ptr) :: get_iclass_basis
        get_iclass_basis = c_loc(iclass_basis)
    end function get_iclass_basis

    function get_shift_basis() bind(c)
        type(c_ptr) :: get_shift_basis
        get_shift_basis = c_loc(shift_basis)
    end function get_shift_basis

    function get_norm_basis() bind(c)
        type(c_ptr) :: get_norm_basis
        get_norm_basis = c_loc(norm_basis)
    end function get_norm_basis

    function get_has_parity_basis() bind(c)
        type(c_ptr) :: get_has_parity_basis
        get_has_parity_basis = c_loc(has_parity_basis)
    end function get_has_parity_basis

    !--------------------------------------------------------------------------
    ! get the translated eigenvector
    !--------------------------------------------------------------------------
    function get_psi(n) bind(c)
        implicit none
        integer(C_INT) :: n
        type(c_ptr) :: get_psi
        call translatePsi(Ham(:, n), psi)
        get_psi = c_loc(psi)
    end function get_psi

END MODULE symmetryED
