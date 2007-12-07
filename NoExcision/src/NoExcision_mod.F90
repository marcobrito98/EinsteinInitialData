! $Header$

#include "cctk.h"

module NoExcision_mod

  implicit none

  CCTK_REAL, parameter :: zero = 0.0
  integer, parameter :: wp = kind(zero)
  CCTK_REAL, parameter :: one = 1.0_wp
  CCTK_REAL, parameter :: two = 2.0_wp
  CCTK_INT :: infnorm_handle, sum_handle
  CCTK_INT :: sym_selector
  CCTK_INT :: loop_counter
  CCTK_REAL, dimension(16) :: infnormresid
  CCTK_REAL, dimension(16) :: delta_new, delta_old, alpha, beta, absres
  character(len=18), dimension(16) :: red_names = (/ 'NoExcision::redgxx', &
                                                     'NoExcision::redgxy', &
                                                     'NoExcision::redgxz', &
                                                     'NoExcision::redgyy', &
                                                     'NoExcision::redgyz', &
                                                     'NoExcision::redgzz', &
                                                     'NoExcision::redkxx', &
                                                     'NoExcision::redkxy', &
                                                     'NoExcision::redkxz', &
                                                     'NoExcision::redkyy', &
                                                     'NoExcision::redkyz', &
                                                     'NoExcision::redkzz', &
                                                     'NoExcision::red   ', &
                                                     'NoExcision::redx  ', &
                                                     'NoExcision::redy  ', &
                                                     'NoExcision::redz  ' /)
  character(len=5), dimension(16) :: var_names = (/ '  gxx', '  gxy', &
                                                    '  gxz', '  gyy', &
                                                    '  gyz', '  gzz', &
                                                    '  kxx', '  kxy', &
                                                    '  kxz', '  kyy', &
                                                    '  kyz', '  kzz', &
                                                    '  alp', 'betax', &
                                                    'betay', 'betaz' /)
  integer :: nx, ny, nz
  integer :: ii
  logical, dimension(16) :: cont = (/ (.true., ii=1, 16) /)

  contains

    ! This is the basic residual calculation routine...

    subroutine residual (val, mask, res, si, order )
      CCTK_REAL, dimension(:,:,:), intent(in) :: val
      CCTK_INT, dimension(:,:,:), intent(in) :: mask
      CCTK_REAL, dimension(:,:,:), intent(out) :: res
      CCTK_REAL, intent(in) :: si
      CCTK_INT, intent(in) :: order
      CCTK_REAL, parameter :: six = 6.0_wp
      CCTK_REAL, parameter :: sixth = 1.0_wp / six
      CCTK_REAL, parameter :: four = 4.0_wp
      CCTK_REAL, parameter :: eighteen = 18.0_wp
      CCTK_REAL, parameter :: eighteenth = 1.0_wp / eighteen
      CCTK_REAL, parameter :: fifteen = 15.0_wp
      CCTK_REAL, parameter :: sixty = 60.0_wp
      CCTK_REAL, parameter :: sixtieth = 1.0_wp / sixty
      integer :: i, j, k
  
      select case (order)
      case (2) 

!$OMP PARALLEL DO PRIVATE(j,i) SCHEDULE(GUIDED)
        do k = 2, nz-1
          do j = 2, ny-1
            do i = 2, nx-1
    
              if (mask(i,j,k) > 0) then
  
                res(i,j,k) = val(i+1,j,k) + val(i-1,j,k) + val(i,j+1,k) + &
                             val(i,j-1,k) + val(i,j,k+1) + val(i,j,k-1) - &
                             six * val(i,j,k)
                res(i,j,k) = sign(one,si) * res(i,j,k)
  
              end if
    
            end do
          end do
        end do
!$OMP END PARALLEL DO

      case (4)

!$OMP PARALLEL DO PRIVATE(j,i) SCHEDULE(GUIDED)
        do k = 3, nz-2
          do j = 3, ny-2
            do i = 3, nx-2

              if (mask(i,j,k) > 0) then
  
                res(i,j,k) = - ( val(i+2,j,k) + val(i-2,j,k) + &
                                 val(i,j+2,k) + val(i,j-2,k) + &
                                 val(i,j,k+2) + val(i,j,k-2) ) + &
                        four * ( val(i+1,j,k) + val(i-1,j,k) + &
                                 val(i,j+1,k) + val(i,j-1,k) + &
                                 val(i,j,k+1) + val(i,j,k-1) ) - &
                      eighteen * val(i,j,k)
                res(i,j,k) = sign(one,si) * res(i,j,k)

              end if

            end do
          end do
        end do
!$OMP END PARALLEL DO

      case (6)

!$OMP PARALLEL DO PRIVATE(j,i) SCHEDULE(GUIDED)
        do k = 4, nz-3
          do j = 4, ny-3
            do i = 4, nx-3

              if (mask(i,j,k) > 0) then

                res(i,j,k) =   ( val(i+3,j,k) + val(i-3,j,k) + &
                                 val(i,j+3,k) + val(i,j-3,k) + &
                                 val(i,j,k+3) + val(i,j,k-3) ) - &
                         six * ( val(i+2,j,k) + val(i-2,j,k) + &
                                 val(i,j+2,k) + val(i,j-2,k) + &
                                 val(i,j,k+2) + val(i,j,k-2) ) + &
                     fifteen * ( val(i+1,j,k) + val(i-1,j,k) + &
                                 val(i,j+1,k) + val(i,j-1,k) + &
                                 val(i,j,k+1) + val(i,j,k-1) ) - &
                         sixty * val(i,j,k)
                res(i,j,k) = sign(one,si) * res(i,j,k)

              end if

            end do
          end do
        end do
!$OMP END PARALLEL DO

      case default

        call CCTK_WARN ( 0, 'Internal error: Order out of range' )
  
      end select

    end subroutine residual


    ! This is the routine calculating residuals for all variables, that
    ! have not yet converged.

    subroutine residual_all ( v1, v2, v3, v4, v5, v6, v7, v8, &
                               v9, v10, v11, v12, v13, v14, v15, v16, &
                               r1, r2, r3, r4, r5, r6, r7, r8, &
                               r9, r10, r11, r12, r13, r14, r15, r16, &
                               mask, si, order )

      CCTK_REAL, dimension(:,:,:), intent(in) :: v1, v2, v3, v4, v5, v6, &
                                                 v7, v8, v9, v10, v11, v12, &
                                                 v13, v14, v15, v16
      CCTK_REAL, dimension(:,:,:), intent(out) :: r1, r2, r3, r4, r5, r6, &
                                                  r7, r8, r9, r10, r11, r12, &
                                                  r13, r14, r15, r16
      CCTK_INT, dimension(:,:,:), intent(in) :: mask
      CCTK_REAL, intent(in) :: si
      CCTK_INT, intent(in) :: order

      if ( cont(1) ) call residual ( v1, mask, r1, si, order )
      if ( cont(2) ) call residual ( v2, mask, r2, si, order )
      if ( cont(3) ) call residual ( v3, mask, r3, si, order )
      if ( cont(4) ) call residual ( v4, mask, r4, si, order )
      if ( cont(5) ) call residual ( v5, mask, r5, si, order )
      if ( cont(6) ) call residual ( v6, mask, r6, si, order )
      if ( cont(7) ) call residual ( v7, mask, r7, si, order )
      if ( cont(8) ) call residual ( v8, mask, r8, si, order )
      if ( cont(9) ) call residual ( v9, mask, r9, si, order )
      if ( cont(10) ) call residual ( v10, mask, r10, si, order )
      if ( cont(11) ) call residual ( v11, mask, r11, si, order )
      if ( cont(12) ) call residual ( v12, mask, r12, si, order )
      if ( cont(13) ) call residual ( v13, mask, r13, si, order )
      if ( cont(14) ) call residual ( v14, mask, r14, si, order )
      if ( cont(15) ) call residual ( v15, mask, r15, si, order )
      if ( cont(16) ) call residual ( v16, mask, r16, si, order )

    end subroutine residual_all


    ! This routine multiplies two variables. Only does it according
    ! to the mask and only for variables that have not converged yet.

    subroutine multiply ( u1, u2, u3, u4, u5, u6, u7, u8, &
                          u9, u10, u11, u12, u13, u14, u15, u16, &
                          v1, v2, v3, v4, v5, v6, v7, v8, &
                          v9, v10, v11, v12, v13, v14, v15, v16, &
                          r1, r2, r3, r4, r5, r6, r7, r8, &
                          r9, r10, r11, r12, r13, r14, r15, r16, &
                          mask )

      CCTK_REAL, dimension(:,:,:), intent(in) :: u1, u2, u3, u4, u5, u6, &
                                                 u7, u8, u9, u10, u11, u12, &
                                                 u13, u14, u15, u16
      CCTK_REAL, dimension(:,:,:), intent(in) :: v1, v2, v3, v4, v5, v6, &
                                                 v7, v8, v9, v10, v11, v12, &
                                                 v13, v14, v15, v16
      CCTK_REAL, dimension(:,:,:), intent(out) :: r1, r2, r3, r4, r5, r6, &
                                                  r7, r8, r9, r10, r11, r12, &
                                                  r13, r14, r15, r16
      CCTK_INT, dimension(:,:,:), intent(in) :: mask

!$OMP PARALLEL

      if ( cont(1) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r1 = u1*v1
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(2) ) then
!$OMP WORKSHARE
         where ( mask > 0 ) r2 = u2*v2
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(3) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r3 = u3*v3
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(4) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r4 = u4*v4
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(5) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r5 = u5*v5
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(6) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r6 = u6*v6
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(7) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r7 = u7*v7
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(8) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r8 = u8*v8
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(9) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r9 = u9*v9
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(10) ) then
!$OMP WORKSHARE
         where ( mask > 0 ) r10 = u10*v10
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(11) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r11 = u11*v11
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(12) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r12 = u12*v12
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(13) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r13 = u13*v13
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(14) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r14 = u14*v14
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(15) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r15 = u15*v15
!$OMP END WORKSHARE NOWAIT
      end if

      if ( cont(16) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) r16 = u16*v16
!$OMP END WORKSHARE NOWAIT
      end if

!$OMP END PARALLEL 

    end subroutine multiply
                  

    ! This routine takes the product of a scalar and a variable and
    ! adds it to the product of another scalar and variable. Does it
    ! according to the mask and only for variables that have not yet
    ! converged.

    subroutine multiply_sum ( u1, u2, u3, u4, u5, u6, u7, u8, &
                              u9, u10, u11, u12, u13, u14, u15, u16, &
                              v1, v2, v3, v4, v5, v6, v7, v8, &
                              v9, v10, v11, v12, v13, v14, v15, v16, &
                              c1, c2, mask )

      CCTK_REAL, dimension(:,:,:), intent(inout) :: u1, u2, u3, u4, u5, u6, &
                                                    u7, u8, u9, u10, u11, u12, &
                                                    u13, u14, u15, u16
      CCTK_REAL, dimension(:,:,:), intent(in) :: v1, v2, v3, v4, v5, v6, &
                                                 v7, v8, v9, v10, v11, v12, &
                                                 v13, v14, v15, v16
      CCTK_REAL, dimension(16), intent(in) :: c1, c2
      CCTK_INT, dimension(:,:,:), intent(in) :: mask

      
!$OMP PARALLEL

      if ( cont(1) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u1 = c1(1)*u1 + c2(1)*v1
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(2) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u2 = c1(2)*u2 + c2(2)*v2
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(3) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u3 = c1(3)*u3 + c2(3)*v3
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(4) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u4 = c1(4)*u4 + c2(4)*v4
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(5) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u5 = c1(5)*u5 + c2(5)*v5
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(6) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u6 = c1(6)*u6 + c2(6)*v6
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(7) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u7 = c1(7)*u7 + c2(7)*v7
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(8) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u8 = c1(8)*u8 + c2(8)*v8
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(9) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u9 = c1(9)*u9 + c2(9)*v9
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(10) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u10 = c1(10)*u10 + c2(10)*v10
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(11) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u11 = c1(11)*u11 + c2(11)*v11
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(12) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u12 = c1(12)*u12 + c2(12)*v12
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(13) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u13 = c1(13)*u13 + c2(13)*v13
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(14) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u14 = c1(14)*u14 + c2(14)*v14
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(15) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u15 = c1(15)*u15 + c2(15)*v15
!$OMP END WORKSHARE NOWAIT
      endif

      if ( cont(16) ) then
!$OMP WORKSHARE
        where ( mask > 0 ) u16 = c1(16)*u16 + c2(16)*v16
!$OMP END WORKSHARE NOWAIT
      endif

!$OMP END PARALLEL 

    end subroutine multiply_sum


end module NoExcision_mod
