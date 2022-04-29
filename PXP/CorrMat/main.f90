!===============================================================================
! Description
! -----------
!   The purpose of this simple program is to test the mod_symmetryED.f90 and
!   make sure the initialize() function inside can correctly run.
!
!                    gfortran -Wall main.f90 -llapack -lblas
!
!   Version:  1.0
!   Created:
!    Author:  Zhiyuan Yao, zhiyuan.yao@me.com
! Institute:  Insitute for Advanced Study, Tsinghua University
!===============================================================================
! Change Log:
!   *
!===============================================================================
! todo:
!   *
!===============================================================================
include "./mod_symmetryED.f90"
PROGRAM main
    USE, INTRINSIC:: ISO_FORTRAN_ENV
    USE symmetryED
    IMPLICIT NONE

    call initialize()
END PROGRAM main

