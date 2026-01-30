module ABINIT_PES_WRAPPER
      use PES_WRAPPER
      implicit none

      character(80), parameter :: fparams='H2CO_ABINIT_PES'
end module ABINIT_PES_WRAPPER

subroutine abinit_initialize(flag)
      use ABINIT_PES_WRAPPER
      implicit none
      integer, intent(out) :: flag

      call initialize_from_file(fparams)
      flag = initialized
end subroutine abinit_initialize

subroutine abinit_pot(coords, energy)
      use ABINIT_PES_WRAPPER
      implicit none

      double precision, intent(in) :: coords(6)
      double precision, intent(out) :: energy

      call calcpot(coords, energy)
end subroutine abinit_pot