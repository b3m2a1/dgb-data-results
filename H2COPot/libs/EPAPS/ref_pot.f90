module REF_PES_WRAPPER
      use PES_WRAPPER

      implicit none
      character(80), parameter :: fparams='H2CO_REF_PES_2011'

end module REF_PES_WRAPPER

subroutine ref_initialize(flag)
      use REF_PES_WRAPPER
      implicit none
      integer, intent(out) :: flag

      call initialize_from_file(fparams)
      flag = initialized
end subroutine ref_initialize

subroutine ref_pot(coords, energy)
      use REF_PES_WRAPPER
      implicit none

      double precision, intent(in) :: coords(6)
      double precision, intent(out) :: energy
      
      call calcpot(coords, energy)
end subroutine ref_pot