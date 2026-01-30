module PES_WRAPPER
      use PES_DATA
      
      implicit none

      double precision, private :: param(maxnparams)
      integer, private :: ideg(6,maxnparams), nparams
      integer :: initialized

contains
      subroutine initialize_from_file(fparams)
            implicit none
            character(80), intent(in) :: fparams

            call  check_initialized(initialized, fparams, param, ideg, nparams)
      end subroutine initialize_from_file

      subroutine calcpot(coords, energy)
            implicit none
      
            double precision, intent(in) :: coords(6)
            double precision, intent(out) :: energy

            call calcpot_fitted(coords, energy, param, ideg, nparams)

      end subroutine calcpot

end module PES_WRAPPER