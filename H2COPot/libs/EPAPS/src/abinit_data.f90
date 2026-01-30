
module PES_DATA

implicit none
character(80), parameter    :: fparams='H2CO_ABINIT_PES'
integer, parameter :: maxnparams = 1000
double precision :: param(maxnparams)

SUBROUTINE load_params()

      write(*,'(a)') 'Read ab initio PES/MEP parameters'

      iunit = 1

      open(unit=iunit,action='read',status='old',file=trim(fparams),iostat=info)

      if (info/=0) then
            write(*,'(/a,1x,a,1x,a)') 'File', trim(fparams), 'not found'
            stop
      endif

      iparam = 0

      do

            iparam = iparam + 1

            read(iunit,*,end=10,err=101), label, ideg(1:6,iparam), param(iparam)

            write(*,'(1x,i3,3x,a10,3x,6i3,3x,es16.8)') iparam, trim(label), ideg(1:6,iparam), param(iparam)

            cycle
      10    write(*,'(a)') 'end'; exit

      enddo

      nparams = iparam-1

      close(iunit)
END SUBROUTINE load_params

end module PES_DATA


!======================================================================================================================

SUBROUTINE calcpot(r, Area)

IMPLICIT NONE
REAL, INTENT(IN) :: r
REAL, INTENT(OUT) :: Area

! Declare local constant Pi
REAL, PARAMETER :: Pi = 3.1415927

Area = Pi * r * r

END SUBROUTINE Compute_Area


character(80), parameter    :: fparams='H2CO_ABINIT_PES'
double precision, parameter :: Pi=3.1415926535897932e0
integer, parameter          :: maxnparams = 1000

double precision :: param(maxnparams), point(6), point_(6), cosrho, xieq(6), y(6), f, xi(6)
integer          :: ideg(6,maxnparams), kdeg(6), iparam, nparams, iunit, ipoint, info
character(30)    :: label


!===========================


write(*,'(a,1x,a)') 'Calculated energies (distances in Angstrom, angles in degrees, energy in cm-1)'

ipoint = 0

do

      ipoint = ipoint + 1

      if (ipoint==1) then

            point(1) = param(1)
            point(2) = param(6)
            point(3) = param(6)
            point(4) = param(11)
            point(5) = param(11)
            point(6) = Pi

            point_(1:3) = point(1:3)
            point_(4:6) = point(4:6)*180.d0/Pi

      else

            read(*,*,end=20,err=201) point_(1:6)

            point(1:3) = point_(1:3)
            point(4:6) = point_(4:6)*Pi/180.d0

      endif


      cosrho = cos(point(6)) + 1.0d0

      xieq(1) = sum(param(1:5)*cosrho**ideg(1,1:5))
      xieq(2) = sum(param(6:10)*cosrho**ideg(2,6:10))
      xieq(3) = xieq(2)
      xieq(4) = sum(param(11:15)*cosrho**ideg(4,11:15))
      xieq(5) = xieq(4)

      y(1:3) = 1.0d0-exp(-(point(1:3)-xieq(1:3)))
      y(4:5) = point(4:5)-xieq(4:5)
      y(6)   = cosrho


      f = 0.0d0

      do iparam = 16, nparams

            xi(1:6) = y(1:6)**ideg(1:6,iparam)

            f = f + param(iparam)*product(xi(1:6))

            if (ideg(2,iparam)/=ideg(3,iparam).or.ideg(4,iparam)/=ideg(5,iparam)) then

                  kdeg(1) = ideg(1,iparam)
                  kdeg(2) = ideg(3,iparam)
                  kdeg(3) = ideg(2,iparam)
                  kdeg(4) = ideg(5,iparam)
                  kdeg(5) = ideg(4,iparam)
                  kdeg(6) = ideg(6,iparam)

                  xi(1:6) = y(1:6)**kdeg(1:6)

                  f = f + param(iparam)*product(xi(1:6))

            endif

      enddo

      write(*,'(1x,i6,1x,6(es16.8),3x,es16.8)') ipoint, point_(1:6), f

      cycle
20    write(*,'(a)') 'end'; exit

enddo


stop


!===========================


101 write(*,'(/a,1x,a)') 'Error has occured when read file', trim(fparams); stop
201 write(*,'(/a)')      'Error has occured when read points'; stop


end program pes_h2co
