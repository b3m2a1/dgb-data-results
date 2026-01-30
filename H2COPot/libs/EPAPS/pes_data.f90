module PES_DATA
      implicit none
      integer, parameter :: maxnparams = 1000

      contains 
      subroutine check_initialized(initialized, fparams, param, ideg, nparams)
            implicit none

            integer, intent(inout):: initialized
            character(80), intent(in) :: fparams
            double precision, intent(out) :: param(maxnparams)
            integer, intent(out) :: ideg(6,maxnparams), nparams

            if (initialized/=1) then
                  call load_params(fparams, param, ideg, nparams)
                  initialized = 1
            endif
      end subroutine check_initialized

      subroutine load_params(fparams, param, ideg, nparams)
            implicit none

            character(80), intent(in) :: fparams
            double precision, intent(out) :: param(maxnparams)
            integer, intent(out) :: ideg(6,maxnparams), nparams

            ! just a slight repackaging of what had been in the previous version

            integer          :: iparam, iunit, info, status
            character(30)    :: label
            character(80)    :: errmsg = "No error caught"

            write(*,'(a)') 'Read PES/MEP parameters'

            iunit = 1

            open(unit=iunit,action='read',status='old',file=trim(fparams),iostat=info)

            if (info/=0) then
                  write(*,'(/a,1x,a,1x,a)') 'File', trim(fparams), 'not found'
                  stop
            endif

            iparam = 0
            do
                  iparam = iparam + 1

                  read(iunit,*,iostat=status,iomsg=errmsg) label, ideg(1:6,iparam), param(iparam)

                  if (status.gt.0) then
                        write(*,'(/a,1x,a)') 'Error: ', errmsg, 'in reading file', trim(fparams); stop
                  else if (status.lt.0) then 
                        write(*,'(a)') errmsg;
                        exit
                  endif

                  write(*,'(1x,i3,3x,a10,3x,6i3,3x,es16.8)') iparam, trim(label), ideg(1:6,iparam), param(iparam)

            enddo

            nparams = iparam-1

            write(*,'(a,i3,a)') "Read ", nparams, " parameters"

            close(iunit)
      end subroutine load_params

      subroutine calcpot_fitted(coords, energy, param, ideg, nparams)
            implicit none
      
            double precision, intent(in) :: coords(6)
            double precision, intent(in) :: param(maxnparams)
            integer, intent(in) :: ideg(6,maxnparams), nparams
            double precision, intent(out) :: energy
      
            double precision :: cosrho, xieq(6), y(6), xi(6)
            integer          :: kdeg(6), iparam
      
            cosrho = cos(coords(6)) + 1.0d0
      
            xieq(1) = sum(param(1:5)*cosrho**ideg(1,1:5))
            xieq(2) = sum(param(6:10)*cosrho**ideg(2,6:10))
            xieq(3) = xieq(2)
            xieq(4) = sum(param(11:15)*cosrho**ideg(4,11:15))
            xieq(5) = xieq(4)
      
            y(1:3) = 1.0d0-exp(-(coords(1:3)-xieq(1:3)))
            y(4:5) = coords(4:5)-xieq(4:5)
            y(6)   = cosrho
      
            energy = 0.0d0
      
            do iparam = 16, nparams
      
                  xi(1:6) = y(1:6)**ideg(1:6,iparam)
      
                  energy = energy + param(iparam)*product(xi(1:6))
      
                  if (ideg(2,iparam)/=ideg(3,iparam).or.ideg(4,iparam)/=ideg(5,iparam)) then
      
                        kdeg(1) = ideg(1,iparam)
                        kdeg(2) = ideg(3,iparam)
                        kdeg(3) = ideg(2,iparam)
                        kdeg(4) = ideg(5,iparam)
                        kdeg(5) = ideg(4,iparam)
                        kdeg(6) = ideg(6,iparam)
      
                        xi(1:6) = y(1:6)**kdeg(1:6)
      
                        energy = energy + param(iparam)*product(xi(1:6))
      
                  endif
      
            enddo
      
      end subroutine calcpot_fitted

end module PES_DATA