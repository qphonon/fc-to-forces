! 
!     Copyright (C) 2021  Chee Kwan Gan (ihpcganck@gmail.com)
!     Copyright (C) 2020  Chee Kwan Gan (ihpcganck@gmail.com)
! 
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
! 

module extra
  use commod
  implicit none

  integer,parameter :: pu           =11
  integer,parameter :: su           =12
  integer,parameter :: sposcaru     =13
  integer,parameter :: fconstu      =14
  integer,parameter :: forcesu      =15

  character(len=*),parameter :: sposcarfile='sposcar.vasp'

  character(len=*),parameter :: pposcarfile='p.vasp'
  character(len=*),parameter :: supercposcarfile='supercell.vasp'

  character(len=*),parameter :: fconstfile='FORCE_CONSTANTS'

  character(len=*),parameter :: forcesfile='forces.dat'

  type fblk
    real(double) :: m(3,3)
  end type fblk

end module extra

program sample
  use extra
  implicit none
  integer :: i,j,ns1,ns2
  type(supercell) :: pposcar,supercposcar,sposcar
  character(len=DSL) :: filestem,inputformat
  integer :: ati,atj
  type(fblk),allocatable :: m(:,:)
  real(double) :: frac1(3),pos1(3),frac2(3),pos2(3),diff
  integer,allocatable :: m13(:),m23(:)
  integer :: ix
  character(len=DSL) :: abspposcarfile
  character(len=DSL) :: abssupercposcarfile
  character(len=DSL) :: abssposcarfile
  character(len=DSL) :: absfconstfile
  character(len=DSL) :: absforcesfile
  character(len=DSL) :: dir
  real(double) :: amplitude
  integer,parameter :: CD=1
  integer,parameter :: FD=2
  integer :: FCD

  write(*,'(A)') './b.out dir difference-scheme amplitude'
  write(*,'(A)') ' Required files: '//trim(sposcarfile)//', '//trim(pposcarfile)//', '//trim(supercposcarfile)

  call fargn(3)

  call fargv(1,dir)

  call fargv(2,FCD)

  call fargv(3,amplitude)

  abspposcarfile = trim(dir)//'/'//trim(pposcarfile)

  abssupercposcarfile = trim(dir)//'/'//trim(supercposcarfile)

  abssposcarfile = trim(dir)//'/'//trim(sposcarfile)

  absfconstfile = trim(dir)//'/'//trim(fconstfile)

  absforcesfile = trim(dir)//'/'//trim(forcesfile)

  write(*,*) ' --------------------'
  write(*,*) 'Reading pposcar...'
  call read_struc(pu,trim(abspposcarfile),filestem,inputformat,pposcar)
  write(*,*)

  write(*,*) ' --------------------'
  write(*,*) 'Reading supercposcar...'
  call read_struc(su,trim(abssupercposcarfile),filestem,inputformat,supercposcar)
  write(*,*)

  write(*,*) ' --------------------'
  write(*,*) 'Reading sposcar...'
  call read_struc(sposcaru,trim(abssposcarfile),filestem,inputformat,sposcar)
  write(*,*)

  open(unit=fconstu,file=trim(absfconstfile),action='read')

  read(fconstu,*) ns2

  write(*,*) 'ns2,sposcar%n,supercposcar%n=',ns2,sposcar%n,supercposcar%n
  if(ns2 /= sposcar%n .or. ns2 /= supercposcar%n) then
    write(*,'(A)') 'err: wrong number of atoms'
    stop 1
  endif
  ns1 = pposcar%n
  write(*,*) 'number of atoms in pposcar = ns1=',ns1
  write(*,*) 'number of atoms in supercellposcar = ns2=',ns2

  allocate(m(ns2,ns2))
  do i = 1, ns2
    do j = 1, ns2
      read(fconstu,*) ati,atj
      if(i /= ati .and. j /= atj) then
        write(*,'(A)') 'err: index problem'
        stop 1
      endif
      read(fconstu,*)  m(i,j)%m(1,1:3)
      read(fconstu,*)  m(i,j)%m(2,1:3)
      read(fconstu,*)  m(i,j)%m(3,1:3)
    enddo
  enddo
  close(fconstu)

  allocate(m13(ns1))
  allocate(m23(ns2))
  m13(1:ns1) = -1

  do i = 1, ns1
    frac1 = pposcar%at(i)%f
    call frac2abs(pposcar%a,frac1,pos1)
    do j = 1, ns2
      frac2 = sposcar%at(j)%f
      call frac2abs(sposcar%a,frac2,pos2)
      diff = vecmag3(pos1-pos2)
      if(diff < 1.0d-8) then
        if(pposcar%at(i)%z /= sposcar%at(j)%z) then
          write(*,*) 'i,j,pposcar%at(i)%z,sposcar%at(j)%z=',i,j,pposcar%at(i)%z,sposcar%at(j)%z
          write(*,'(A)') 'err: Serious identity problem.'
          stop 1
        endif

        m13(i) = j
        exit
      endif
    enddo
    if(m13(i) == -1) then
      write(*,'(A)') 'err: m13 error.'
      stop 1
    endif
  enddo
  write(*,*)
  do i = 1, ns1
    write(*,*) 'i,m13(i)=',i,m13(i)
  enddo
  write(*,*)
  m23(1:ns2) = -1
  do i = 1, ns2
    frac1 = supercposcar%at(i)%f
    call frac2abs(supercposcar%a,frac1,pos1)
    do j = 1, ns2
      frac2 = sposcar%at(j)%f
      call frac2abs(sposcar%a,frac2,pos2)
      diff = vecmag3(pos1-pos2)
      if(diff < 1.0d-8) then
        if(supercposcar%at(i)%z /= sposcar%at(j)%z) then
          write(*,*) 'i,j,supercposcar%at(i)%z,sposcar%at(j)%z=',i,j,supercposcar%at(i)%z,sposcar%at(j)%z
          write(*,'(A)') 'err: Serious identity problem.'
          stop 1
        endif
        m23(i) = j
      endif
    enddo
    if(m23(i) == -1) then
      write(*,'(A)') 'err: m23 error.'
      stop 1
    endif
  enddo
  do i = 1, ns2

  enddo

  open(unit=forcesu,file=trim(absforcesfile),status='replace')

  do i = 1, ns2
    write(forcesu,*) zero,zero,zero
  enddo

  do i = 1, ns1
    do ix = 1, 3
      do j = 1, ns2

        write(forcesu,*) -m(m13(i),m23(j))%m(ix,1:3)*amplitude
      enddo
      do j = 1, ns2
        if(FCD == CD) then
          write(forcesu,*) m(m13(i),m23(j))%m(ix,1:3)*amplitude
        else if(FCD == FD) then
          write(forcesu,*) 0.0, 0.0, 0.0
        endif
      enddo
    enddo
  enddo
  close(forcesu)
  write(*,*)
  write(*,'(A)') trim(absforcesfile)//' is produced'
  write(*,*)
  write(*,'(A,F20.10)') 'amplitude for displacement is ',amplitude
  write(*,*)
  if(FCD == CD) then
    write(*,'(A)') 'We can use option 1 (for central diff) in qphonon.par with forces.dat'
  else if(FCD == FD) then
    write(*,'(A)') 'We can use option 2 (for forward diff) in qphonon.par with forces.dat'
  endif
  write(*,*)

end program sample
