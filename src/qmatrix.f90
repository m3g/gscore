!
! Q-matrix calculator
!
! Program qmatrix
!
! This program computes the correlation matrix of a series
! of PDB structures, considering the Q-score, which is the
! number of contacts shared by both structures.
!
! L. Martinez
! Institute of Chemistry - University of Campinas
! Nov 11, 2016
! http://leandro.iqm.unicamp.br
!

program qmatrix

  use types
  use file_operations
  implicit none
  type contact_list
    integer, allocatable :: i(:), j(:)
  end type
  integer :: narg
  integer :: ioerr, nmodels, fdomain, ldomain, ndomain, &
             ires, iat, jat, nres, imodel, jmodel, ipair, i, icount, ic, jc, &
             maxcontacts, npairs
  real :: dcontact, contact_square, d 
  real, allocatable :: x(:,:,:), correlation(:)
  character(len=200) :: pdblist, record, output, format
  logical, allocatable :: hascontact(:)
  type(model_type), allocatable :: model(:)
  type(contact_list), allocatable :: contact(:)

  write(*,"(a)") "#" 
  write(*,"(a)") "# Contact correlation" 
  call title()
  write(*,"(a)") "# L. Martinez - Institute of Chemistry, University of Campinas" 
  write(*,"(a)") "# http://leandro.iqm.unicamp.br" 
  write(*,"(a)") "#" 

  ! Input parameters

  narg = iargc()
  call getarg(1,pdblist)
  call getarg(2,output)
  call getarg(3,record)
  read(record,*,iostat=ioerr) dcontact
  if ( ioerr /= 0 ) call argerror()
  contact_square = dcontact*dcontact
  if ( narg == 5 ) then
    call getarg(4,record)
    read(record,*,iostat=ioerr) fdomain
    if ( ioerr /= 0 ) call argerror()
    call getarg(5,record)
    read(record,*,iostat=ioerr) ldomain
    if ( ioerr /= 0 ) call argerror()
  else
    fdomain = 0
    ldomain = 0
  end if
  if ( fdomain == 0 .and. ldomain == 0 ) ndomain = 0

  ! Just check if output file exists and if the user wants to overwrite it

  open(10,file=output,status="new",iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,"(a,a,a)",advance="no") &
         '# Warning: File ', trim(adjustl(output)), ' exists. Overwrite it? (Y/N) '
    read(*,*) record
    if ( record /= "Y" ) then
      close(10)
      write(*,*) ' Quitting ... '
      stop
    end if
  end if
  close(10)

  ! Get number of PDB files and number of residues from input file

  open(10,file=pdblist,status='old',action='read',iostat=ioerr)
  nmodels = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    nmodels = nmodels + 1
    if ( nmodels == 1 ) then
      nres = 0
      open(20,file=record,iostat=ioerr,action="read",status="old")
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Could not open PDB file: ', trim(adjustl(record))
        close(10)
        stop
      end if
      do
        read(20,"(a200)",iostat=ioerr) record 
        if ( ioerr /= 0 ) exit
        if ( record(1:4) == "ATOM" .and. &
             trim(adjustl(record(13:16))) == "CA" ) then
          nres = nres + 1
          if ( ndomain == 0 ) then
            read(record(23:26),*,iostat=ioerr) ires
            if ( nres == 1 ) fdomain = ires
            ldomain = ires
          end if
        end if
      end do
      close(20)
    end if
  end do
  write(*,"(a,i8)") '# Number of PDB files: ', nmodels
  write(*,"(a,i8)") '# Number of residues in structure: ', nres
  write(*,"(a,i8)") '# First residue in domain: ', fdomain
  write(*,"(a,i8)") '# Last residue in domain: ', ldomain

  ! Allocate model vector to sort names

  allocate(model(nmodels))
  rewind(10)
  imodel = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    imodel = imodel + 1
    model(imodel)%file = record
    model(imodel)%name = basename(record)
  end do
  close(10)

  !
  ! Sort model names according to string comparisons
  !

  call sort_by_name(nmodels,model)

  ndomain = ldomain-fdomain+1
  write(*,"(a,i8)") '# Number of residues in domain: ', ndomain
  allocate(x(nmodels,ndomain,3))

  ! Reading the coordinates to memory

  write(*,"(a)") "# Reading coordinates ... "
  do imodel = 1, nmodels
    call progress(imodel,1,nmodels)

    ! Reading coordinates of this PDB file

    open(20,file=model(imodel)%file,status="old",action="read",iostat=ioerr)
    iat = 0
    do
      read(20,"(a200)",iostat=ioerr) record 
      if ( ioerr /= 0 ) exit
      if ( record(1:4) == "ATOM" .and. &
           trim(adjustl(record(13:16))) == "CA" ) then
        read(record(23:26),*,iostat=ioerr) ires
        if ( ioerr /= 0 ) then
          write(*,*) ' ERROR: Could not read residue number for a CA atom. '
          write(*,*) '        File: ', trim(adjustl(model(imodel)%name))
          write(*,*) '        Line: ', trim(adjustl(record))
          close(20)
          stop
        end if
        if ( ires >= fdomain .and. ires <= ldomain ) then
          iat = iat + 1
          read(record(31:38),*,iostat=ioerr) x(imodel,iat,1)
          if ( ioerr /= 0 ) then
            write(*,*) ' ERROR: Failed reading atom coordintes in file: ', &
                        trim(adjustl(model(imodel)%name))
            close(20) ; stop
          end if
          read(record(39:46),*,iostat=ioerr) x(imodel,iat,2)
          if ( ioerr /= 0 ) then
            write(*,*) ' ERROR: Failed reading atom coordintes in file: ', &
                        trim(adjustl(model(imodel)%name))
            close(20) ; stop
          end if
          read(record(47:54),*,iostat=ioerr) x(imodel,iat,3)
          if ( ioerr /= 0 ) then
            write(*,*) ' ERROR: Failed reading atom coordintes in file: ', &
                        trim(adjustl(model(imodel)%name))
            close(20) ; stop
          end if
        end if
      end if
    end do
    close(20)
  end do

  ! Computing the number of contacts and maximum number of contacts in a model

  write(*,"(a)") "# Computing number of contacts in each model ... "
  maxcontacts = 0
  do imodel = 1, nmodels
    call progress(imodel,1,nmodels)
    icount = 0
    do iat = 1, ndomain - 1
      do jat = iat + 3, ndomain
        ipair = ipair + 1
        d = (x(imodel,iat,1) - x(imodel,jat,1))**2 + &
            (x(imodel,iat,2) - x(imodel,jat,2))**2 + &
            (x(imodel,iat,3) - x(imodel,jat,3))**2
        if ( d <= contact_square ) icount = icount + 1
      end do
    end do
    model(imodel)%ncontacts = icount
    maxcontacts = max(icount,maxcontacts)
  end do

  write(*,"(a,i10)") "# Maximum number of contacts in a structure: ", maxcontacts
  write(*,"(a,f12.5,a)") "# Estimated memory requirement: ",&
                          8*(real(maxcontacts)/1024)*real(nmodels/1024)/1024," GB"

  allocate(contact(nmodels))
  do imodel = 1, nmodels
    allocate(contact(imodel)%i(model(imodel)%ncontacts),&
             contact(imodel)%j(model(imodel)%ncontacts))
  end do

  ! Compute contact matrices for each structure and store
  
  write(*,"(a)") "# Computing the contact matrix for each model ... "
  do imodel = 1, nmodels
    call progress(imodel,1,nmodels)
    ic = 0
    do iat = 1, ndomain - 1
      do jat = iat + 3, ndomain
        ipair = ipair + 1
        d = (x(imodel,iat,1) - x(imodel,jat,1))**2 + &
            (x(imodel,iat,2) - x(imodel,jat,2))**2 + &
            (x(imodel,iat,3) - x(imodel,jat,3))**2
        if ( d <= contact_square ) then
          ic = ic + 1
          contact(imodel)%i(ic) = iat
          contact(imodel)%j(ic) = jat
        end if
      end do
    end do
  end do
  deallocate(x)

  ! Write output file preamble 

  open(10,file=output,iostat=ioerr)
  write(10,"(a)") '# This a compact contact-correlation file'
  write(10,"(a)") '# Computed with Qcorrelation'
  write(10,"(a,a)") '# PDB list: ', trim(adjustl(pdblist))
  write(10,"(a,f12.5)") '# Score type: Contact-correlation - dcontact = ', dcontact
  write(10,*) nmodels, maxcontacts
  do imodel = 1, nmodels
    write(10,*) imodel, trim(adjustl(basename(model(imodel)%name))), model(imodel)%ncontacts
  end do
  write(format,*) nmodels
  format = "("//trim(adjustl(format))//"(tr1,f10.3))"

  ! Actually computing the correlation matrix and writting output, by line

  write(*,"(a)") "# Computing the correlation matrix and writting output ... "

  i = 0
  npairs = ndomain*(ndomain-1)/2
  allocate(hascontact(npairs))
  allocate(correlation(nmodels))
  do imodel = 1, nmodels - 1

    ! Compute logical contact vector for imodel

    do ic = 1, npairs
      hascontact(ic) = .false.
    end do
    do ic = 1, model(imodel)%ncontacts
      iat = contact(imodel)%i(ic)
      jat = contact(imodel)%j(ic)
      ! Counting for j > i, without the diagonal, thus the lower triangle
      ! and the diagonal must be discounted from the list -(iat*(iat-1)/2+iat)=-(iat*(iat+1)/2)
      jc = (iat-1)*ndomain+jat - iat*(iat+1)/2
      hascontact(jc) = .true.
    end do

    ! Check the number of contacts present in both models
      
    do jmodel = imodel + 1, nmodels
      i = i + 1
      call progress(i,1,nmodels*(nmodels-1)/2)
      icount = 0
      do ic = 1, model(jmodel)%ncontacts
        iat = contact(jmodel)%i(ic)
        jat = contact(jmodel)%j(ic)
        jc = (iat-1)*ndomain+jat - iat*(iat+1)/2 
        if ( hascontact(jc) ) then
          icount = icount + 1
        end if
        correlation(jmodel) = icount
      end do
    end do

    ! Write this line of the output file

    write(10,format) (correlation(jmodel),jmodel=imodel+1,nmodels)
  end do
  close(10)
  
  write(*,"(a)") "#"
  write(*,"(3a)") "# Output file: ", trim(adjustl(output)), " written."
  write(*,"(a)") "#"
  write(*,"(a)") "# END "

end program qmatrix

subroutine argerror()

  write(*,*) ' Run with: qmatrix pdblist.txt qmatrix.dat [dcontact] [fdomain] [ldomain] '
  write(*,*) ' Where: pdblist.txt is the list of PDB files. '
  write(*,*) '        qmatrix.dat is the output correlation matrix. '
  write(*,*) '        dcontact is the distance that defines a contact between CA atoms. '
  write(*,*) '        fdomain and ldomain (optional ) are the number of first and last residues '
  write(*,*) '             to be considered. '
  write(*,*)
  stop

end subroutine argerror
