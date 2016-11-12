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

  use file_operations
  implicit none
  integer :: narg
  integer :: ioerr, npdb, npairs, fdomain, ldomain, ndomain, &
             ires, iat, jat, nres, ipdb, jpdb, ipair, i, icount
  integer, allocatable :: ncontacts(:), fcontact(:), nextcontact(:,:)
  real :: dcontact, contact_square, d 
  real, allocatable :: x(:,:), correlation(:,:)
  character(len=200), allocatable :: pdb(:)
  character(len=200) :: pdblist, record, output, format
  logical, allocatable :: contact(:,:)


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

  ! Get number of PDB files and number of residues from input file

  open(10,file=pdblist,status='old',action='read',iostat=ioerr)
  npdb = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    npdb = npdb + 1
    if ( npdb == 1 ) then
      nres = 0
      open(20,file=record,iostat=ioerr)
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
          if ( fdomain == 0 .and. ldomain == 0 ) then
            read(record(23:26),*,iostat=ioerr) ires
            if ( nres == 1 ) fdomain = ires
            ldomain = ires
          end if
        end if
      end do
      close(20)
    end if
  end do
  write(*,"(a,i8)") '# Number of PDB files: ', npdb
  write(*,"(a,i8)") '# Number of residues in structure: ', nres
  write(*,"(a,i8)") '# First residue in domain: ', fdomain
  write(*,"(a,i8)") '# Last residue in domain: ', ldomain
  ndomain = ldomain-fdomain+1
  npairs = ndomain*(ndomain-1)/2
  write(*,"(a,i8)") '# Number of residues in domain: ', ndomain
  allocate(x(ndomain,3),pdb(npdb),contact(npdb,npairs))

  ! Computing contact vectors

  rewind(10)
  write(*,"(a)") "# Computing the contact vector for all structures ... "
  ipdb = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    ipdb = ipdb + 1
    call progress(ipdb,1,npdb)
    pdb(ipdb) = record

    ! Reading coordinates of this PDB file

    open(20,file=pdb(ipdb),iostat=ioerr)
    iat = 0
    do
      read(20,"(a200)",iostat=ioerr) record 
      if ( ioerr /= 0 ) exit
      if ( record(1:4) == "ATOM" .and. &
           trim(adjustl(record(13:16))) == "CA" ) then
        read(record(23:26),*,iostat=ioerr) ires
        if ( ioerr /= 0 ) then
          write(*,*) ' ERROR: Could not read residue number for a CA atom. '
          write(*,*) '        File: ', trim(adjustl(pdb(ipdb)))
          write(*,*) '        Line: ', trim(adjustl(record))
          close(20)
          close(10)
          stop
        end if
        if ( ires >= fdomain .and. ires <= ldomain ) then
          iat = iat + 1
          read(record(31:38),*,iostat=ioerr) x(iat,1)
          if ( ioerr /= 0 ) then
            write(*,*) ' ERROR: Failed reading atom coordintes in file: ', trim(adjustl(pdb(ipdb)))
            stop
          end if
          read(record(39:46),*,iostat=ioerr) x(iat,2)
          if ( ioerr /= 0 ) then
            write(*,*) ' ERROR: Failed reading atom coordintes in file: ', trim(adjustl(pdb(ipdb)))
            stop
          end if
          read(record(47:54),*,iostat=ioerr) x(iat,3)
          if ( ioerr /= 0 ) then
            write(*,*) ' ERROR: Failed reading atom coordintes in file: ', trim(adjustl(pdb(ipdb)))
            stop
          end if
        end if
      end if
    end do
    close(20)
   
    ! Computing the contact logical vectors

    ipair = 0
    do iat = 1, ndomain - 1
      do jat = iat + 1, ndomain
        ipair = ipair + 1
        d = (x(iat,1) - x(jat,1))**2 + &
            (x(iat,2) - x(jat,2))**2 + &
            (x(iat,3) - x(jat,3))**2
        if ( d <= contact_square ) then
          contact(ipdb,ipair) = .true.
        else
          contact(ipdb,ipair) = .false.
        end if
      end do
    end do

  end do
  close(10)

  ! Setting up the linked vectors for running only on existing contacts

  write(*,"(a)") "# Setting up linked contact vectors ... "
  allocate(fcontact(npdb),nextcontact(npdb,npairs),ncontacts(npdb))
  do ipdb = 1, npdb
    call progress(ipdb,1,npdb)
    ncontacts(ipdb) = 0
    fcontact(ipdb) = 0
    do ipair = 1, npairs
      if ( contact(ipdb,ipair) ) then
        nextcontact(ipdb,ipair) = fcontact(ipdb)
        fcontact(ipdb) = ipair
        ncontacts(ipdb) = ncontacts(ipdb) + 1
      end if
    end do
  end do
  
  write(*,"(a)") "# Computing the correlation matrix ... "
  allocate(correlation(npdb,npdb))
  i = 0
  do ipdb = 1, npdb - 1
    do jpdb = ipdb + 1, npdb
      i = i + 1
      call progress(i,1,npdb*(npdb-1)/2)
      icount = 0
      if ( ncontacts(ipdb) <= ncontacts(jpdb) ) then
        ipair = fcontact(ipdb)
        do while( ipair > 0 ) 
          if ( contact(ipdb,ipair) .and. contact(jpdb,ipair) ) icount = icount + 1
          ipair = nextcontact(ipdb,ipair)
        end do
      else
        ipair = fcontact(jpdb)
        do while( ipair > 0 ) 
          if ( contact(ipdb,ipair) .and. contact(jpdb,ipair) ) icount = icount + 1
          ipair = nextcontact(jpdb,ipair)
        end do
      end if
      correlation(ipdb,jpdb) = real(icount)
      correlation(jpdb,ipdb) = real(icount)
    end do
  end do

  !
  ! Writting results
  !

  ! Compact log file:

  open(10,file=output)
  write(10,"(a)") '# This a compact contact-correlation file'
  write(10,"(a)") '# Computed with Qcorrelation'
  write(10,"(a,a)") '# PDB list: ', trim(adjustl(pdblist))
  write(10,"(a,a)") '# Score type: Contact-correlation'
  write(10,*) npdb
  do ipdb = 1, npdb
    write(10,*) ipdb, trim(adjustl(basename(pdb(ipdb))))
  end do
  write(format,*) npdb
  format = "("//trim(adjustl(format))//"(tr1,f8.3))"
  do ipdb = 1, npdb - 1
    write(10,format) (correlation(ipdb,jpdb),jpdb=ipdb+1,npdb)
  end do
  close(10)

end program qmatrix

subroutine argerror()

  write(*,*) ' Run with: gscore-q pdblist.txt qmatrix.dat [dcontact] [fdomain] [ldomain] '
  write(*,*) ' Where: pdblist.txt is the list of PDB files. '
  write(*,*) '        qmatrix.dat is the output correlation matrix. '
  write(*,*) '        dcontact is the distance that defines a contact between CA atoms. '
  write(*,*) '        fdomain and ldomain (optional ) are the number of first and last residues '
  write(*,*) '             to be considered. '
  write(*,*)
  stop

end subroutine argerror
