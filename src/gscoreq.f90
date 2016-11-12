!
! Program Gscore-Q
!

program gscoreq

  use file_operations
  implicit none
  integer :: ioerr, npdb, npairs, fdomain, ldomain, ndomain, &
             ires, iat, jat, nres, ipdb, jpdb, ipair, i, icount
  integer, allocatable :: ncontacts(:), fcontact(:), nextcontact(:,:)
  real :: dcontact, contact_square, d, ccut
  real, allocatable :: x(:,:), correlation(:,:), degree(:)
  character(len=200), allocatable :: pdb(:)
  character(len=200) :: pdblist, record, compactlog, gscorelog, format
  logical, allocatable :: contact(:,:)

  ! Input parameters

  pdblist='small.txt'
  fdomain = 55
  ldomain = 306
  ndomain = ldomain-fdomain+1
  npairs = ndomain*(ndomain-1)/2

  dcontact = 8.d0
  contact_square = dcontact*dcontact
  
  ccut = 0.5

  compactlog='compactlog.dat'
  gscorelog='gscoreq.dat'

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
        end if
      end do
      close(20)
    end if
  end do
  write(*,"(a,i8)") '# Number of PDB files: ', npdb
  write(*,"(a,i8)") '# Number of residues in structure: ', nres
  write(*,"(a,i8)") '# Number of residues in domain: ', ndomain
  rewind(10)
  allocate(x(ndomain,3),pdb(npdb),contact(npdb,npairs))

  ! Computing contact vectors

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

  ! Computing, for each model, the fraction of structures which have a correlation
  ! greater than ccut (the degree)

  allocate(degree(npdb))
 
  write(*,"(a)") "# Computing the gscores ... "
  do ipdb = 1, npdb
    degree(ipdb) = 1.e0
  end do
  do ipdb = 1, npdb - 1
    call progress(ipdb,1,npdb)
    do jpdb = ipdb + 1, npdb
      if ( correlation(ipdb,jpdb) > ccut ) then
        degree(ipdb) = degree(ipdb) + 1.e0
        degree(jpdb) = degree(jpdb) + 1.e0
      end if
    end do
  end do
  do ipdb = 1, npdb
    degree(ipdb) = degree(ipdb) / npdb
  end do
  call progress(npdb,1,npdb)

  !
  ! Writting results
  !

  ! Compact log file:

  open(10,file=compactlog)
  write(10,"(a)") '# This a compact contact-correlation file'
  write(10,"(a)") '# Computed with GscoreQ'
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

  ! Gscore file:

  open(10,file='gscoreq.dat')
  write(10,"(a)") "# Output of GscoreQ"
  write(10,"(a,a)") "# Associated compact log file: ", trim(adjustl(compactlog))
  write(10,"(a)") "# Score type: Contact-correlation"
  write(10,"(a,e12.5)") "# Score cutoff:", ccut
  write(10,"(a)") "#"
  write(10,"(a)") "#    G-score  Model"
  do ipdb = 1, npdb
    write(10,"(f12.5,tr2,a)") degree(ipdb), trim(adjustl(pdb(ipdb)))
  end do
  close(10)

end program gscoreq




