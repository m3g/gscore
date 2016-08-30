!
! Module with functions to operate on file names and strings
!

module file_operations

  contains

    !
    ! Function that determines the basename of a file,
    ! removing the path and the extension
    !
    
    character(len=200) function basename(filename)
    
      integer :: i
      character(len=200) :: filename
    
      basename = trim(adjustl(filename))
      i = length(basename)
      idot = i+1
      do while(basename(i:i) /= "/")
        if ( basename(i:i) == "." ) then
          idot = i
        end if
        i = i - 1
        if ( i == 0 ) exit
      end do
      i = i + 1
      basename = basename(i:idot-1)
      do i = idot, 200
        basename(i:i) = achar(32)
      end do
    
    end function basename
    
    !
    ! Function that determines the length of a string
    !
    
    integer function length(string)
    
      implicit none
      character(len=200) :: string
      length = 200
      do while( empty_char(string(length:length)) ) 
        length = length - 1
      end do
    
    end function length
    
    !
    ! Function that determines if a character is empty
    !
    
    logical function empty_char(char)
    
      implicit none
      character :: char
      empty_char = .false.
      if ( char == achar(9) .or. &
           char == achar(32) .or. &
           char == '' ) then
        empty_char = .true.
      end if 
    
    end function empty_char
    
    !
    ! Function that checks if a line is a comment line
    !
    
    logical function comment(string)
      
      implicit none
      integer :: i
      character(len=200) :: string
      i = 1
      do while( empty_char(string(i:i)) .and. i < 200 ) 
        i = i + 1
      end do
      comment = .false.
      if ( string(i:i) == "#" .or. i == 200 ) comment = .true.
    
    end function comment

    !
    ! Subroutine that tests if file exists, and tries to open it.
    ! Asks the user if he/she wants the file to be overwritten 
    !

    subroutine checkfile(file)
    
      implicit none
      integer :: ioerr
      character(len=200) :: file
      character(len=1) :: char
    
      open(10,file=file,status='new',action='write',iostat=ioerr)
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Trying to create file: ', trim(adjustl(file)),' but file already exists '
        write(*,"(a,$)") '  Overwrite it? (Y/N): '
        read(*,*) char
        if ( char == "Y" ) then
          open(10,file=file,action='write',iostat=ioerr)
          if ( ioerr /= 0 ) then
            write(*,*) ' Could not open file. Quitting. '
            stop
          end if
          close(10)
        else
          write(*,*) ' Quitting. '
          stop
        end if
      end if
      close(10)
    
    end subroutine checkfile

end module file_operations









