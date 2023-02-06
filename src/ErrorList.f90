!######################################################################################
module ErrorList
!######################################################################################
	use CommonParam
	implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
Contains

!**************************************************************************************
subroutine ERROR(n)
!**************************************************************************************
    integer(4) n
	
	select case (n)
        case(1)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. TRACK already exists!'
        case(2)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. COORD is not found!'
        case(3)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. First line in COORD is not correct!'
        case(4)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. Input data is invalid or does not exist!'
        case(5)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. First line in VELOC is not correct!'
        case(6)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. BONDS is not found'
        case(7)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. First line in BONDS is not correct'
        case(8)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. CONTR is not found!'
        case(9)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. CONTR is not correct!'
        case(10)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. FIELD is not found!'
        case(11)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. FIELD is not correct!'
        case(12)
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. FIXED is not correct!'
        case default
            write(n_infor,'(a,I0,a)') 'ERROR ',n,'. Unknown error!'
    end select
    
    
    close(n_infor)
    
    stop
end subroutine ERROR
!**************************************************************************************


!######################################################################################
end module ErrorList
!######################################################################################
