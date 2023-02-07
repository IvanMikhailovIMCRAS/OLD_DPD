!######################################################################################
program DPD ! The program for Dissipative Particle Dynamic simulation 
!######################################################################################
	use CommonParam
	use ErrorList
	implicit none
	integer(4) N, NB, num_types ! number of beads, number of bonds, number of the bead types
	real(4) BOX(3) ! box size
	real(4), allocatable, dimension(:)    :: x, y, z, vx, vy, vz
	real(4), allocatable, dimension(:,:)  :: FF          ! "force-field"
	integer(4), allocatable, dimension(:) :: b1, b2      ! bond lists
	integer(4), allocatable, dimension(:) :: typ         ! list of the bead types
	integer(4), allocatable, dimension(:) :: list_fix    ! list of fixed beads
	integer(4), allocatable, dimension(:,:) :: list_ang  ! list of "angles"
	integer(4) i, ioer, tmp_int, num_fix, num_ang
	character(1) tmp_char
	character(255) :: tmp_string
	character(2) :: flag
	
	path = ' '
	call get_command_argument(1, flag)
	if (flag == '-p') call get_command_argument(2, path)
	
	! create INFOR-file to report about the code execution and the fatal errors	
	open(n_infor,file=(trim(path) // name_infor))
	! TRACK-file must be not found in current folder to save old TRACK
	open(n_track,file=(trim(path) // name_track),iostat=ioer,status='new')
	if (ioer.ne.0) call ERROR(1)
	close(n_track)
	! open COORD-file with coordinates (it must be exist in folder)
	open(n_coord,file=(trim(path) // name_coord),iostat=ioer,status='old')
	if (ioer.ne.0) call ERROR(2)
	! read the numbers of beads and the box size
	read(n_coord,'(a)',advance='NO',iostat=ioer) tmp_string
	read(tmp_string,*,iostat=ioer) tmp_char, N, tmp_char, BOX(1), BOX(2), BOX(3)
	if (ioer.ne.0) then
		read(tmp_string,*,iostat=ioer) tmp_char, N, tmp_char, BOX(1)
		BOX(2) = BOX(1); BOX(3) = BOX(1)              
	endif
	if (ioer.ne.0.or.N.lt.3.or.BOX(1).le.0.0) call ERROR(3)
	write(n_infor,'(a,I0,a,F16.8)') 'N = ', N, '  BOX = ', BOX(1)
	! read coordinates of beads
	allocate(x(1:N)); allocate(y(1:N)); allocate(z(1:N)) ! allocate memory for x,y,z
	allocate(typ(1:N))
	do i = 1, N
		read(n_coord,*,iostat = ioer) tmp_int, x(i), y(i), z(i), typ(i)
		if (ioer.ne.0 &
		& .or.abs(x(i)).gt.0.5*BOX(1).or.abs(y(i)).gt.0.5*BOX(2)) then  
			write(n_infor,'(a,I0)') 'COORD file: string: ', i+1
			call ERROR(4)
		endif
	enddo
	close(n_coord)
	write(n_infor,'(a)') 'Coordinates have been read successfully.' 
	! open VELOC-file with velocities (if it absents, put all velosities being equal 0)
	open(n_veloc,file=(trim(path) // name_veloc),iostat=ioer,status='old')
	allocate(vx(1:N)); allocate(vy(1:N)); allocate(vz(1:N)) ! allocate memory for vx,xy,vz
	if (ioer.ne.0) then
		vx(:) = 0.0; vy(:) = 0.0; vz(:) = 0.0
		write(n_infor,'(a)') 'Warning! File VELOC does not exist! Velocities are set equal to zero.'
		write(n_infor,'(a)') 'Velocities are set equal to zero.'
	else
		read(n_veloc,*,iostat=ioer) tmp_char
		if (ioer.ne.0) call ERROR(5)
		do i = 1, N
			read(n_veloc,*,iostat = ioer) tmp_int, vx(i), vy(i), vz(i)
			if (ioer.ne.0) then 
				write(n_infor,'(a,I0)') 'VELOC file: string: ', i+1
				call ERROR(4)
			endif
		enddo
		close(n_veloc)
		write(n_infor,'(a)') 'Velocities have been read successfully.' 
	endif
	
	! open BONDS-file containing a list of bonds (it must be exist in folder)
	open(n_bonds,file=(trim(path) // name_bonds),iostat=ioer,status='old')
	if (ioer.ne.0) call ERROR(6)
	! read number of bonds
	read(n_bonds,*,iostat=ioer) tmp_char, NB
	if (ioer.ne.0.or.NB.lt.0) call ERROR(7)
	write(n_infor,'(a,I0)') 'NB = ', NB
	allocate(b1(1:NB)); allocate(b2(1:NB)) ! allocate memory for the bond lists
	do i = 1, NB
		read(n_bonds,*,iostat=ioer) b1(i), b2(i)
		if (ioer.ne.0.or.b1(i).ge.b2(i) &
		& .or.b1(i).lt.1.or.b1(i).gt.N.or.b2(i).lt.1.or.b2(i).gt.N) then 
			write(n_infor,'(a,I0)') 'BONDS file: string: ', i+1
			call ERROR(4)
		endif
	enddo
	close(n_bonds)
	write(n_infor,'(a)') 'Bonds have been read successfully.' 
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(n_field,file=(trim(path) // name_field),iostat=ioer,status='old')
	if (ioer.ne.0) call ERROR(10)
	! read the number of the bead types
	read(n_field,*,iostat=ioer) tmp_char, num_types
	if (ioer.ne.0.or.num_types.lt.1) call ERROR(11)
	write(n_infor,'(a,I0)') 'num_types = ', num_types
	allocate(FF(1:num_types,1:num_types)) ! allocate memory for "force-field" matrix
	read(n_field,*,iostat=ioer) FF
	if (ioer.ne.0) call ERROR(11)
	close(n_field)
	write(n_infor,'(a)') 'Force-field have been read successfully.' 
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(n_fixed,file=(trim(path) // name_fixed),iostat=ioer,status='old')
	if (ioer.ne.0) then
		num_fix = 0
	else
		! read list of the fixed beads
		read(n_fixed,*,iostat=ioer) tmp_char, num_fix
		if (ioer.ne.0.or.num_fix.lt.1) call ERROR(12)
		write(n_infor,'(a,I0)') 'num_fixed = ', num_fix
		allocate(list_fix(1:N)); list_fix(:) = 0
		do i = 1, num_fix
			read(n_fixed,*,iostat=ioer) list_fix(i)
			if (ioer.ne.0.or.list_fix(i).gt.N.or.list_fix(i).lt.1) call ERROR(12)
		enddo
	endif
	
	close(n_fixed)
	write(n_infor,'(a)') 'Fixed points have been read successfully.' 
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(n_angls,file=(trim(path) // name_angls),iostat=ioer,status='old')
	if (ioer.ne.0) then
		num_ang = 0
	else
		! reading numbers of beads representing valence angles
		read(n_angls,*,iostat=ioer) tmp_char, num_ang
		if (ioer.ne.0.or.num_ang.lt.1) call ERROR(12)
		write(n_infor,'(a,I0)') 'num_ang = ', num_ang
		allocate(list_ang(1:3,1:N)); list_ang(:,:) = 0
		do i = 1, num_ang
			read(n_angls,*,iostat=ioer) list_ang(:,i)
			if (ioer.ne.0) call ERROR(12)
		enddo
	endif
	close(n_angls)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	write(n_infor,'(a)') 'Fixed points have been read successfully.' 
	
	
	call control(N, NB, BOX, x, y, z, vx, vy, vz, typ, b1, b2, num_types, FF, &
&   			num_fix, list_fix, num_ang, list_ang)
	
	deallocate(x); deallocate(y); deallocate(z)
	deallocate(vx); deallocate(vy); deallocate(vz)
	deallocate(typ); deallocate(b1); deallocate(b2)
	deallocate(FF)
	if (num_fix.gt.0) deallocate(list_fix)
	write(n_infor,'(a)') 'Awesome! That worked out!' 
	close(n_infor)
	
	stop
end program DPD
!######################################################################################

!**************************************************************************************
subroutine control(N, NB, BOX, x, y, z, vx, vy, vz, typ, b1, b2, num_types, FF, &
&					num_fix, list_fix, num_ang, list_ang)
!**************************************************************************************
	use CommonParam
	use ErrorList
	use Engine
	use Lib
	implicit none
	integer(4), intent(in) :: N, NB, num_types, num_fix, num_ang
	real(4), intent(in) :: BOX(3)
	real(4), intent(inout) :: x(1:N), y(1:N), z(1:N), vx(1:N), vy(1:N), vz(1:N)
	real(4), intent(inout) :: FF(1:num_types,1:num_types)
	integer(4), intent(inout) :: typ(1:N), list_fix(1:N), list_ang(1:3,1:num_ang)
	integer(4), intent(in) :: b1(1:NB), b2(1:NB)
	integer(4) num_step, num_snapshot, num_track, num_statis
	real(4) dt, l_bond, k_bond, k_ang, wall
	integer(4) i, ioer, i_tmp
	integer(4) bond_list(1:N)                                             
	integer(4) nx, ny, nz, num_cell, N_dblist, N_dbcell 
	integer(4), allocatable, dimension(:) :: list1, list2, main_list, att_list
	
	open(n_contr,file=(trim(path) // name_contr),iostat=ioer,status='old')
	if (ioer.ne.0) call ERROR(8)
	
	read(n_contr,*,iostat=ioer) dt           ! size of the time step
	if (ioer.ne.0.or.dt.le.0.0) call ERROR(9)
	read(n_contr,*,iostat=ioer) l_bond       ! the bond length
	if (ioer.ne.0.or.l_bond.le.0.0) call ERROR(9)
	read(n_contr,*,iostat=ioer) k_bond       ! coefficient of the bond rigidity
	if (ioer.ne.0.or.k_bond.lt.0.0) call ERROR(9)
	read(n_contr,*,iostat=ioer) k_ang        ! coefficient of the angle rigidity
	if (ioer.ne.0.or.k_ang.lt.0.0) call ERROR(9)
	read(n_contr,*,iostat=ioer) num_step     ! number of steps
	if (ioer.ne.0.or.num_step.lt.0) call ERROR(9)
	read(n_contr,*,iostat=ioer) num_snapshot ! how often make snapshots of system
	if (ioer.ne.0.or.num_snapshot.lt.1) call ERROR(9)
	read(n_contr,*,iostat=ioer) num_track    ! how often print the coordinates
	if (ioer.ne.0.or.num_track.lt.1) call ERROR(9)
	read(n_contr,*,iostat=ioer) wall         ! interaction parameter between wall and beads
	if (ioer.ne.0.or.wall.lt.0.0) call ERROR(9)
	read(n_contr,*,iostat=ioer) i_tmp        ! exist wall normal to z-direction (0 or 1)
	if (ioer.ne.0) call ERROR(9)
	wall_on = (i_tmp.eq.1)                   ! logical label about existence of wall
	read(n_contr,*,iostat=ioer) solvent, i_tmp ! type of the solvent beads and (0, 1) - to show the solvent
	if (ioer.ne.0) call ERROR(9)
	solvent_on = (i_tmp.eq.1) 
	read(n_contr,*,iostat=ioer) i_tmp
	if (ioer.ne.0) call ERROR(9)
	interact_on = (i_tmp.eq.1)
	read(n_contr,*,iostat=ioer) UL, UH, dH   ! attractive potential
	if (ioer.ne.0.or.dH.le.0.0) call ERROR(9)
		    	
	close(n_contr)
	
	num_solvent = 0
	do i = 1, N
		if (typ(i).eq.solvent) num_solvent = num_solvent + 1
	enddo
	
	write(n_infor,'(a)') 'CONTR-file have been read successfully.' 
	
	bond_list(:) = 0
	do i = 1, NB
		bond_list(b2(i)) = b1(i)     
	enddo

	nx = int(BOX(1)); ny = int(BOX(2)); nz = int(BOX(3))
	num_cell = nx * ny * nz   ! number of cells
	N_dbcell = 13 * num_cell  ! number of cell pairs
	N_dblist = N * 256        ! paars numbers of the interacting beads
	
	allocate(list1(0:N_dbcell));   allocate(list2(1:N_dbcell)) 
    allocate(main_list(num_cell)); allocate(att_list(N))
    list1(:) = 0; list2(:) = 0; main_list(:) = 0; att_list(:) = 0
    
    ! one time calculate list of the interacting cells
    call DB_LIST(N_dbcell,BOX,list1,list2)
    
    call main(N, NB, BOX, x, y, z, vx, vy, vz, typ, b1, b2, bond_list, dt, FF, &
	&        l_bond, k_bond, num_step, num_snapshot, num_track, N_dblist, num_fix, list_fix, &
	&        num_types, num_cell, N_dbcell, main_list, att_list, list1, list2, nx,ny,nz,wall, &
	&        num_ang, list_ang, k_ang)

	
	deallocate(list1); deallocate(list2); deallocate(main_list); deallocate(att_list)
	
	
	
	return
end subroutine control
!**************************************************************************************
