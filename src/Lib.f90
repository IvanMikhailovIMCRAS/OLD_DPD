!######################################################################################
module Lib
!######################################################################################
	use CommonParam
    implicit none
		
Contains



!***************************************************************************************
subroutine periodic_test(x, Box) ! test of periodic boundary conditions
!***************************************************************************************
	real(4),intent(inout) :: x
	real(4),intent(in) :: Box

	if (abs(x).gt.0.5*Box) x = x - x/abs(x)*Box

return
end subroutine periodic_test
!***************************************************************************************

!**************************************************************************************
subroutine scan_list(N, num_cell, N_dblist, N_dbcell, main_list, att_list, list1,list2, &
&                    bond_list, ub1, ub2)
!**************************************************************************************
	integer(4),intent(in):: N, num_cell, N_dblist, N_dbcell
	integer(4), intent(in) :: main_list(1:num_cell), att_list(1:N)
	integer(4), intent(in) :: list1(0:N_dbcell), list2(1:N_dbcell), bond_list(1:N)
	integer(4),intent(out) :: ub1(0:N_dblist),ub2(1:N_dblist)
	integer(4) i, j, num_ub, n1, n2
	real(4) dx, dy, dz, r2
	
	num_ub = 0
	
	! pairs of interacting beads in the same cell
	do i = 1, num_cell
		n1 = main_list(i)
		do while (n1.ne.0)
			n2 = att_list(n1)
			do while (n2.ne.0)
				if (n1.lt.n2) then
					if (bond_list(n2).ne.n1) then
						num_ub = num_ub + 1
						ub1(num_ub) = n1; ub2(num_ub) = n2
					endif
				else
					if (bond_list(n1).ne.n2) then
						num_ub = num_ub + 1
						ub1(num_ub) = n2; ub2(num_ub) = n1
					endif
				endif
				n2 = att_list(n2)
			enddo
			n1 = att_list(n1)
		enddo
	enddo
		
	! pairs of interacting beads in neighboring cells
	do i = 1, list1(0)
		n1 = main_list(list1(i))
		do while (n1.ne.0)
			n2 = main_list(list2(i))
			do while (n2.ne.0)
				if (n1.lt.n2) then
					if (bond_list(n2).ne.n1) then
						num_ub = num_ub + 1
						ub1(num_ub) = n1; ub2(num_ub) = n2
					endif
				else
					if (bond_list(n1).ne.n2) then
						num_ub = num_ub + 1
						ub1(num_ub) = n2; ub2(num_ub) = n1
					endif
				endif
				n2 = att_list(n2)
			enddo
			n1 = att_list(n1)
		enddo
	enddo

	ub1(0) = num_ub
return
end subroutine scan_list
!**************************************************************************************


!**************************************************************************************
subroutine attached_list(num_cell,num_atom,x,y,z,nx,ny,nz,BOX,main_list,att_list)
!**************************************************************************************
	integer(4), intent(in) :: num_cell, num_atom, nx, ny, nz
	real(4), intent(in) :: x(1:num_atom), y(1:num_atom), z(1:num_atom), BOX(3) 
	integer(4), intent(inout) :: main_list(1:num_cell), att_list(1:num_atom)
	integer(4) i, n, na, kx, ky, kz
	
	main_list(:) = 0; att_list(:) = 0
	
	do i = 1, num_atom
		! calculate the "coordinate" of a cell
		! Warning if the bead is on the border of the cell,
		! the "coordinate" of the cell cannot be calculated correctly,
		! that's why the magic "0.4999998" is used
		kx = int((x(i)/BOX(1) + 0.499998)*nx) + 1
		ky = int((y(i)/BOX(2) + 0.499998)*ny) + 1
		kz = int((z(i)/BOX(3) + 0.499998)*nz) + 1
		if (kz.lt.1) kz = 1
		if (kz.gt.nz) kz = nz
		n = ((kx-1)*ny + ky-1)*nz + kz  
		if (main_list(n).eq.0) then
			main_list(n) = i
		else
			na = main_list(n)
			do
				if (att_list(na).eq.0) then
					att_list(na) = i; exit
				else
					na = att_list(na)
				endif
			enddo
		endif		
	enddo

	return
end subroutine attached_list
!**************************************************************************************

!**************************************************************************************
subroutine DB_LIST(N_dbcell,BOX,list1,list2) ! Getting list of interacting beads
!**************************************************************************************
	integer(4), intent(in) :: N_dbcell
	real(4), intent(in) :: BOX(3)
	integer(4), intent(inout) :: list1(0:N_dbcell), list2(1:N_dbcell) 
	integer(4) i, j, k, nx, ny, nz, n, n_connect
	real(4) Lx, Ly, Lz, diagonal, dist, dx, dy, dz
	real(4), allocatable, dimension(:) :: x, y, z
	
	nx = int(BOX(1)); ny = int(BOX(2)); nz = int(BOX(3)) 
	n = nx * ny * nz
	Lx = BOX(1)/nx;   Ly = BOX(2)/ny;   Lz = BOX(3)/nz
	
	diagonal = Lx**2 + Ly**2 + Lz**2
	diagonal = diagonal * 1.1
	
	allocate(x(1:N)); allocate(y(1:N)); allocate(z(1:N))
	x(:) = 0.0;       y(:) = 0.0;       z(:) = 0.0
	
	n = 0
	do i = 1, nx
		do j = 1, ny
			do k = 1, nz
				n = n + 1
				x(n) = Lx * i - 0.5 * Lx - 0.5 * BOX(1)
				y(n) = Ly * j - 0.5 * Ly - 0.5 * BOX(2)
				z(n) = Lz * k - 0.5 * Lz - 0.5 * BOX(3)
			enddo
		enddo
	enddo
	
	n_connect = 0
	do i = 1, n-1
		do j = i+1, n
			dx = x(j) - x(i); dy = y(j) - y(i); dz = z(j) - z(i)
			call periodic_test(dx, BOX(1))
			call periodic_test(dy, BOX(2))
			if (.not. wall_on) call periodic_test(dz, BOX(3))
			dist = dx**2 + dy**2 + dz**2
			if (dist.lt.diagonal) then
				n_connect = n_connect + 1
				list1(n_connect) = i; list2(n_connect) = j; 
			endif
		enddo
	enddo
	
	list1(0) = n_connect

return
end subroutine DB_LIST
!**************************************************************************************

!**************************************************************************************
subroutine angular_projections(xa,  xb,  xc,  ya,  yb,  yc,  za,  zb,  zc, &
&                              fxa, fxb, fxc, fya, fyb, fyc, fza, fzb, fzc) 
!**************************************************************************************
	real(4), intent(in) ::   xa,  xb,  xc,  ya,  yb,  yc,  za,  zb,  zc
	real(4), intent(out) :: fxa, fxb, fxc, fya, fyb, fyc, fza, fzb, fzc
	real(4) x_ba, x_cb, y_ba, y_cb, z_ba, z_cb, x_abc, y_abc, z_abc
	real(4) x1, x2, y1, y2, z1, z2, r1, r2, A
	
	x_ba = xb - xa ; x_cb = xc - xb ; x_abc = xa - 2.0*xb + xc
	y_ba = yb - ya ; y_cb = yc - yb ; y_abc = ya - 2.0*yb + yc 
	z_ba = zb - za ; z_cb = zc - zb ; z_abc = za - 2.0*zb + zc
	
	x1 = x_ba**2; y1 = y_ba**2; z1 = z_ba**2
	x2 = x_cb**2; y2 = y_cb**2; z2 = z_cb**2
	
	r1 = sqrt(x1 + y1 + z1); r2 = sqrt(x2 + y2 + z2); 
	A = x_ba*x_cb + y_ba*y_cb + z_ba*z_cb
	
	fxa = (-x_cb + x_ba*A/r1**2)/r1/r2
	fxc = ( x_ba - x_cb*A/r2**2)/r1/r2
	fxb = (x_abc - (x_ba/r1**2 - x_cb/r2**2)*A)/r1/r2   ! ???
	
	fya = (-y_cb + y_ba*A/r1**2)/r1/r2
	fyc = ( y_ba - y_cb*A/r2**2)/r1/r2
	fyb = (y_abc - (y_ba/r1**2 - y_cb/r2**2)*A)/r1/r2   ! ???
	
	fza = (-z_cb + z_ba*A/r1**2)/r1/r2
	fzc = ( z_ba - z_cb*A/r2**2)/r1/r2
	fzb = (z_abc - (z_ba/r1**2 - z_cb/r2**2)*A)/r1/r2  ! ???

return
end subroutine angular_projections
!**************************************************************************************

    
!######################################################################################
end module Lib
!######################################################################################
