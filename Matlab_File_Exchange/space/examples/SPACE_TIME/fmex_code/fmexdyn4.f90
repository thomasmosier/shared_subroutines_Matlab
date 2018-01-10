!C==============================================================
!C	Subroutine tnn(xcoord,ycoord,m,n,istart,nnmat)
!C	real(kind=8), intent(in):: xcoord(:),ycoord(:)
!C	real(kind=8) xyc2(size(xcoord,1)),d2(size(xcoord,1))
!C	integer, intent(out)::nnmat(:,:)
!C	integer, intent(in)::n,istart
!C	integer, intent(inout)::m
!C	integer irngt(size(xcoord,1))
!C	integer i,j


!C     The gateway subroutine
		subroutine mexfunction(nlhs, plhs, nrhs, prhs)

		integer plhs(*), prhs(*)

		integer pr1_in 
   		integer pr2_in 
		integer pr3_in
		integer pr4_in
		integer pr5_in

		integer m1_in, n1_in, size1
		integer m2_in, n2_in, size2
		integer m3_in, n3_in, size3
		integer m4_in, n4_in, size4
		integer m5_in, n5_in, size5

		integer pr_out, m_out, n_out, istart
    
		integer mxGetPr, mxcreatedoublematrix
		integer nlhs, nrhs, mxGetM, mxGetN

		real(kind=8) mreal, istartreal

		real(kind=8), allocatable::xcoord(:), ycoord(:)
		integer, allocatable:: nnmat(:,:)
	
!xcoord	        
      m1_in = mxGetM(prhs(1))
      n1_in = mxGetN(prhs(1))
      size1 = m1_in * n1_in
      pr1_in = mxGetPr(prhs(1))

!ycoord	
      m2_in = mxGetM(prhs(2))
      n2_in = mxGetN(prhs(2))
      size2 = m2_in * n2_in
      pr2_in = mxGetPr(prhs(2))

!m	
      m3_in = mxGetM(prhs(3))
      n3_in = mxGetN(prhs(3))
      size3 = m3_in * n3_in
      pr3_in = mxGetPr(prhs(3))
	
!n	
      m4_in = mxGetM(prhs(4))
      n4_in = mxGetN(prhs(4))
      size4 = m4_in * n4_in
      pr4_in = mxGetPr(prhs(4))

!istart	
      m5_in = mxGetM(prhs(5))
      n5_in = mxGetN(prhs(5))
      size5 = m5_in * n5_in
      pr5_in = mxGetPr(prhs(5))

	call mxCopyPtrToReal8(pr3_in, mreal, size3)
	m_out=int(mreal)

!	OPEN(UNIT=666,FILE='m_out.ASC')
!	WRITE(666,*) mreal, m_out, m3_in, n3_in, size3
!	CLOSE(666)

	n_out=m1_in

	call mxCopyPtrToReal8(pr5_in, istartreal, size5)
	istart=int(istartreal)

!	OPEN(UNIT=666,FILE='istart_out.ASC')
!  	WRITE(666,*) istartreal,istart,m5_in,n5_in,size5
!	CLOSE(666)
	
	  allocate(xcoord(m1_in), ycoord(m2_in), nnmat(m_out, n_out))

	  call mxCopyPtrToReal8(pr1_in,xcoord,n_out)
	  call mxCopyPtrToReal8(pr2_in,ycoord,n_out)

  

!	Call the computational routine.


!	Subroutine tnn(xcoord,ycoord,mvec,n,istart,nnmat)
      
	call tnmex4(xcoord, ycoord, m_out, n_out, istart, nnmat)

	plhs(1) = mxCreateDoubleMatrix(m_out, n_out, 0)
      	pr_out = mxGetPr(plhs(1))
	call mxcopyreal8toptr(dble(nnmat), pr_out, n_out*m_out)

!      	return
	end subroutine mexfunction

!	  tnmex3(xcoord,ycoord,m,n,istart,nnmat)

