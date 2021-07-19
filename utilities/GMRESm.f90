!----------------------------------------------------------------------
! Openpipeflow.org.  If used in your work, please cite
! Willis, A. (2017) SoftwareX 6, 124-127.
! https://doi.org/10.1016/j.softx.2017.05.003 (open access)
!                                      Thanks in advance! Ashley 2019.
!----------------------------------------------------------------------
! solve A x = b for x ;  
! minimise |Ax-b| subject to constraint |x| < delta .
! requires lapack routines dgelsy, dgesvd.
!----------------------------------------------------------------------
! m	  gmres dimension
! n 	  dimension of x
! x	  on input:  guess for x, can be 0
!         on exit:  solution x, subject to constraint if del>0
! b	  input b
! matvec  performs y := A x, call matvec(x, y)
! dotprd  dot product, d = dotprd(a,b)
! h       Hessian matrix,  size (m+1)*m
! v       Krylov subspace, size n*(m+1)
! res	  on input: |Ax-b|/|b|<res; 
!         on exit:  residual reached
! del     on input: if(del>0) then the x returned is the hookstep
!         on exit:  norm of next b predicted by hook
! its	  on input: max num its; 
!         on exit:  number of its taken
! info	  on input: if(info==1) print* residuals
!                   if(info==2) recalc hookstep with new del
! 	  on exit:  0 sucessful, 1 method breakdown, 2 max its
!							A.P.Willis 2008
!----------------------------------------------------------------------

 subroutine gmresm(m,n,x,b,matvec,dotprd,h,v,res,del,its,info,out)
   
   integer(i4),          intent(in)    :: m
   integer(i4),          intent(in)    :: n
   real(dp), intent(inout) :: x(n)
   real(dp), intent(in)    :: b(n)
   external :: matvec
   real(dp), external      :: dotprd
   real(dp), intent(inout) :: h(m+1,m)
   real(dp), intent(inout) :: v(n,m+1)
   real(dp), intent(inout) :: res
   real(dp), intent(inout) :: del
   integer(i4),          intent(inout) :: its
   integer(i4),          intent(inout) :: info
   integer(i4), intent(in) :: out
   real(dp) :: tol,res_,stgn, w(n), z(n)
   real(dp) :: h_(m+1,m), y(m+1), p(m+1), work(4*m+1)
   integer(i4) :: imx, piv(m), rank, i
   real(dp), save :: beta
   integer(i4), save :: j
   logical :: done   

   if(info==2) then
      call hookstep(j,h,m,beta,del,y,out)
      z = matmul(v(:,1:j),y(1:j))
      x = z
      info = 0
      return
   end if     

   tol = res
   imx = its
   its = 0
   v   = 0

 1 continue
   res_ = 1.0e99_dp
   stgn = 1.0_dp - 1.0e-14_dp
 
   beta = sqrt(dotprd(x,x)) 
   if(abs(beta)<small)  w = 0
   if(abs(beta)>=small) call matvec(x, w)
   w = b - w
   beta = sqrt(dotprd(w,w)) 
   v(:,1) = w / beta
     
   h = 0
   do j = 1, m
      its = its + 1
      z = v(:,j)      
      call matvec(z, w)
      do i = 1, j
        ! Fixed a possible bug here, size of v(1,i) and w do not match
        ! And the Gram-Schmidt doesn't make sense otherwise
         h(i,j) = dotprd(w,v(:,i))
         w = w - h(i,j)*v(:,i)
      end do
      h(j+1,j) = sqrt(dotprd(w,w))
      v(:,j+1) = w / h(j+1,j)
          
      p(1) = beta
      p(2:j+1) = 0
      h_(1:j+1,1:j) = h(1:j+1,1:j)
      call dgelsy(j+1,j,1,h_,m+1,p,m+1,piv,m,rank,work,4*m+1,i)
      if(i/=0) then
        write(out, *) 'gmresm: dgelsy'
        flush(out)
        stop
      end if
      y = p

      p(1:j+1) = - matmul(h(1:j+1,1:j),y(1:j))
      p(1) = p(1) + beta
      res = sqrt(dot_product(p(1:j+1),p(1:j+1)))
      if(info==1) write(out, *)  'gmresm: it=', its,' res=', res
      
      done = (res<=tol .or. its==imx .or. res>res_)
      if(done .or. j==m) then
         if(del>0.0_dp)  call hookstep(j,h,m,beta,del,y,out)
         z = matmul(v(:,1:j),y(1:j))
         x = x + z
         if(its==imx) info = 2
         if(res>res_) info = 1
         if(res<=tol) info = 0
         if(done)     return
         if(del>0.0_dp)  write(out, *)  'gmres: warning! restart affects hookstep'
         goto 1       ! (j==m) restart
      end if
      res_ = res*stgn
      
   end do   
 
 end subroutine gmresm
 
!-----------------------------------------------------------------
! replace y with a vector that generates a hookstep
! c.f. Viswanath (2008) arXiv:0809.1498
!-----------------------------------------------------------------
 subroutine hookstep(j,h,m,beta,del,y,out)
   
   integer(i4),          intent(in)    :: j, m
   real(dp), intent(in)    :: h(m+1,j), beta
   real(dp), intent(inout) :: del
   real(dp), intent(out)   :: y(j)
   integer(i4), intent(in) :: out
   real(dp) :: a(j+1,j), s(j), u(j+1,j+1), vt(j,j), work(5*(j+1))
   real(dp) :: p(j+1), q(j), mu, qn
   integer(i4) :: info
   
   a = h(1:j+1,1:j)
   
   call dgesvd('A','A',j+1,j,a,j+1,s,u,j+1,vt,j,work,5*(j+1),info)
   if(info/=0) then
        write(out,*) 'hookstep: dgesvd'
        flush(out)
        stop
   end if
   
   p(1:j) = beta * u(1,1:j)   

   mu = max(s(j)*s(j)*1.0e-6_dp,1.0e-99_dp)
   qn = 1.0e99_dp
   do while(qn>del)
      mu = mu * 1.1_dp
      q = p(1:j)*s/(mu+s*s)
      qn = sqrt(dot_product(q,q))
   end do

   y = matmul(q,vt)

   p = - matmul(h(1:j+1,1:j),y(1:j))
   p(1) = p(1) + beta
   del = sqrt(dot_product(p,p))
 
 end subroutine hookstep
