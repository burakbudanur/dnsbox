!----------------------------------------------------------------------
! Openpipeflow.org.  If used in your work, please cite
! Willis, A. (2017) SoftwareX 6, 124-127.
! https://doi.org/10.1016/j.softx.2017.05.003 (open access)
!                                      Thanks in advance! Ashley 2019.
!----------------------------------------------------------------------
! Find root f(x) = 0 by Newton method with gmres-hookstep:
!    invert  df(x_n)/dx . s = f(x_n)  
! for s subject to constraint |s| < del  (size of 'trust region'), 
!    then     x_{n+1} = x_n - s
! To define subroutines f,df below, write as
!    df . s = f
! where f is a function to be minimised and operator df provides the
! action of the Jacobian, but may also contain any constraints on the 
! update s (e.g. no shifts in a homogenous direction), arranged such 
! that the rhs of the constraint (in f) is 0 when converged.
!----------------------------------------------------------------------
!
! Allocate(new_x(1:n),new_fx(1:n)), then put initial guess in new_x(1:n).  
! This variable is updated with the current best solution each iteration.
!
! f       evaluates y=f(x):  call f(x, y)
! df      evaluates y=Jacobian.x, see GMRESm.f90:  call df(x, y)
! sub     parameterless subroutine called at end of each iteration,
!         so that data could be saved etc.
! inprod  dot product:  d = inprod(a,b)
! m	  gmres dimension (also max num gmres its)
! n 	  dimension of x
! gtol    tolerence for gmres, typically 0.001
! tol     request |f(x)|<tol
! del     initial size of trust region for usefulness of df
!         set = 0d0 for no hookstep, <0d0 for del = |f(x0)|/10.
! mndl    min size of trust region
! mxdl    max size of trust region, if del<0 set no upper limit
! nits	  max num Newton its 
! info	  on input if =1 print* details
!         on exit if =0 sucessful
!                    =1 gmres failed
!                    =2 reached max iterations, nits
!                    =3 trust region got too small
!                    =4 unable to allocate memory
!							A.P.Willis 2008
!----------------------------------------------------------------------

 subroutine newtonhook(f, df, sub, inprod, m, n, gtol, &
                       tol, del, mndl, mxdl, nits, info,out)
   
   external                        :: f, df, sub
   real(dp), external      :: inprod
   integer(i4),          intent(in)    :: m, n
   real(dp), intent(in)    :: gtol
   real(dp), intent(inout) :: tol, del, mndl, mxdl
   integer(i4),          intent(in)    :: nits
   integer(i4),          intent(inout) :: info
   integer(i4), intent(in) :: out
   real(dp), allocatable :: v(:,:)
   real(dp) :: tol_,x_(n),fx_(n), tol__,x__(n),fx__(n),del__
   real(dp) :: mxdl_, ared, pred, snrm, s(n)
   real(dp) :: gres, gdel, h((m+1)*m)
   integer(i4) :: ginfo

   real(dp) :: normOld, normNew
   
   write(out, *) 'beginning newton'
   flush(out)
   allocate(v(n,m+1), stat=ginfo)
   if(ginfo/=0) then 
      if(info==1) write(out, *)  'newton: unable to allocate memory'
      info = 4
      return
   end if

   new_nits = 0
   new_gits = 0
   new_del  = del
   mxdl_    = mxdl
   ginfo    = info
   call f(new_x, new_fx)
   new_tol  = sqrt(inprod(new_fx,new_fx))
   if(del<0.0_dp)  new_del = new_tol / 10.0_dp
   if(del<0.0_dp)  mxdl_   = 1.0e99_dp
    if(info==1)  then
        write(out, *) 'newton: nits=',new_nits,' res=', new_tol
    end if
   call sub()
   x_   = new_x
   fx_  = new_fx
   tol_ = new_tol
   tol__ = 1.0e99_dp

   if(new_tol<tol) then
    if(info==1)  then
        write(out, *) 'newton: input already converged'
        flush(out)
    end if
      info = 0
      deallocate(v)
      return
   end if

   normOld = sqrt(inprod(new_x,new_x))
   del__ = new_del
    ! - - - - Start main loop - - - - -  -
    do while(.true.)

        if(new_del<mndl) then
            if(info==1)  then
                write(out, *)  'newton: trust region too small'
                flush(out)
            end if
            info = 3
            deallocate(v)
            return      
        end if
                     ! find hookstep s and update x
      s        = 0
      gres     = gtol * new_tol
      gdel     = new_del
      if(ginfo/=2) new_gits = m
      if(abs(del)<small) new_gits = 9999
      call gmresm(m,n,s,fx_,df,inprod,h,v,gres,gdel,new_gits,ginfo,out)
      ginfo = info
      new_x = x_ - s
                     ! calc new norm, compare with prediction
      call f(new_x, new_fx)
      new_tol = sqrt(inprod(new_fx,new_fx))
      snrm = sqrt(inprod(s,s))
      ared = tol_ - new_tol
      pred = tol_ - gdel
         
      if(info==1) then 
         write(out, *) 'newton: nits=',new_nits,' res=', new_tol
         write(out, *) 'newton: gits=',new_gits,' del=', new_del
         write(out, *) 'newton: |s|=', snrm,' pred=', pred
         write(out, *) 'newton: ared/pred=', ared/pred
      end if

      if(abs(del)<small) then
        if(info==1) write(out, *) 'newton: took full newton step'
      else if(new_tol>tol__) then
         if(info==1) write(out, *)  'newton: accepting previous step'
         new_x   = x__
         new_fx  = fx__
         new_tol = tol__
         new_del = del__
      else if(ared<0.0_dp) then
         if(info==1) write(out, *)  'newton: norm increased, try smaller step'
         new_del = snrm * 0.5_dp
         ginfo   = 2
      else if(ared/pred<0.75_dp) then
         if(info==1)  write(out, *)  'newton: step ok, trying smaller step'
         x__     = new_x
         fx__    = new_fx
         tol__   = new_tol
         if(ared/pred> 0.1_dp)  del__ = snrm
         if(ared/pred<=0.1_dp)  del__ = snrm*0.5_dp
         new_del = snrm * 0.7_dp
         ginfo   = 2
      else if(snrm<new_del*0.9_dp) then
         if(info==1) write(out, *)  'newton: step good, took full newton step'
         new_del = min(mxdl_,snrm*2.0_dp)
      else if(new_del<mxdl_*0.9_dp) then
         if(info==1) write(out, *)  'newton: step good, trying larger step'
         x__     = new_x
         fx__    = new_fx
         tol__   = new_tol
         del__   = new_del
         new_del = min(mxdl_,snrm*2.0_dp)
         ginfo   = 2
      end if      
                         ! check if need to try another s
      if(ginfo==2) cycle
                         ! end of iteration
      new_nits = new_nits + 1
      call sub()
      x_   = new_x
      fx_  = new_fx
      tol_ = new_tol
      tol__ = 1.0e99_dp

      ! Re-compute parameters after each call of sub()
      ! Otherwise when Newton stops is not consistent with the relative error
      ! requirement in general
      normNew = sqrt(inprod(new_x,new_x))
      if(abs(normNew)<small) normNew=1.0_dp
      tol = (tol / normOld) * normNew
      del = (del / normOld) * normNew
      mndl = (mndl / normOld) * normNew
      mxdl = (mxdl / normOld) * normNew
      normOld = normNew

      new_del = del
      mxdl_ = mxdl
      if(del<0.0_dp)  new_del = new_tol / 10.0_dp
      if(del<0.0_dp)  mxdl_   = 1.0e99_dp

      if(new_tol<tol) then
         if(info==1) write(out, *)  'newton: converged'
         info = 0
         deallocate(v)
         return
      else if(new_nits==nits) then
         if(info==1) write(out, *)  'newton: reached max its'
         info = 2
         deallocate(v)
         return
      end if
      
   end do

 end subroutine newtonhook

 