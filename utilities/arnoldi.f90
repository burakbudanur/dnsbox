!*********************************************************************
!  n        dimension of system, leading dimension of sv,q,b
!  k          set k==0 for initialisation, 
!         on exit, number of completed calls/iterations.
!  kmax        max number of iterations, leading dim of h
!  ncgd        number of converged eigenvalues requested
!  inprod      dot product: d = inprod(0,a,b)
!  tol        tolerence of converged eigvals from H_{k-1} and H_k
!  sv        on entry, A.q_k ,
!         on exit, latest basis vector q_{k+1}.
!  h        Hessenburg matrix H
!  q        Krylov subspace Q
!  b        converged eigenvectors
!  wr2, wi2    converged eigenvalues.
!  ifail == -1  converged!!!
!            0  ok, continue
!            1  reached max iterations
!            2  invariant subspace
!          > 2  error from Lapack
!                         A. P. Willis 2002
!*********************************************************************

      subroutine arnold(n,k,kmax,ncgd,inprod,tol,sv,h,q,b,wr2,wi2,ifail)

! arguments
      integer(i4) ::       n,k,kmax,ncgd
      real(dp) :: inprod
      external          inprod
      real(dp) :: tol, sv(n), q(n,kmax), b(n,ncgd*2),  &
                  wr2(ncgd), wi2(ncgd), h(kmax,kmax)
      integer(i4) ::       ifail

! local variables
      integer(i4) ::       i,j,m,p
      real(dp) ::   norm,nm

! lapack stuff
      real(dp) ::   scale(kmax), tau(kmax),  &
                 wr(kmax), wi(kmax), work((kmax+2)*kmax),  &
                 h2(kmax,kmax), h3(kmax,kmax), vr(kmax,ncgd*2)
      logical ::       select(kmax)
      integer(i4) ::       ilo, ihi, ifailr(ncgd*2)


!---------------------------------------------------------------------
      ifail = 0

!     too many iterations
      if(k>=kmax) then
         ifail = 1
         return
      end if

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Arnoldi Process
!
!     copy soln vect to krylov subspace at k+1
      do i = 1, n
         q(i,k+1) = sv(i)
      end do

!     modified gram-schmidt
      do j = 1, k
!        get inner product --> k^th col of kxk matrix h
         h(j,k) = inprod(0,q(1,k+1),q(1,j))
!        orthogonalise against previous k vectors
         do i = 1, n
            q(i,k+1) = q(i,k+1) - h(j,k)*q(i,j)
         end do
      end do
!     reorthogonolise
      do j = 1, k
         norm = inprod(0,q(1,k+1),q(1,j))
         do i = 1, n
            q(i,k+1) = q(i,k+1) - norm*q(i,j)
         end do
         h(j,k) = h(j,k) + norm
      end do

!     normalise the k+1 basis vector
      norm = inprod(0,q(1,k+1),q(1,k+1))
      norm = sqrt(norm)
      if(norm<small) then
         ifail = 2
         return
      end if
      do i = 1, n
         q(i,k+1) = q(i,k+1) / norm
         sv(i) = q(i,k+1)
      end do

!     for starting vector just copy to subspace and normalise,
!     zero the hessenburg matrix
      if(k==0) then
         do j = 1, kmax
            do i = 1, kmax
               h(i,j) = 0
            end do
         end do
         goto 199
      end if

!     subdiagonal of hessenburg
      h(k+1,k) = norm
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Eigenvalues / Eigenvectors
!
!     need large enough matrix to get eigvecs
      if(k<ncgd) goto 199

! get eigvals and eigvecs of kxk hessenberg matrix h by QR
      do j = 1, k
         do i = 1, kmax
            h2(i,j) = h(i,j)
         end do
      end do

!     balance the matrix
      call dgebal('B',k,h2,kmax,ilo,ihi,scale,ifail)      
      if(ifail/=0) then
         ifail = ifail + 1000
         return
      end if

!     reduce
      m = (kmax+2)*kmax
      call dgehrd(k,ilo,ihi,h2,kmax,tau,work,m,ifail)
      if(ifail/=0) then
         ifail = ifail + 2000
         return
      end if

!     QR
      do j = 1, k
         do i = 1, kmax
            h3(i,j) = h2(i,j)
         end do
      end do
      call dhseqr('E','N',k,ilo,ihi,h3,kmax,wr,wi,0.0_dp,1,work,m,ifail)      
      if(ifail/=0) then
         ifail = ifail + 3000
         return
      end if
      
!     order the eigvals, largest magnitude
!     compare with eigvals from k-1
!     copy new values to old
      do i = 1, k
         select(i) = .false.
      end do
      p = 0
      do i = 1, ncgd
         norm = -1.0_dp
         do j = 1, k
            if(.not.select(j)) then
               nm = wr(j)*wr(j) + wi(j)*wi(j)
               if(nm>norm) then
                  norm = nm
                  m = j
               end if            
            end if
         end do
         select(m) = .true.
         if(abs(wi(m))>small) select(m+1)=.true.
         nm = sqrt((wr(m)-wr2(i))**2 + (wi(m)-wi2(i))**2)
         if(nm < tol) p = p+1
         wr2(i) = wr(m)
         wi2(i) = wi(m)
      end do

!     check if reqiured num converged...
      if(p<ncgd) goto 199
      
!     Enough converged!! 
!     Get eigvecs of hessenburg by inverse iteration.
      do i = 1, k
         if(select(i) .and. abs(wi(i))>small) select(i+1)=.false.
      end do
      call dhsein('R','Q','N',select,k,h2,kmax,wr,wi,  &
                  0.0_dp,1,vr,kmax,ncgd*2,m, work, 1,ifailr,ifail )
      if(ifail/=0) then
         ifail = ifail + 4000
         return
      end if
      
      call dgebak('B','R',k,ilo,ihi,scale,m,vr,kmax,ifail)
      if(ifail/=0) then
         ifail = ifail + 5000
         return
      end if

!     ? check correct number of cols
      p = 0
      do i = 1, k
         if(select(i)) then
            p = p+1
            if(abs(wi(i))>small) p = p+1
         endif
      end do
      if(p/=m) then
         ifail = m + p*100 + 6000
      end if

!     order the cols of vr
      m = 0
      do j = 1, k
         if(select(j)) then
            i = 0
            p = 0
 101        continue
            i = i+1
            if(abs(wr2(i)-wr(j))>small .or. abs(wi2(i)-wi(j))>small) then
               p = p+1
               if(abs(wi2(i))>small) p = p+1
               goto 101
            end if
            p = p+1
            m = m+1
            do i = 1, k
               work((p-1)*k+i) = vr(i,m)
            end do
            if(abs(wi(j))>small) then
               p = p+1
               m = m+1
               do i = 1, k
                  work((p-1)*k+i) = vr(i,m)
               end do            
            end if
         end if
      end do
      do p = 1, m
         do i = 1, k
            vr(i,p) = work((p-1)*k+i)
         end do
      end do

!     Get eigvecs, m cols
      do p = 1, m
         do i = 1, n
            b(i,p) = 0
            do j = 1, k
               b(i,p) = b(i,p) + q(i,j)*vr(j,p)
            end do
         end do
      end do

!     converged !!!
      ifail = -1
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 199  continue
      k = k+1

      end subroutine

!*********************************************************************
!  wr + i wi = e^{sigma T}
!  overwrite wr,wi with sigma.
!  do nothing to eigvecs.
!  if  w v = A v  with w,v complex, then equating re and imag parts
!     wr vr - wi vi = A vr  and  wi vr + wr vi = A vi
!
!*********************************************************************
      subroutine logeig(ncgd,wr,wi,t,ifail)

! arguments
      integer(i4) ::       ncgd
      real(dp) :: wr(ncgd), wi(ncgd), t
      integer(i4) ::       ifail

! local variables
      integer(i4) ::       i
      real(dp) ::   norm,re,im,th
      
!---------------------------------------------------------------------
      ifail = 0

!     log of the eigvals
      do i = 1, ncgd
         re = wr(i)
         im = wi(i)
         if( abs(re)<small ) then
            if( abs(im)<small ) then
               ifail = i
               return
            else if( im>0 ) then
               th =  0.5_dp*pi
            else 
               th = -0.5_dp*pi
            end if
         else if( re>0 ) then
            th = atan(im/re)
         else
            if( im>0 ) then
               th =  pi + atan(im/re)
            else
               th = -pi + atan(im/re)
            end if
         end if
         norm = sqrt( re*re + im*im )
         wr(i) = log( norm ) / t
         wi(i) = th / t
      end do

      end subroutine
!*********************************************************************
