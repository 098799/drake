module ggint
  use precision, only : prec
  use file_OUT, only : LOUT
  implicit none
  integer, parameter :: newhn=120
  integer, parameter :: ggn=400,hn=200,ghn=100,tho=1000,RECB=32,RECS=16!RECB=40,RECS=24
  integer, parameter :: ROIsum=70
  real(prec) :: hermiteh_gh(0:hn,ggn),ghabsciss(ghn),ghweights(ghn)!,ggabsciss(:),ggweights(:),hermiteh_gg(4,0:hn,ggn)
  real(prec) :: dzejmu(4,-5:tho),norm(0:tho)
  real(prec) :: factorials(0:tho),binomials(0:tho,0:tho),coo(0:hn,0:hn)
  real(prec), parameter :: sqrtpi&
       &=1.77245385090551602729816748334114518_prec
  real(prec), parameter :: sqrttwo&
       &=1.41421356237309504880168872420969807857_prec

contains

  subroutine make_vecP(m,n,vecp,intinterm12)
    implicit none
    integer    :: m,n,i,j,a,b,al,ij
    real(prec) :: beta,temp
    real(prec) :: vecp(:,:,:), intinterm12(0:newhn,0:newhn,0:newhn,0:newhn)
    vecP = 0._prec
    beta = 2.5_prec
    do i = 0, m-1
       do j = 0, m-1
          do a = 0, n-1
             do b = 0, a
                if ((modulo(i+j+a+b,2).EQ.0)) then
                   temp = 0._prec
                   ij = i*m+j+1
                   do al = 0, n-1
                      temp = temp + inttildef12(al,b,a,j,al,i,beta,intinterm12)
                   end do
                   vecP(ij,a,b) = temp
                   if (a.NE.b) vecP(ij,b,a) = temp
                end if
             end do
          end do
       end do
    end do
  end subroutine make_vecP

  subroutine make_p1(m,n,matps,matpt,eorb,matpes,matpet)
    implicit none
    integer    :: m,n
    integer    :: i,j,k,l,al,p,ij,ji,kl,lk
    real(prec) :: temp,temp2,tempe,tempe2
    real(prec) :: tempp
    real(prec) :: matps(:,:), matpt(:,:), matpes(:,:), matpet(:,:)
    real(prec) :: eorb(n)
    matps = 0._prec
    matpt = 0._prec
    do l = 0, m-1
       do k = 0, m-1
          kl=k*m+l+1
          lk=l*m+k+1
          do j = 0, m-1
             do i = 0, m-1
                ij=i*m+j+1
                ji=j*m+i+1
                temp2 = 0._prec
                tempe2 = 0._prec
                ! tempp = 0._prec
                do p = 0, ROIsum
                   temp = 0._prec
                   tempe = 0._prec
                   do al = 0, n-1
                      temp = temp + intf12(i,j,al,p,3._prec)*intf12(al,p,k,l,3._prec)
                      tempe = tempe + temp*eorb(al+1)
                   end do
                   temp2 = temp2 + temp
                   tempe2 = tempe2 + tempe
                   ! ! ROI check
                   ! tempp = tempp + intf12(p,2,1,2,3._prec)*intf12(3,p,3,1,3._prec)
                   ! print*, p,tempp
                end do
                ! print*, tempp
                ! stop
                matps(ij,kl) = matps(ij,kl) + temp2 !there will be problem with which indices come from
                ! !!!tau_i and which from tau_j... fuck
                matps(ji,kl) = matps(ji,kl) + temp2
                matps(ij,lk) = matps(ij,lk) + temp2
                matps(ji,lk) = matps(ji,lk) + temp2
                matpt(ij,kl) = matpt(ij,kl) + temp2
                matpt(ji,kl) = matpt(ji,kl) - temp2
                matpt(ij,lk) = matpt(ij,lk) + temp2
                matpt(ji,lk) = matpt(ji,lk) - temp2
                matpes(ij,kl) = matpes(ij,kl) + tempe2
                matpes(ji,kl) = matpes(ji,kl) + tempe2
                matpes(ij,lk) = matpes(ij,lk) + tempe2
                matpes(ji,lk) = matpes(ji,lk) + tempe2
                matpet(ij,kl) = matpet(ij,kl) + tempe2
                matpet(ji,kl) = matpet(ji,kl) - tempe2
                matpet(ij,lk) = matpet(ij,lk) + tempe2
                matpet(ji,lk) = matpet(ji,lk) - tempe2
             end do
          end do
       end do
    end do
  end subroutine make_p1

  subroutine make_mats(m,n,int6f122,matjs,matjt,matks,matkt,matms,matmt)
    implicit none
    integer :: m,n
    integer :: i,j,k,l,al,be
    integer :: ii,jj,kk,ll,all,bee
    real(prec) :: val1,val2,val3,val4
    real(prec) :: int6f122(0:m-1,0:m-1,0:n-1,0:n-1,0:m-1,0:m-1)
    real(prec) :: matjs(:,:,:,:), matjt(:,:,:,:)
    real(prec) :: matks(:,:,:,:), matkt(:,:,:,:)
    real(prec) :: matms(:,:,:,:), matmt(:,:,:,:)
    matjs=0._prec
    matjt=0._prec
    matks=0._prec
    matkt=0._prec
    matms=0._prec
    matmt=0._prec
    do ll = 0, m-1
       do kk = 0, m-1
          do bee = 0, n-1
             do all = 0, n-1
                do jj = 0, m-1
                   do ii = 0, m-1
                      if (jj.GT.ii) then
                         j=jj
                         i=ii
                      else
                         j=ii
                         i=jj
                      end if
                      if (ll.GT.kk) then
                         l=ll
                         k=kk
                      else
                         l=kk
                         k=ll
                      end if
                      if (bee.GT.all) then
                         be=bee
                         al=all
                      else
                         be=all
                         al=bee
                      end if
                      val1=int6f122(i,j,al,be,k,l)
                      val2=int6f122(k,l,al,be,i,j)
                      val3=int6f122(i,l,al,be,k,j)
                      val4=int6f122(k,j,al,be,i,l)
                      matjs(i*m+j+1,k*m+l+1,al+1,be+1)=val1+val2+val3+val4
                      matjt(i*m+j+1,k*m+l+1,al+1,be+1)=val1+val2-val3-val4
                      matks(i*m+j+1,k*m+l+1,al+1,be+1)=val1+val2+val3+val4
                      matkt(i*m+j+1,k*m+l+1,al+1,be+1)=val1+val2-val3-val4
                      matms(i*m+j+1,k*m+l+1,al+1,be+1)=val1-val2+val3-val4
                      matmt(i*m+j+1,k*m+l+1,al+1,be+1)=val1-val2-val3+val4
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine make_mats

  subroutine read_int6f122_file(m,n,g,int6f122,fun)
    implicit none
    integer :: m,n,fun
    integer :: number,numberab,ii,i,j,k,l,al,be
    real(prec) :: lolzix,g
    real(prec) :: int6f122(0:m-1,0:m-1,0:n-1,0:n-1,0:m-1,0:m-1)
    open(11,file='bas/file_int6f122.F',status='unknown',form='unformatted',access='direct',RECL=RECB)
    ! if (modulo(m,2).EQ.1) then
    !    print*, "odd is not good" !check if this is proper? if this shit works for odd...
    !    stop
    ! else
    !    fun=int(0.25_prec*((n/2._prec)*(n/2+1)*m**2+(n*(n+1))/2._prec*m**3+(n*(n+1))/4._prec*m**4))
    ! end if
    int6f122=0.0_prec
    do ii = 1, fun
       read(11,rec=ii) number,numberab,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       k=iand(ishft(number,-16),255)
       l=ishft(number,-24)
       al=iand(numberab,255)
       be=iand(ishft(numberab,-8),255)
       lolzix = lolzix*g
       int6f122(i,j,al,be,k,l) = lolzix
       int6f122(j,i,al,be,k,l) = lolzix
       int6f122(i,j,be,al,k,l) = lolzix
       int6f122(i,j,al,be,l,k) = lolzix
       int6f122(j,i,be,al,k,l) = lolzix
       int6f122(j,i,al,be,l,k) = lolzix
       int6f122(i,j,be,al,l,k) = lolzix
       int6f122(j,i,be,al,l,k) = lolzix
    end do
    close(11)
  end subroutine read_int6f122_file

  subroutine create_int6f122_file(m,n,intinterm22,caunta)
    implicit none
    integer :: al,be,i,j,k,l,ij,kl,n
    integer, intent(in) :: m
    integer, intent(out) :: caunta
    real(prec) :: integral,beta
    real(prec) :: intinterm22(0:newhn,0:newhn,0:newhn,0:newhn)
    beta = 4._prec
    open(11,file='bas/file_int6f122.F',status='unknown',form='unformatted',access='direct',RECL=RECB)
    caunta = 0
    do l = 0, m-1
       do k = 0, l
          do be = 0, n-1
             do al = 0, be
                do j = 0, m-1
                   do i = 0, j
                      if ((modulo(i+j+al+be+k+l,2).EQ.0)) then
                         integral=inttildef122(i,j,al,be,k,l,beta,intinterm22)
                         caunta = caunta + 1
                         write(11,rec=caunta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),al+ishft(be,8),integral
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do
    close(11)
  end subroutine create_int6f122_file

  subroutine read_int2_vec_file(twe,m,n)
    implicit none
    integer :: m,n,fun
    integer :: number,ii,i,j,k,l,jj
    real(prec) :: lolzix,mm,nn
    real(prec) :: twe(:,:,:)
    open(11,file='bas/file_int2_vec.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    nn=real(n,prec)
    if (modulo(m,2).EQ.1) then
       print*, "odd is not good" !check if this is proper? if this shit works for odd...
       stop
    else
       fun=int(mm*(mm+1._prec)*nn*nn/4._prec)
    end if
    twe=0.0_prec
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       k=iand(ishft(number,-16),255)
       l=ishft(number,-24)
       jj=k*m+l+1
       twe(jj,i+1,j+1) = lolzix
       jj=l*m+k+1
       twe(jj,i+1,j+1) = lolzix
    end do
    close(11)
  end subroutine read_int2_vec_file

  subroutine create_int2_vec_file(m,n)
    implicit none
    integer :: i,j,k,l,ij,kl,counter,nm,n
    integer, intent(in) :: m
    real(prec) :: integral,beta
    beta = 3._prec
    open(11,file='bas/file_int2_vec.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do l = 0, nm
       do k = 0, l
          kl = k*(k+1)/2+l
          do j = 0, n-1
             do i = 0, n-1
                ij = i*(i+1)/2+j
                if ((modulo(modulo(i,2)+modulo(j,2)+modulo(k,2)+modulo(l,2),2).EQ.0)) then
                   integral=intf12chem(i,j,k,l,beta) !check if this is proper?
                   counter = counter + 1
                   write(11,rec=counter) i+ishft(j,8)+ishft(k,16)+ishft(l,24),integral
                end if
             end do
          end do
       end do
    end do
    close(11)
  end subroutine create_int2_vec_file

  subroutine read_matH_file(ABtype,m,twe)
    implicit none
    character(1),intent(in) :: ABtype
    integer :: m,fun
    integer :: number,ii,i,j,k,l,jj,kk
    real(prec) :: lolzix,mm
    real(prec) :: twe(:,:)
    open(11,file='bas/fileH.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    if (modulo(m,2).EQ.1) then !check if this is proper? if this shit works for odd...
       fun=int((1._prec/8._prec)*((mm)**4._prec+2._prec*(mm)**3._prec+2._prec*(mm)**2._prec+2._prec*(m)+1))
    else
       fun=int((1._prec/8._prec)*((mm)**4._prec+2._prec*(mm)**3._prec+2._prec*(mm)**2._prec))
    end if
    twe=0.0_prec
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       if (ABtype.eq."A") then
          i=iand(number,255)
          j=iand(ishft(number,-8),255)
          k=iand(ishft(number,-16),255)
          l=ishft(number,-24)
       else if (ABtype.eq."B") then
          i=iand(number,255)
          l=iand(ishft(number,-8),255)
          k=iand(ishft(number,-16),255)
          j=ishft(number,-24)
       else
          print*, "Something's gone shitty with the ABtype"
          call kurwout
       end if
       jj=i*m+j+1
       kk=k*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+j+1
       kk=l*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=j*m+i+1
       kk=k*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=j*m+i+1
       kk=l*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
    end do
    close(11)
  end subroutine read_matH_file

  subroutine create_matH_file(m)
    implicit none
    integer :: i,j,k,l,ij,kl,counter,nmd,nm
    integer, intent(in) :: m
    real(prec) :: integral,beta
    beta = 5._prec !check this shit in the future. Is it really 5?
    open(11,file='bas/fileH.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    nmd = (nm+1)/2
    counter = 0
    do l = 0, nm
       do k = 0, l
          kl = k*(k+1)/2+l
          do j = 0, nm
             do i = 0, j
                ij = i*(i+1)/2+j
                if ((modulo(modulo(i,2)+modulo(j,2)+modulo(k,2)+modulo(l,2),2).EQ.0)) then
                   integral=0.5_prec*(intdf122chem(i,j,k,l,beta)+(i+j+k+l+2)*intf122chem(i,j,k,l,beta))
                   !check derivation to be sure
                   counter = counter + 1
                   write(11,rec=counter) i+ishft(j,8)+ishft(k,16)+ishft(l,24),integral
                end if
             end do
          end do
       end do
    end do
    close(11)
  end subroutine create_matH_file

  subroutine read_matS_file(ABtype,m,twe)
    implicit none
    character(1),intent(in) :: ABtype
    integer :: m,fun
    integer :: number,ii,i,j,k,l,jj,kk
    real(prec) :: lolzix,mm
    real(prec) :: twe(:,:)
    open(11,file='bas/fileS.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    if (modulo(m,2).EQ.1) then
       fun=int((1._prec/8._prec)*((mm)**4._prec+2._prec*(mm)**3._prec+2._prec*(mm)**2._prec+2._prec*(m)+1))
    else
       fun=int((1._prec/8._prec)*((mm)**4._prec+2._prec*(mm)**3._prec+2._prec*(mm)**2._prec))
    end if
    twe=0.0_prec
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       if (ABtype.eq."A") then
          i=iand(number,255)
          j=iand(ishft(number,-8),255)
          k=iand(ishft(number,-16),255)
          l=ishft(number,-24)
       else if (ABtype.eq."B") then
          i=iand(number,255)
          l=iand(ishft(number,-8),255)
          k=iand(ishft(number,-16),255)
          j=ishft(number,-24)
       else
          print*, "Something's gone shitty with the ABtype"
          call kurwout
       end if
       jj=i*m+j+1
       kk=k*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+j+1
       kk=l*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=j*m+i+1
       kk=k*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=j*m+i+1
       kk=l*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
    end do
    close(11)
  end subroutine read_matS_file

  subroutine create_matS_file(m)
    implicit none
    integer :: i,j,k,l,ij,kl,counter,nmd,nm
    integer, intent(in) :: m
    real(prec) :: integral,beta
    beta = 5._prec
    open(11,file='bas/fileS.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    nmd = (nm+1)/2
    counter = 0
    do l = 0, nm
       do k = 0, l
          kl = k*(k+1)/2+l
          do j = 0, nm
             do i = 0, j
                ij = i*(i+1)/2+j
                if ((modulo(modulo(i,2)+modulo(j,2)+modulo(k,2)+modulo(l,2),2).EQ.0)) then
                   integral=intf122chem(i,j,k,l,beta)
                   counter = counter + 1
                   write(11,rec=counter) i+ishft(j,8)+ishft(k,16)+ishft(l,24),integral
                end if
             end do
          end do
       end do
    end do
    close(11)
  end subroutine create_matS_file

  subroutine read_SCF_file(m,twe,gie)
    implicit none
    integer :: m,fun
    integer :: number,ii,i,j,k,l,jj,kk
    real(prec) :: lolzix,gie
    real(prec) :: twe(:,:)
    open(10,file='bas/file2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    fun=int(1._prec/96._prec*(81._prec+15._prec*(-1._prec)**(m-1._prec)+2._prec*(65._prec+3._prec*&
         (-1._prec)**(m-1._prec))*(m-1._prec)+76._prec*(m-1._prec)**2._prec+20._prec*(m-1._prec)**3._prec&
         +2._prec*(m-1._prec)**4._prec))
    do ii = 1, fun
       read(10,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       k=iand(ishft(number,-16),255)
       l=ishft(number,-24)
       lolzix = lolzix * gie
       jj=i*m+j+1
       kk=k*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=j*m+i+1
       kk=k*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+j+1
       kk=l*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=j*m+i+1
       kk=l*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+k+1
       kk=j*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=k*m+i+1
       kk=j*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+k+1
       kk=l*m+j+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=k*m+i+1
       kk=l*m+j+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+l+1
       kk=k*m+j+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+l+1
       kk=j*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=l*m+i+1
       kk=k*m+j+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=l*m+i+1
       kk=j*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
    end do
    close(10)
  end subroutine read_SCF_file

  subroutine create_SCF_file(m)
    implicit none
    integer :: i,j,k,l,ij,kl,counter,nmd,nm
    integer, intent(in) :: m
    open(10,file='bas/file2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    nmd = (nm+1)/2
    counter = 0
    do l = 0, nm
       do k = 0, l
          kl = k*(k+1)/2+l
          do j = 0, k
             do i = 0, j
                ij = i*(i+1)/2+j
                if ((modulo(modulo(i,2)+modulo(j,2)+modulo(k,2)+modulo(l,2),2).EQ.0)) then
                   counter = counter + 1
                   write(10,rec=counter) i+ishft(j,8)+ishft(k,16)+ishft(l,24),intgh(i,j,k,l,2._prec)
                end if
             end do
          end do
       end do
    end do
    close(10)
  end subroutine create_SCF_file

  real(prec) function coeff(s,i,k)
    implicit none
    integer :: i,k,s
    coeff = factorials(s)*2._prec**real(s,prec)*binomials(i,s)*binomials(k,s)
  end function coeff

  ! subroutine import_gg
  !   implicit none
  !   integer            :: i
  !   open(10,file='dat/ggweights.dat')
  !   ggweights=0._prec
  !   do i = 1, ggn
  !      read(10,*) ggweights(i)
  !   end do
  !   close(10)

  !   open(10,file='dat/ggabsciss.dat')
  !   ggabsciss=0._prec
  !   do i = 1, ggn
  !      read(10,*) ggabsciss(i)
  !   end do
  !   close(10)
  ! end subroutine import_gg

  subroutine import_gh
    implicit none
    integer            :: i
    open(10,file='dat/ghabsciss.dat')
    ghabsciss=0._prec
    do i = 1, ghn
       read(10,*) ghabsciss(i)
    end do
    close(10)

    open(10,file='dat/ghweights.dat')
    ghweights=0._prec
    do i = 1, ghn
       read(10,*) ghweights(i)
    end do
    close(10)
  end subroutine import_gh

  ! subroutine gener_hermiteh_gg(beta)
  !   implicit none
  !   integer         :: i,j,b
  !   real(prec)      :: beta
  !   b = bet(beta)
  !   hermiteh_gg(b,:,:)=0._prec
  !   do i = 1, ggn
  !      hermiteh_gg(b,0,i)=1._prec
  !      hermiteh_gg(b,1,i)=2._prec*ggabsciss(i)/sqrt(beta)
  !   end do
  !   do j = 2, hn
  !      do i = 1, ggn
  !         hermiteh_gg(b,j,i) = 2._prec*ggabsciss(i)/sqrt(beta)*hermiteh_gg(b,j-1,i)-2._prec*(j-1)*hermiteh_gg(b,j-2,i)
  !      end do
  !   end do
  ! end subroutine gener_hermiteh_gg

  subroutine gener_hermiteh_gh(beta)
    implicit none
    integer         :: i,j
    real(prec)      :: beta
    hermiteh_gh=0._prec
    do i = 1, ghn
       hermiteh_gh(0,i)=1._prec
       hermiteh_gh(1,i)=2._prec*ghabsciss(i)/sqrt(beta)
    end do
    do j=2,hn
       do i = 1, ghn
          hermiteh_gh(j,i) = 2._prec*ghabsciss(i)/sqrt(beta)*hermiteh_gh(j-1,i)-2._prec*(j-1)*hermiteh_gh(j-2,i)
       end do
    end do
  end subroutine gener_hermiteh_gh

  function sortt(vec) result(vecet)
    implicit none
    integer               :: n,j,temp,bubble
    integer, dimension(4) :: vec,vecet
    n = 4
    do while (n.GT.1)
       bubble = 0
       do j = 1, (n-1)
          if (vec(j)>vec(j+1)) then
             temp = vec(j)
             vec(j) = vec(j+1)
             vec(j+1) = temp
             bubble = j
          end if
       end do
       n = bubble
    end do
    vecet = vec
  end function sortt

  subroutine make_norm()
    implicit none
    integer :: m
    ! norm(:) = 0._prec
    do m = 0, tho
       norm(m) = fnorm(m)
    end do
  end subroutine make_norm

  real(prec) function fnorm(m)
    implicit none
    integer :: m
    fnorm = 1._prec/sqrt(2._prec**m*factorials(m)*sqrtpi)
  end function fnorm

  real(prec) function abn(m,n,pm)
    implicit none
    integer    :: m,n,i
    ! real(prec) :: abn
    logical    :: pm
    abn = 1._prec
    if (n.LE.m) then
       if (pm) then
          do i = 1, n
             abn = abn * real(m+i,prec)
          end do
          abn = sqrt(abn) * sqrt(2._prec**real(n,prec))
       else
          do i = 1, n
             abn = abn / real(m-i+1,prec)
          end do
          abn = sqrt(abn) / sqrt(2._prec**real(n,prec))
       end if
    end if
  end function abn

  subroutine import_coo
    implicit none
    integer :: iter,i,n
    open(10,file='dat/coo.dat')
    coo=0._prec
    do iter = 1, 22951
       read(10,*) i,n,coo(i,n)
    end do
    close(10)
  end subroutine import_coo

  subroutine import_factorials
    implicit none
    integer            :: i
    open(10,file='dat/factorials.dat')
    factorials=0._prec
    do i = 0, tho
       read(10,*) factorials(i)
    end do
    close(10)
  end subroutine import_factorials

  subroutine import_binomials
    implicit none
    integer            :: i,j,counter
    open(10,file='dat/binomials.dat')
    binomials=0._prec
    do i = 0, tho
       do j = 0, i
          read(10,*) binomials(i,j)
       end do
    end do
    close(10)
  end subroutine import_binomials

  real(prec) function intgh(mm,oo,pp,rr,beta)
    implicit none
    integer    :: m,o,p,r,i
    integer    :: mm,oo,pp,rr,vec(4)
    real(prec) :: beta
    intgh = 0._prec
    if (modulo(modulo(mm,2)+modulo(oo,2)+modulo(pp,2)+modulo(rr,2),2).EQ.0) then
       vec(1) = mm
       vec(2) = oo
       vec(3) = pp
       vec(4) = rr
       vec = sortt(vec)
       m = vec(1)
       o = vec(2)
       p = vec(3)
       r = vec(4)
       do i = 1, ghn
          intgh = intgh + hermiteh_gh(m,i)*hermiteh_gh(o,i)*hermiteh_gh(p,i)*hermiteh_gh(r,i) * ghweights(i)
       end do
       intgh = intgh*norm(m)*norm(o)*norm(p)*norm(r)/sqrt(beta)
    end if
  end function intgh

  real(prec) function intf12(mm,oo,pp,rr,beta)
    implicit none
    integer    :: m,o,p,r,s,t
    integer    :: mm,oo,pp,rr
    real(prec) :: tempf12,beta
    intf12 = 0._prec
    if (modulo(m+o+p+r,2).EQ.0) then
       if (pp.lt.mm) then
          m=pp
          p=mm
       else
          m=mm
          p=pp
       end if
       if (rr.lt.oo) then
          r=oo
          o=rr
       else
          r=rr
          o=oo
       end if
       do s = 0, m
          tempf12 = 0._prec
          do t = 0, o
             tempf12 = tempf12 + coeff(t,o,r) * imunuf12norm(m+p-2*s,o+r-2*t,beta)/norm(o+r-2*t)
          end do
          intf12 = intf12 + coeff(s,m,p) * tempf12/norm(m+p-2*s)
       end do
       intf12 = norm(mm)*norm(oo)*norm(pp)*norm(rr)*intf12
    end if
  end function intf12

  real(prec) function intf12chem(mm,oo,pp,rr,beta)
    implicit none
    integer, intent(in)    :: mm,oo,pp,rr
    real(prec)             :: beta
    intf12chem = intf12(mm,pp,oo,rr,beta)
  end function intf12chem

  real(prec) function intf122(mm,oo,pp,rr,beta)
    implicit none
    integer    :: m,o,p,r,s,t
    integer    :: mm,oo,pp,rr
    real(prec) :: tempf12,beta
    intf122 = 0._prec
    if (modulo(modulo(mm,2)+modulo(oo,2)+modulo(pp,2)+modulo(rr,2),2).EQ.0) then
       if (pp.lt.mm) then
          m=pp
          p=mm
       else
          m=mm
          p=pp
       end if
       if (rr.lt.oo) then
          r=oo
          o=rr
       else
          r=rr
          o=oo
       end if
       do s = 0, m
          tempf12 = 0._prec
          do t = 0, o
             tempf12 = tempf12 + coeff(t,o,r) * imunuf122norm(m+p-2*s,o+r-2*t,beta)/norm(o+r-2*t)
          end do
          intf122 = intf122 + coeff(s,m,p) * tempf12/norm(m+p-2*s)
       end do
       intf122 = norm(mm)*norm(oo)*norm(pp)*norm(rr)*intf122
    end if
  end function intf122

  real(prec) function intf122chem(mm,oo,pp,rr,beta)
    implicit none
    integer, intent(in)    :: mm,oo,pp,rr
    real(prec)             :: beta
    intf122chem = intf122(mm,pp,oo,rr,beta)
  end function intf122chem

  real(prec) function intdf122(mm,oo,pp,rr,beta) !check this mofo -- 1,1,2,2 works, but 1,2,1,2 also???
    implicit none
    integer    :: m,o,p,r,s,t
    integer    :: mm,oo,pp,rr
    real(prec) :: tempf12,beta
    intdf122 = 0._prec
    if (modulo(mm+oo+pp+rr,2).EQ.0) then
       if (pp.lt.mm) then
          m=pp
          p=mm
       else
          m=mm
          p=pp
       end if
       if (rr.lt.oo) then
          r=oo
          o=rr
       else
          r=rr
          o=oo
       end if
       do s = 0, m
          tempf12 = 0._prec
          do t = 0, o
             tempf12 = tempf12 + coeff(t,o,r) * imunudf122norm(m+p-2*s,o+r-2*t,beta)/norm(o+r-2*t)
          end do
          intdf122 = intdf122 + coeff(s,m,p) * tempf12/norm(m+p-2*s)
       end do
       intdf122 = norm(mm)*norm(oo)*norm(pp)*norm(rr)*intdf122
    end if
  end function intdf122

  real(prec) function intdf122chem(mm,oo,pp,rr,beta)
    implicit none
    integer, intent(in)    :: mm,oo,pp,rr
    real(prec), intent(in) :: beta
    ! real(prec)             :: intdf122chem
    intdf122chem = intdf122(mm,pp,oo,rr,beta)
  end function intdf122chem

  real(prec) function oldimunudf122(m,n,beta)
    integer    :: m,n,b
    real(prec) :: beta
    oldimunudf122 = 0._prec
    if (mod(m+n,2).NE.1) then
       b = bet(beta)
       oldimunudf122 = 2._prec*sqrtpi*(-1._prec)**n*2._prec**(-(m+n)/2._prec)*(&
            &dzejmu(b,m+n)&
            &-(real(b,prec)-1._prec)/2._prec*( dzejmu(b,m+n+2)+&
            &real(4._prec*(m+n)+2,prec)*dzejmu(b,m+n)+&
            &real(4._prec*(m+n)*(m+n-1),prec)*dzejmu(b,m+n-2) )&
            &+((real(b,prec)-1._prec)/4._prec)**2*&
            &( dzejmu(b,m+n+4)+&
            & real(8._prec*(m+n)+12._prec,prec)*dzejmu(b,m+n+2)+&
            & real(24._prec*(m+n)**2+24._prec*(m+n)+12._prec,prec)*dzejmu(b,m+n)+&
            & real(32._prec*(m+n)**3-48._prec*(m+n)**2+16._prec*(m+n),prec)*dzejmu(b,m+n-2)+&
            & real(16._prec*(m+n)*(m+n-1)*(m+n-2)*(m+n-3),prec)*dzejmu(b,m+n-4)&
            &))
    end if
  end function oldimunudf122

  real(prec) function imunudf122(m,n,beta)
    integer    :: m,n,b
    real(prec) :: beta
    imunudf122 = 0._prec
    if (mod(m+n,2).NE.1) then
       b = bet(beta)
       imunudf122 = 2._prec*sqrtpi*(-1._prec)**n*2._prec**(-(m+n)/2._prec)*(dzejmu(b,m+n)-&
            & 2._prec*(beta-1._prec)/4._prec*(dzejmu(b,m+n+2)+&
            & real(4._prec*(m+n)+2,prec)*dzejmu(b,m+n)+&
            & real(4._prec*(m+n)**2-4._prec*(m+n),prec)*dzejmu(b,m+n-2))+&
            & ((beta-1._prec)/4._prec)**2*(dzejmu(b,m+n+4)+&
            & real(8._prec*(m+n)+12._prec,prec)*dzejmu(b,m+n+2)+&
            & real(24._prec*(m+n)**2+24._prec*(m+n)+12._prec,prec)*dzejmu(b,m+n)+&
            & real(32._prec*(m+n)**3-48._prec*(m+n)**2+16._prec*(m+n),prec)*dzejmu(b,m+n-2)+&
            & real(16._prec*(m+n)*(m+n-1)*(m+n-2)*(m+n-3),prec)*dzejmu(b,m+n-4)))
    end if
  end function imunudf122

  real(prec) function imunudf122norm(m,n,beta)
    integer    :: m,n,b
    real(prec) :: beta,mn
    !some super bogus code...
    if (m.EQ.0.and.n.eq.0) then
       imunudf122norm = 0.5390125124744174978149239120563448922870409268698824918081870965
    else if (m.EQ.2.and.n.eq.0.or.m.EQ.0.and.n.eq.2) then
       imunudf122norm = -0.1345197891935502887834707784128512214525037706050587983590004814
    else if (m.EQ.1.and.n.eq.1) then
       imunudf122norm = 0.1902397102850885286405613807257687855130732683070173500499483870
    else
       !/bogus
       mn = real(m+n,prec)
       imunudf122norm = 0._prec
       if (mod(m+n,2).NE.1) then
          b = bet(beta)
          imunudf122norm = norm(n)*norm(m)/norm(m+n)*&
               & 2._prec*sqrtpi*(-1._prec)**n*2._prec**(-(m+n)/2._prec)*&
               &(dzejmu(b,m+n)-&
               & 2._prec*(beta-1._prec)/4._prec*(dzejmu(b,m+n+2)*abn(m+n,2,.TRUE.)+&
               & (4._prec*(m+n)+2._prec)*dzejmu(b,m+n)+&
               & (4._prec*mn**2-4._prec*mn)*dzejmu(b,m+n-2)*abn(m+n,2,.FALSE.))+&
               & ((beta-1._prec)/4._prec)**2*(dzejmu(b,m+n+4)*abn(m+n,4,.TRUE.)+&
               & (8._prec*mn+12._prec)*dzejmu(b,m+n+2)*abn(m+n,2,.TRUE.)+&
               & (24._prec*mn**2+24._prec*mn+12._prec)*dzejmu(b,m+n)+&
               & (32._prec*mn**3-48._prec*mn**2+16._prec*mn)*dzejmu(b,m+n-2)*abn(m+n,2,.FALSE.)+&
               & (16._prec*mn*(mn-1._prec)*(mn-2._prec)*(mn-3._prec))*dzejmu(b,m+n-4)*abn(m+n,4,.FALSE.)))
       end if
    end if
  end function imunudf122norm

  real(prec) function imunuf122norm(mu,nu,beta)
    implicit none
    integer    :: mu,nu,b
    real(prec) :: beta
    imunuf122norm = 0._prec
    if (mod(mu+nu,2).NE.1) then
       b = bet(beta)
       imunuf122norm = norm(nu)*norm(mu)/norm(mu+nu)*4._prec*sqrtpi*(-1._prec)**nu*2._prec**(-(mu+nu)/2._prec)&
            &*(0.25_prec*dzejmu(b,mu+nu+2)*sqrt(4._prec*real(mu+nu+2)*real(mu+nu+1))&
            &+(real(mu+nu,prec)+0.5_prec)*dzejmu(b,mu+nu)&
            &+sqrt(real(mu+nu,prec)*real(mu+nu-1,prec))*dzejmu(b,mu+nu-2)/2._prec)
    end if
  end function imunuf122norm

  real(prec) function imunuf12norm(mu,nu,beta)
    implicit none
    integer    :: mu,nu,b
    real(prec) :: beta
    imunuf12norm = 0._prec
    if (mod(mu+nu,2).NE.1) then
       b = bet(beta)
       imunuf12norm = norm(nu)*norm(mu)/norm(mu+nu)*sqrttwo*sqrtpi*(-1._prec)**nu*2._prec**(-real(mu+nu,prec)/2._prec)&
            &*(sqrt(2._prec*real(mu+nu+1,prec))*dzejmu(b,mu+nu+1)+sqrt(2._prec*(real(mu+nu,prec)))*dzejmu(b,mu+nu-1))
    end if
  end function imunuf12norm

  real(prec) function imunutildef12norm(mu,nu,beta)
    implicit none
    integer    :: mu,nu,b
    real(prec) :: beta
    imunutildef12norm = 0._prec
    if (mod(mu+nu,2).NE.1) then
       b = bet(beta)
       imunutildef12norm = norm(nu)*norm(mu)/norm(mu+nu)*&
            &sqrtpi*2._prec**((real(nu,prec)-1._prec)/2._prec)*3._prec**(-(real(mu+nu,prec)-1._prec)/2._prec)*&
            &(-1._prec)**(real(nu,prec))*&
            &(sqrt(real(mu+nu+1,prec))*dzejmu(b,mu+nu+1)+sqrt(real(mu+nu,prec))*dzejmu(b,mu+nu-1))
    end if
  end function imunutildef12norm

  real(prec) function preintildef12(mu,nu,beta)
    implicit none
    integer    :: mu,nu,i
    real(prec) :: beta
    preintildef12 = 0._prec
    do i = 0, mu/2
       preintildef12 = preintildef12 + coo(i,mu)*imunutildef12norm(mu-2*i,nu,beta)/norm(mu-2*i)
    end do
    preintildef12 = preintildef12*norm(mu)
  end function preintildef12

  subroutine create_intildef12_file()
    implicit none
    integer    :: counter,i,j,nm,m
    m=newhn
    open(10,file='bas/file_t1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do j = 0, nm
       do i = 0, nm
          if ((modulo(i+j,2).EQ.0)) then
             counter = counter + 1
             write(10,rec=counter) i+ishft(j,8),preintildef12(i,j,2.5_prec)
          end if
       end do
    end do
    close(10)
  end subroutine create_intildef12_file

  subroutine read_intildef12_file(intildef12)
    implicit none
    integer             :: i,j,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intildef12(0:newhn,0:newhn)
    m=newhn
    intildef12 = 0._prec
    open(11,file='bas/file_t1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**2)/2)
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       intildef12(i,j) = lolzix
       ! print*, i,j,lolzix,preintildef12(i,j,2.5_prec)
    end do
    close(11)
  end subroutine read_intildef12_file

  subroutine create_interm1_file(intildef12)
    implicit none
    integer    :: counter,nm,m,s,t,u,other
    real(prec) :: temp
    real(prec) :: intildef12(0:newhn,0:newhn)
    m=newhn/2
    open(10,file='bas/file_t1i1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do other = 0, (m*2)-1
       do t = 0, nm
          do s = 0, t
             temp = 0._prec
             do u = 0, s
                if ((modulo(other+s+t-2*u,2).EQ.0)) then
                   ! write(10,rec=counter) i+ishft(j,8),preintildef12(i,j,2.5_prec)
                   temp = temp + coeff(u,s,t)*intildef12(other,s+t-2*u)/norm(s+t-2*u)
                end if
             end do
             counter = counter + 1
             write(10,rec=counter) other+ishft(s,8)+ishft(t,16),temp
          end do
       end do
    end do
    close(10)
  end subroutine create_interm1_file

  subroutine read_interm1_file(intinterm1)
    implicit none
    integer             :: other,s,t,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intinterm1(0:newhn,0:newhn,0:newhn)
    m=newhn/2
    intinterm1 = 0._prec
    open(11,file='bas/file_t1i1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**2+mm**3))
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       other=iand(number,255)
       s=iand(ishft(number,-8),255)
       t=iand(ishft(number,-16),255)
       intinterm1(other,s,t) = lolzix
    end do
    close(11)
  end subroutine read_interm1_file

  subroutine create_interm12_file(intinterm1)
    implicit none
    integer    :: counter,nm,m,w,i,j,s,t,summing
    real(prec) :: temp,val
    real(prec) :: intinterm1(0:newhn,0:newhn,0:newhn)
    m=newhn/2
    open(10,file='bas/file_t1i2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do t = 0, nm
       do s = 0, t
          do j = 0, nm
             do i = 0, nm
                temp = 0._prec
                if (i.LT.j) then
                   summing = i
                else
                   summing = j
                endif
                do w = 0, summing
                   ! val = intinterm1(i+j-2*w,s,t)
                   ! if (val.NE.0) then
                   if (modulo(i+j-2*w+s+t,2).EQ.0) then
                      ! write(10,rec=counter) i+ishft(j,8),preintildef12(i,j,2.5_prec)
                      temp = temp + coeff(w,i,j)*intinterm1(i+j-2*w,s,t)/norm(i+j-2*w)
                   end if
                end do
                counter = counter + 1
                write(10,rec=counter) i+ishft(j,8)+ishft(s,16)+ishft(t,24),temp
             end do
          end do
       end do
    end do
    close(10)
    ! print*, counter
  end subroutine create_interm12_file

  subroutine read_interm12_file(intinterm12)
    implicit none
    integer             :: i,j,s,t,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intinterm12(0:newhn,0:newhn,0:newhn,0:newhn)
    m=newhn/2
    intinterm12 = 0._prec
    open(11,file='bas/file_t1i2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**4+mm**3)/2._prec)
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       s=iand(ishft(number,-16),255)
       t=ishft(number,-24)
       intinterm12(i,j,s,t) = lolzix
    end do
    close(11)
  end subroutine read_interm12_file

  real(prec) function inttildef12(mm,nn,oo,pp,ss,tt,beta,intinterm12)
    implicit none
    integer    :: m,n,o,p,s,t,u,w,x,y
    integer    :: mm,nn,oo,pp,ss,tt
    integer    :: summing
    integer    :: v1(4)
    real(prec) :: tempf12,ttempf12,tttempf12,beta
    real(prec) :: intinterm12(0:newhn,0:newhn,0:newhn,0:newhn)
    inttildef12 = 0._prec
    if (modulo(mm+nn+oo+pp+ss+tt,2).EQ.0) then
       if (ss.lt.tt) then
          s=ss
          t=tt
       else
          s=tt
          t=ss
       end if
       v1 = [mm,nn,oo,pp]
       v1 = sortt(v1)
       m = v1(1)
       n = v1(2)
       o = v1(3)
       p = v1(4)
       ! do u = 0, s
       tempf12 = 0._prec
       do y = 0, o
          ttempf12 = 0._prec
          do x = 0, m
             ! tttempf12 = 0._prec
             ! if (m+n-2*x.LT.o+p-2*y) then
             !    summing = m+n-2*x
             ! else
             !    summing = o+p-2*y
             ! endif
             ! do w = 0, summing
                ! tttempf12 = tttempf12 + coeff(w,m+n-2*x,o+p-2*y)*intildef12(m+n-2*x+o+p-2*y-2*w,s+t-2*u)&
                !      &/norm(m+n-2*x+o+p-2*y-2*w)/norm(s+t-2*u)
             !    tttempf12 = tttempf12 + coeff(w,m+n-2*x,o+p-2*y)*intinterm1(m+n-2*x+o+p-2*y-2*w,s,t)&
             !         &/norm(m+n-2*x+o+p-2*y-2*w)
             ! end do
             ttempf12 = ttempf12 + coeff(x,m,n)*intinterm12(m+n-2*x,o+p-2*y,s,t)
          end do
          tempf12 = tempf12 + coeff(y,o,p)*ttempf12
       end do
       ! inttildef12 = inttildef12 + coeff(u,s,t) * tempf12
       ! end do
       inttildef12 = norm(mm)*norm(nn)*norm(oo)*norm(pp)*norm(ss)*norm(tt)*tempf12
    end if
  end function inttildef12

  real(prec) function imunutildef122norm(mu,nu,beta)
    implicit none
    integer    :: mu,nu,b
    real(prec) :: beta
    imunutildef122norm = 0._prec
    if (mod(mu+nu,2).NE.1) then
       b = bet(beta)
       imunutildef122norm = norm(nu)*norm(mu)/norm(mu+nu)*&
            &(sqrtpi*2._prec**((real(nu,prec)-1._prec)/2._prec)*3._prec**(-real(mu+nu,prec)/2._prec+1._prec)*&
            &(-1._prec)**(real(nu,prec)))*&
            &(2._prec*sqrt(real(mu+nu+2,prec)*real(mu+nu+1,prec))*&
            &0.25_prec*dzejmu(b,mu+nu+2)+(real(mu+nu,prec)+0.5_prec)*dzejmu(b,mu+nu)+&
            &0.5_prec*sqrt(real(mu+nu,prec)*(real(mu+nu,prec)-1._prec))*dzejmu(b,mu+nu-2))
    end if
  end function imunutildef122norm

  real(prec) function preintildef122(mu,nu,beta)
    implicit none
    integer    :: mu,nu,i
    real(prec) :: beta
    preintildef122 = 0._prec
    do i = 0, mu/2
       preintildef122 = preintildef122 + coo(i,mu)*imunutildef122norm(mu-2*i,nu,beta)/norm(mu-2*i)
    end do
    preintildef122 = preintildef122*norm(mu)
  end function preintildef122

  subroutine create_intildef122_file()
    implicit none
    integer    :: counter,i,j,nm,m
    m=newhn
    open(10,file='bas/file_t2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do j = 0, nm
       do i = 0, nm
          if ((modulo(i+j,2).EQ.0)) then
             counter = counter + 1
             write(10,rec=counter) i+ishft(j,8),preintildef122(i,j,4._prec)
          end if
       end do
    end do
    close(10)
  end subroutine create_intildef122_file

  subroutine read_intildef122_file(intildef122)
    implicit none
    integer             :: i,j,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intildef122(0:newhn,0:newhn)
    m=newhn
    intildef122 = 0._prec
    open(11,file='bas/file_t2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**2)/2)
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       intildef122(i,j) = lolzix
       ! print*, i,j,lolzix,preintildef12(i,j,2.5_prec)
    end do
    close(11)
  end subroutine read_intildef122_file

  subroutine create_interm2_file(intildef122)
    implicit none
    integer    :: counter,nm,m,s,t,u,other
    real(prec) :: temp
    real(prec) :: intildef122(0:newhn,0:newhn)
    m=newhn/2
    open(10,file='bas/file_t2i1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do other = 0, (m*2)-1
       do t = 0, nm
          do s = 0, t
             temp = 0._prec
             do u = 0, s
                if ((modulo(other+s+t-2*u,2).EQ.0)) then
                   ! write(10,rec=counter) i+ishft(j,8),preintildef12(i,j,2.5_prec)
                   temp = temp + coeff(u,s,t)*intildef122(other,s+t-2*u)/norm(s+t-2*u)
                end if
             end do
             counter = counter + 1
             write(10,rec=counter) other+ishft(s,8)+ishft(t,16),temp
          end do
       end do
    end do
    close(10)
  end subroutine create_interm2_file

  subroutine read_interm2_file(intinterm2)
    implicit none
    integer             :: other,s,t,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intinterm2(0:newhn,0:newhn,0:newhn)
    m=newhn/2
    intinterm2 = 0._prec
    open(11,file='bas/file_t2i1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**2+mm**3))
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       other=iand(number,255)
       s=iand(ishft(number,-8),255)
       t=iand(ishft(number,-16),255)
       intinterm2(other,s,t) = lolzix
    end do
    close(11)
  end subroutine read_interm2_file

  subroutine create_interm22_file(intinterm2)
    implicit none
    integer    :: counter,nm,m,w,i,j,s,t,summing
    real(prec) :: temp
    real(prec) :: intinterm2(0:newhn,0:newhn,0:newhn)
    m=newhn/2
    open(10,file='bas/file_t2i2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do t = 0, nm
       do s = 0, t
          do j = 0, nm
             do i = 0, nm
                temp = 0._prec
                if (i.LT.j) then
                   summing = i
                else
                   summing = j
                endif
                do w = 0, summing
                   if ((modulo(i+j-2*w+s+t,2).EQ.0)) then
                      temp = temp + coeff(w,i,j)*intinterm2(i+j-2*w,s,t)/norm(i+j-2*w)
                   end if
                end do
                counter = counter + 1
                write(10,rec=counter) i+ishft(j,8)+ishft(s,16)+ishft(t,24),temp
             end do
          end do
       end do
    end do
    close(10)
  end subroutine create_interm22_file

  subroutine read_interm22_file(intinterm22)
    implicit none
    integer             :: i,j,s,t,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intinterm22(0:newhn,0:newhn,0:newhn,0:newhn)
    m=newhn/2
    intinterm22 = 0._prec
    open(10,file='bas/file_t2i2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**4+mm**3)/2._prec)
    do ii = 1, fun
       read(10,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       s=iand(ishft(number,-16),255)
       t=ishft(number,-24)
       intinterm22(i,j,s,t) = lolzix
    end do
    close(10)
  end subroutine read_interm22_file

  real(prec) function inttildef122(mm,nn,oo,pp,ss,tt,beta,intinterm22)
    implicit none
    integer    :: m,n,o,p,s,t,u,w,x,y
    integer    :: mm,nn,oo,pp,ss,tt
    integer    :: summing
    integer    :: v1(4)
    real(prec) :: tempf12,ttempf12,tttempf12,beta
    real(prec) :: intinterm22(0:newhn,0:newhn,0:newhn,0:newhn)
    inttildef122 = 0._prec
    if (modulo(mm+nn+oo+pp+ss+tt,2).EQ.0) then
       if (ss.lt.tt) then
          s=ss
          t=tt
       else
          s=tt
          t=ss
       end if
       v1 = [mm,nn,oo,pp]
       v1 = sortt(v1)
       m = v1(1)
       n = v1(2)
       o = v1(3)
       p = v1(4)
       ! do u = 0, s
       tempf12 = 0._prec
       do y = 0, o
          ttempf12 = 0._prec
          do x = 0, m
             ! tttempf12 = 0._prec
             ! if (m+n-2*x.LT.o+p-2*y) then
             !    summing = m+n-2*x
             ! else
             !    summing = o+p-2*y
             ! endif
             ! do w = 0, summing
                ! tttempf12 = tttempf12 + coeff(w,m+n-2*x,o+p-2*y)*intildef12(m+n-2*x+o+p-2*y-2*w,s+t-2*u)&
                !      &/norm(m+n-2*x+o+p-2*y-2*w)/norm(s+t-2*u)
             !    tttempf12 = tttempf12 + coeff(w,m+n-2*x,o+p-2*y)*intinterm1(m+n-2*x+o+p-2*y-2*w,s,t)&
             !         &/norm(m+n-2*x+o+p-2*y-2*w)
             ! end do
             ttempf12 = ttempf12 + coeff(x,m,n)*intinterm22(m+n-2*x,o+p-2*y,s,t)
          end do
          tempf12 = tempf12 + coeff(y,o,p)*ttempf12
       end do
       ! inttildef12 = inttildef12 + coeff(u,s,t) * tempf12
       ! end do
       inttildef122 = norm(mm)*norm(nn)*norm(oo)*norm(pp)*norm(ss)*norm(tt)*tempf12
    end if
  end function inttildef122

  ! subroutine make_dzejmu(beta)
  !   implicit none
  !   real(prec) :: beta
  !   integer    :: i,mu,b
  !   b = bet(beta)
  !   dzejmu(b,:)=0._prec
  !   do mu = 0, hn
  !      do i = 1, ggn
  !         dzejmu(b,mu) = dzejmu(b,mu) + hermiteh_gg(b,mu,i) * ggweights(i)
  !      end do
  !      dzejmu(b,mu) = 1._prec/sqrt(beta) * dzejmu(b,mu) * norm(mu)
  !   end do
  ! end subroutine make_dzejmu

  subroutine import_dzejmu()
    implicit none
    real(prec) :: beta
    integer    :: i,mu,b
    open(10,file='bas/dz_2.5.dat')
    b = bet(2.5_prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
    open(10,file='bas/dz_3.dat')
    b = bet(3._prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
    open(10,file='bas/dz_4.dat')
    b = bet(4._prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
    open(10,file='bas/dz_5.dat')
    b = bet(5._prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
    ! print*, abs(dzejmu(1,1)-0.177042).GT.1.e-5
    ! print*, abs(dzejmu(2,1)-0.106225).GT.1.e-5
    ! stop
  end subroutine import_dzejmu

  integer function bet(beta)
    real(prec) :: beta
    if (abs(beta-3._prec).LT.1e-8) then
       bet = 1
    else if (abs(beta-5._prec).LT.1e-8) then
       bet = 2
    else if (abs(beta-2.5_prec).LT.1e-8) then
       bet = 3
    else if (abs(beta-4._prec).LT.1e-8) then
       bet = 4
    else
       print*, "Unsupported beta. Bombing out."
       call kurwout
    end if
  end function bet

  subroutine check(what)
    implicit none
    character(len=*), intent(in) :: what
    real(prec), parameter :: thr = 1.e-5
    select case (what)
    case ('dzejmu')
       if (abs(dzejmu(1,1)-0.333333).LT.thr) then
          print*, "dzejmu is probably not normalized, issues??"
       else if (abs(dzejmu(1,1)-0.177042).GT.thr.or.abs(dzejmu(2,1)-0.106225).GT.thr) then
          print*, abs(dzejmu(1,1)-0.177042),abs(dzejmu(2,1)-0.106225)
          print*, "normalized dzejmu is wrong"
          call kurwout()
       else
          write(LOUT,'(a)') "    Dzejmu is correct"
       end if
    case ('intf12')
       if (abs(intf12(3,3,3,3,3._prec)-0.126686006).GT.thr.or.abs(intf12(1,2,3,4,3._prec)-0.0280728333).GT.thr) then
          print*, "intf12 doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    Intf12 is correct"
       end if
    case ('intf122')
       if (abs(intf122(3,3,3,3,5._prec)-0.03901133673).GT.thr.or.abs(intf122(1,2,3,4,5._prec)-0.011738373).GT.thr) then
          print*, "intf122 doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    Intf122 is correct"
       end if
    case ('imunudf122norm')
       if (abs(imunudf122norm(3,3,5._prec)-0.086242001).GT.thr) then
          print*, "imunudf122norm doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    Imunudf122norm is correct"
       end if
    case ('intdf122')
       if (abs(intdf122(3,3,3,3,5._prec)-0.1662632851030596).GT.thr.or.&
            &abs(intdf122(1,2,3,4,5._prec)-0.03302934860107639).GT.thr) then
          print*, "intf122 doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    Intdf122 is correct"
       end if
    case ('imunutildef122norm')
       if (abs(imunutildef122norm(2,2,4._prec)+0.007343637523333).GT.thr.or.&
            &abs(imunutildef122norm(1,1,4._prec)-0.013847295710).GT.thr) then
          print*, "imunutildef122norm doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    imunutildef122 is correct"
       end if
    case ('imunutildef12norm')
       if (abs(imunutildef12norm(2,2,2.5_prec)+0.020014809331908905).GT.thr.or.&
            &abs(imunutildef12norm(1,1,2.5_prec)-0.03265986323710979).GT.thr) then
          print*, "imunutildef12norm doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    imunutildef12 is correct"
       end if
    ! case ('intildef12')
    !    if (abs(intildef12(3,3)+0.020425369523330).GT.thr.or.&
    !         &abs(intildef12(0,0)-0.3464101615137748).GT.thr) then
    !       print*, "intildef12 doesn't work"
    !       call kurwout()
    !    else
    !       write(LOUT,'(a)') "    intildef12 is correct"
    !    end if
    ! case ('inttildef12')
    !    if (abs(inttildef12(2,2,2,2,2,2,2.5_prec)-0.03541268534).GT.thr.or.&
    !         &abs(inttildef12(1,5,6,2,7,5,2.5_prec)+0.00090355492).GT.thr) then
    !       print*, "2,2,2,2,2,2",inttildef12(2,2,2,2,2,2,2.5_prec)
    !       print*, "1,5,6,2,7,5",inttildef12(1,5,6,2,7,5,2.5_prec)
    !       print*, "inttildef12 doesn't work"
    !       call kurwout()
    !    else
    !       write(LOUT,'(a)') "    inttildef12 is correct"
    !    end if
    ! case ('intildef122')
    !    if (abs(intildef122(3,3)+0.00826159).GT.thr.or.&
    !         &abs(intildef122(0,0)-0.117498).GT.thr) then
    !       print*, "intildef122 doesn't work"
    !       call kurwout()
    !    else
    !       write(LOUT,'(a)') "    intildef122 is correct"
    !    end if
    ! case ('inttildef122')
    !    if (abs(inttildef122(2,2,2,2,2,2,4._prec)-0.010651100993).GT.thr.or.&
    !         &abs(inttildef122(1,5,6,2,7,5,4._prec)+0.00039196641).GT.thr) then
    !       print*, "inttildef122 doesn't work"
    !       call kurwout()
    !    else
    !       write(LOUT,'(a)') "    inttildef122 is correct"
    !    end if
    end select
  end subroutine check

  ! function dzejmu(mu,beta)
  !   implicit none
  !   integer    :: i,mu
  !   real(prec) :: dzejmu,beta
  !   dzejmu = 0._prec
  !   if (mu.lt.0) then
  !      dzejmu = 0._prec
  !   else
  !      do i = 1, ggn
  !         dzejmu = dzejmu + hermiteh_gg(mu,i) * ggweights(i)
  !      end do
  !      dzejmu = 1._prec/sqrt(beta) * dzejmu * norm(mu)
  !   end if
  ! end function dzejmu

  subroutine kurwout()
    print*, "__/\\\________/\\\__/\\\________/\\\____/\\\\\\\\\______/\\\______________/\\\_____/\\\\\\\\\____       "
    print*, " _\/\\\_____/\\\//__\/\\\_______\/\\\__/\\\///////\\\___\/\\\_____________\/\\\___/\\\\\\\\\\\\\__      "
    print*, "  _\/\\\__/\\\//_____\/\\\_______\/\\\_\/\\\_____\/\\\___\/\\\_____________\/\\\__/\\\/////////\\\_     "
    print*, "   _\/\\\\\\//\\\_____\/\\\_______\/\\\_\/\\\\\\\\\\\/____\//\\\____/\\\____/\\\__\/\\\_______\/\\\_    "
    print*, "    _\/\\\//_\//\\\____\/\\\_______\/\\\_\/\\\//////\\\_____\//\\\__/\\\\\__/\\\___\/\\\\\\\\\\\\\\\_   "
    print*, "     _\/\\\____\//\\\___\/\\\_______\/\\\_\/\\\____\//\\\_____\//\\\/\\\/\\\/\\\____\/\\\/////////\\\_  "
    print*, "      _\/\\\_____\//\\\__\//\\\______/\\\__\/\\\_____\//\\\_____\//\\\\\\//\\\\\_____\/\\\_______\/\\\_ "
    print*, "       _\/\\\______\//\\\__\///\\\\\\\\\/___\/\\\______\//\\\_____\//\\\__\//\\\______\/\\\_______\/\\\_"
    print*, "        _\///________\///_____\/////////_____\///________\///_______\///____\///_______\///________\///_"
    stop
  end subroutine kurwout

end module ggint
