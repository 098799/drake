module ggint
  use precision, only : prec
  implicit none
  !    integer, parameter :: prec=16
  integer, parameter :: ggn=400,hn=200,ghn=100
  real(prec) :: ggabsciss(ggn),ggweights(ggn),hermiteh_gg(2,0:hn,ggn),hermiteh_gh(0:hn,ggn),ghabsciss(ghn),ghweights(ghn)
  real(prec) :: dzejmu(2,-4:hn)
  real(prec) :: factorials(0:hn),binomials(0:hn,0:hn)!sqrtfactorials(0:200)
  real(prec), parameter :: sqrtpi&
       &=1.77245385090551602729816748334114518_prec
  real(prec), parameter :: sqrttwo&
       &=1.41421356237309504880168872420969807857_prec

contains

  subroutine read_matS_file(ABtype,m,twe)
    implicit none
    character(1),intent(in) :: ABtype
    integer :: m,fun
    integer :: number,ii,i,j,k,l,jj,kk
    real(prec) :: lolzix,mm
    real(prec) :: twe(:,:)
    open(11,file='bas/fileS.F',status='unknown',form='unformatted',access='direct',RECL=24)
    mm=real(m,prec)
    if (modulo(m,2).EQ.1) then
       fun=int((1._prec/8._prec)*((mm)**4._prec+2._prec*(mm)**3._prec+2._prec*(mm)**2._prec+2._prec*(m)+1))
    else
       fun=int((1._prec/8._prec)*((mm)**4._prec+2._prec*(mm)**3._prec+2._prec*(mm)**2._prec))
    end if
    twe=13.33_prec
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
    open(11,file='bas/fileS.F',status='unknown',form='unformatted',access='direct',RECL=24)
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
    open(10,file='bas/file2.F',status='unknown',form='unformatted',access='direct',RECL=24)
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
    open(10,file='bas/file2.F',status='unknown',form='unformatted',access='direct',RECL=24)
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

  function coeff(i,k,s)
    implicit none
    integer :: i,k,s
    real(prec) :: coeff
    coeff = factorials(s)*2._prec**real(s,prec)*binomials(i,s)*binomials(k,s)
  end function coeff

  subroutine import_gg
    implicit none
    integer            :: i
    open(10,file='dat/ggweights.dat')  
    ggweights=0._prec
    do i = 1, ggn
       read(10,*) ggweights(i)
    end do
    close(10)

    open(10,file='dat/ggabsciss.dat')  
    ggabsciss=0._prec
    do i = 1, ggn
       read(10,*) ggabsciss(i)
    end do
    close(10)
  end subroutine import_gg

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

  subroutine gener_hermiteh_gg(beta)
    implicit none
    integer         :: i,j,b
    real(prec)      :: beta
    b = bet(beta)
    hermiteh_gg(b,:,:)=0._prec
    do i = 1, ggn
       hermiteh_gg(b,0,i)=1._prec
       hermiteh_gg(b,1,i)=2._prec*ggabsciss(i)/sqrt(beta)
    end do
    do j = 2, hn
       do i = 1, ggn
          hermiteh_gg(b,j,i) = 2._prec*ggabsciss(i)/sqrt(beta)*hermiteh_gg(b,j-1,i)-2._prec*(j-1)*hermiteh_gg(b,j-2,i)
       end do
    end do
  end subroutine gener_hermiteh_gg

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

  function norm(m)
    implicit none
    integer    :: m
    real(prec) :: norm
    norm = 1._prec/sqrt(2._prec**m*factorials(m)*sqrtpi)
  end function norm

  subroutine import_factorials
    implicit none
    integer            :: i
    open(10,file='dat/factorials.dat')  
    factorials=0._prec
    do i = 0, hn
       read(10,*) factorials(i)
    end do
    close(10)
  end subroutine import_factorials

  subroutine import_binomials
    implicit none
    integer            :: i,j
    open(10,file='dat/binomials.dat')  
    binomials=0._prec
    do i = 0, hn
       do j = 0, i
          read(10,*) binomials(i,j)
       end do
    end do
    close(10)
  end subroutine import_binomials

  function intgh(mm,oo,pp,rr,beta)
    implicit none
    integer    :: m,o,p,r,i
    integer    :: mm,oo,pp,rr,vec(4)
    real(prec) :: intgh,beta
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

  function intf12(mm,oo,pp,rr,beta)
    implicit none
    integer    :: m,o,p,r,s,t
    integer    :: mm,oo,pp,rr
    real(prec) :: intf12,tempf12,beta
    intf12 = 0._prec
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
             tempf12 = tempf12 + coeff(o,r,t) * imunuf12norm(m+p-2*s,o+r-2*t,beta)/norm(o+r-2*t)
          end do
          intf12 = intf12 + coeff(m,p,s) * tempf12/norm(m+p-2*s)
       end do
       intf12 = norm(mm)*norm(oo)*norm(pp)*norm(rr)*intf12
    end if
  end function intf12

  function intf12chem(mm,oo,pp,rr,beta)
    implicit none
    integer, intent(in)    :: mm,oo,pp,rr
    real(prec)             :: intf12chem,beta
    intf12chem = intf12(mm,pp,oo,rr,beta)
  end function intf12chem

  function intf122(mm,oo,pp,rr,beta)
    implicit none
    integer    :: m,o,p,r,s,t
    integer    :: mm,oo,pp,rr
    real(prec) :: intf122,tempf12,beta
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
             tempf12 = tempf12 + coeff(o,r,t) * imunuf122norm(m+p-2*s,o+r-2*t,beta)/norm(o+r-2*t)
          end do
          intf122 = intf122 + coeff(m,p,s) * tempf12/norm(m+p-2*s)
       end do
       intf122 = norm(mm)*norm(oo)*norm(pp)*norm(rr)*intf122
       ! else
       !    intf122 = 14.33_prec
    end if
  end function intf122

  function intf122chem(mm,oo,pp,rr,beta)
    implicit none
    integer, intent(in)    :: mm,oo,pp,rr
    real(prec)             :: intf122chem,beta
    intf122chem = intf122(mm,pp,oo,rr,beta)
  end function intf122chem

  function intdf122(mm,oo,pp,rr,beta)
    integer    :: m,o,p,r,s,t
    integer    :: mm,oo,pp,rr
    real(prec) :: intdf122,tempf12,beta
    intdf122 = 0._prec
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
             tempf12 = tempf12 + coeff(o,r,t) * imunudf122(m+p-2*s,o+r-2*t,beta)
          end do
          intdf122 = intdf122 + coeff(m,p,s) * tempf12
       end do
       intdf122 = norm(mm)*norm(oo)*norm(pp)*norm(rr)*intdf122
    end if
  end function intdf122

  function intdf122chem(mm,oo,pp,rr,beta)
    implicit none
    integer, intent(in)    :: mm,oo,pp,rr
    real(prec), intent(in) :: beta
    real(prec)             :: intdf122chem
    intdf122chem = intdf122(mm,pp,oo,rr,beta)
  end function intdf122chem

  function imunudf122(m,n,beta)
    integer    :: m,n,b
    real(prec) :: imunudf122,beta
    b = bet(beta)
    imunudf122 = 2._prec*sqrtpi*(-1._prec)**n*2._prec**(-(m+n)/2._prec)*(&
         &dzejmu(b,m+n)&
         &-(real(b,prec)-1._prec)/2._prec*&
         &( dzejmu(b,m+n+2)+real(4*(m+n)+2,prec)*dzejmu(b,m+n)+real(4*(m+n)*(m+n-1),prec)*dzejmu(b,m+n-2) )&
         &+((real(b,prec)-1._prec)/4._prec)**2*&
         &( dzejmu(b,m+n+4)+real(8*(m+n)+12,prec)*dzejmu(b,m+n+2)+real(24*(m+n)**2+24*(m+n)+12,prec)*dzejmu(b,m+n)+&
         &real(32*(m+n)**3-48*(m+n)**2+16*(m+n),prec)*dzejmu(b,m+n-2)+real(16*(m+n)*(m+n-1)*(m+n-2)*(m+n-3),prec)*&
         &dzejmu(b,m+n-4) )  )
  end function imunudf122

  function imunuf122norm(mu,nu,beta)
    implicit none
    integer    :: mu,nu,b
    real(prec) :: imunuf122norm,beta
    b = bet(beta)
    imunuf122norm = norm(nu)*norm(mu)/norm(mu+nu)*4._prec*sqrtpi*(-1._prec)**nu*2._prec**(-(mu+nu)/2._prec)&
         &*(0.25_prec*dzejmu(b,mu+nu+2)*sqrt(4._prec*real(mu+nu+2)*real(mu+nu+1))&
         &+(real(mu+nu,prec)+0.5_prec)*dzejmu(b,mu+nu)&
         &+sqrt(real(mu+nu,prec)*real(mu+nu-1,prec))*dzejmu(b,mu+nu-2)/2._prec)
  end function imunuf122norm

  function imunuf12norm(mu,nu,beta)
    implicit none
    integer    :: mu,nu,b
    real(prec) :: imunuf12norm,beta
    b = bet(beta)
    imunuf12norm = norm(nu)*norm(mu)/norm(mu+nu)*sqrttwo*sqrtpi*(-1._prec)**nu*2._prec**(-real(mu+nu,prec)/2._prec)&
         &*(sqrt(2._prec*real(mu+nu+1,prec))*dzejmu(b,mu+nu+1)+sqrt(2._prec*(real(mu+nu,prec)))*dzejmu(b,mu+nu-1))
  end function imunuf12norm

  subroutine make_dzejmu(beta)
    implicit none
    real(prec) :: beta
    integer    :: i,mu,b
    b = bet(beta)
    dzejmu(b,:)=0._prec
    print*, "making dzejmu of b=",b
    do mu = 0, hn
       do i = 1, ggn
          dzejmu(b,mu) = dzejmu(b,mu) + hermiteh_gg(b,mu,i) * ggweights(i)
       end do
       dzejmu(b,mu) = 1._prec/sqrt(beta) * dzejmu(b,mu)! * norm(mu)
    end do
  end subroutine make_dzejmu

  integer function bet(beta)
    real(prec) :: beta
    if (abs(beta-3._prec).LT.1e-8) then
       bet = 1
    else if (abs(beta-5._prec).LT.1e-8) then
       bet = 2
    else
       print*, "Unsupported beta. Bombing out."
       call kurwout
    end if
  end function bet

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

  subroutine kurwout
    print*, "__/\\\________/\\\__/\\\________/\\\____/\\\\\\\\\______/\\\______________/\\\_____/\\\\\\\\\____        "
    print*, " _\/\\\_____/\\\//__\/\\\_______\/\\\__/\\\///////\\\___\/\\\_____________\/\\\___/\\\\\\\\\\\\\__       "
    print*, "  _\/\\\__/\\\//_____\/\\\_______\/\\\_\/\\\_____\/\\\___\/\\\_____________\/\\\__/\\\/////////\\\_      "
    print*, "   _\/\\\\\\//\\\_____\/\\\_______\/\\\_\/\\\\\\\\\\\/____\//\\\____/\\\____/\\\__\/\\\_______\/\\\_     "
    print*, "    _\/\\\//_\//\\\____\/\\\_______\/\\\_\/\\\//////\\\_____\//\\\__/\\\\\__/\\\___\/\\\\\\\\\\\\\\\_    "
    print*, "     _\/\\\____\//\\\___\/\\\_______\/\\\_\/\\\____\//\\\_____\//\\\/\\\/\\\/\\\____\/\\\/////////\\\_   "
    print*, "      _\/\\\_____\//\\\__\//\\\______/\\\__\/\\\_____\//\\\_____\//\\\\\\//\\\\\_____\/\\\_______\/\\\_  "
    print*, "       _\/\\\______\//\\\__\///\\\\\\\\\/___\/\\\______\//\\\_____\//\\\__\//\\\______\/\\\_______\/\\\_ "
    print*, "        _\///________\///_____\/////////_____\///________\///_______\///____\///_______\///________\///__"
    stop
  end subroutine kurwout

end module ggint
