MODULE pf
  USE com
  USE l63
  IMPLICIT NONE
  PUBLIC

  !------------particle members-------------
  !
  ! x_da_k(t,:) = x_da(t) + x_prtb(t,:)
  !
  !-----------------------------------------
  REAL(dp) :: x_da_k(0:nt_asm,m)
  REAL(dp) :: y_da_k(0:nt_asm,m)
  REAL(dp) :: z_da_k(0:nt_asm,m)

  !-------------particle weights------------
  REAL(dp) :: W(0:nt_asm/obs_interval,m)
  !-----------------------------------------
  ! LIKELIHOOD OF EACH PARTICLES
  !-----------------------------------------
!  REAL(dp) :: L(0:nt_asm/obs_interval,m)
!  REAL(dp) :: Neff(0:nt_asm/obs_interval)
!  REAL(dp) :: N_thres=REAL(m/2)

CONTAINS
  SUBROUTINE sir(m,p,Y,Win,H,R,x_da_k,y_da_k,z_da_k,x_da,y_da,z_da,Wout)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: m!particles
    INTEGER,INTENT(IN) :: p !observation dimension
    REAL(dp),INTENT(IN) :: Y(p)
    REAL(dp),INTENT(IN) :: Win(m)
    REAL(dp),INTENT(IN) :: H(p,n)
    REAL(dp),INTENT(IN) :: R(p,p)
    REAL(dp),INTENT(INOUT) :: x_da_k(m),y_da_k(m),z_da_k(m)
    REAL(dp),INTENT(OUT) :: x_da,y_da,z_da
    REAL(dp),INTENT(OUT) :: Wout(m)
    REAL(dp) :: Wtmp(m)
    REAL(dp) :: L(m)
    REAL(dp) :: Neff
    REAL(dp) :: N_thres
    REAL(dp) :: Pf(n,n)
    INTEGER :: k,k9
    REAL(dp) :: X_k(n)
    REAL(dp) :: innov(p)
    REAL(dp) :: innov_k
    REAL(dp) :: W_sq(m)
    REAL(dp) :: noise1,noise2
    REAL(dp) :: gnoise!,gnoise(3,1)
    REAL(dp) :: Wcum(0:m) ! cummulative weight

    REAL(dp) :: X_mean(1:n)
    REAL(dp) :: X_prtb(1:n,1:m)

    X_mean(1)=SUM(x_da_k(1:m))/m
    X_mean(2)=SUM(y_da_k(1:m))/m
    X_mean(3)=SUM(z_da_k(1:m))/m

    DO k=1,m
      X_prtb(1,k)=x_da_k(k)-X_mean(1)
      X_prtb(2,k)=y_da_k(k)-X_mean(2)
      X_prtb(3,k)=z_da_k(k)-X_mean(3)
    ENDDO

    Pf=1.0d0/(m-1.0d0)*MATMUL(X_prtb,TRANSPOSE(X_prtb))

    N_thres=REAL(m/2)

    Wtmp(1:m)=Win(1:m)

    WRITE(*,*)
    WRITE(*,*) '--------Weights--------'
    WRITE(*,'(100F10.7)') Wtmp(1:m)

    Wcum(0)=0.0d0
    DO k=1,m
      Wcum(k)=Wcum(k-1)+Wtmp(k)
    ENDDO

    WRITE(*,*)
    WRITE(*,*) '-------Cummulative Weight-------'
    WRITE(*,'(100F10.7)') Wcum(1:m)

    DO k=1,m
      X_k(1)=x_da_k(k)
      X_k(2)=y_da_k(k)
      X_k(3)=z_da_k(k)
      innov(1:p)=Y(1:p)-MATMUL(H,X_k)
      innov_k=SUM(innov**2)/n
      L(k)=EXP(-(innov_k/(2*(R(1,1)+R(2,2)+R(3,3))**2)))
    ENDDO

    WRITE(*,*)
    WRITE(*,*) '---------Likelihoods--------'
    WRITE(*,'(100F10.7)') L(1:m)

    DO k=1,m
      W_sq(k)=Wtmp(k)**2
    ENDDO

    Neff=1.0d0/SUM(W_sq(1:m))
    WRITE(*,*)
    WRITE(*,*) '---------Effective sample size--------'
    WRITE(*,'(F10.5)') Neff

    IF (Neff < N_thres) THEN
      DO k=1,m
        !--------------------------------------------------
        ! SIR algorithm
        ! sample x_da_k(t,k9) etc. with probability W
        !--------------------------------------------------
        CALL RANDOM_NUMBER(noise1)

        k9=0
        DO WHILE (Wcum(k9) < noise1)
          k9=k9+1
        ENDDO

        WRITE(*,*)
        WRITE(*,*)'----------resampled particle----------'
        WRITE(*,'(I5)') k9
        CALL RANDOM_NUMBER(noise1)
        CALL RANDOM_NUMBER(noise2)
        gnoise=SQRT(0.01d0*Pf(1,1))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
          &*COS(2.0d0*PI*noise2)
        x_da_k(k)=x_da_k(k9)+gnoise
        CALL RANDOM_NUMBER(noise1)
        CALL RANDOM_NUMBER(noise2)
        gnoise=SQRT(0.01d0*Pf(2,2))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
          &*COS(2.0d0*PI*noise2)
        y_da_k(k)=y_da_k(k9)+gnoise
        CALL RANDOM_NUMBER(noise1)
        CALL RANDOM_NUMBER(noise2)
        gnoise=SQRT(0.01d0*Pf(3,3))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
          &*COS(2.0d0*PI*noise2)
        z_da_k(k)=z_da_k(k9)+gnoise
      ENDDO
      Wtmp(1:m)=1.0d0/m
    ENDIF

    DO k=1,m
      Wout(k)=L(k)*Wtmp(k)
    ENDDO
    Wout(1:m)=Wout(1:m)/SUM(Wout(1:m))
    WRITE(*,*)
    WRITE(*,*) '--------updated Weight--------'
    WRITE(*,'(100F10.7)') Wout(1:m)

    x_da=DOT_PRODUCT(x_da_k(1:m),Wout(1:m))
    y_da=DOT_PRODUCT(y_da_k(1:m),Wout(1:m))
    z_da=DOT_PRODUCT(z_da_k(1:m),Wout(1:m))
  END SUBROUTINE sir

  SUBROUTINE dbf(H,x_da,y_da,z_da,Y,Prin,Prout)
    IMPLICIT NONE

    REAL(dp),PARAMETER :: q=5
    REAL(dp),INTENT(IN) :: H(p,n)
    REAL(dp),INTENT(IN) :: x_da,y_da,z_da
    REAL(dp),INTENT(IN) :: Y(p)
    REAL(dp),INTENT(IN) :: Prin(q)
    REAL(dp),INTENT(OUT) :: Prout(q)
    REAL(dp),INTENT(OUT) :: m(q)

    REAL(dp) :: W(q)
    !-------------------------
    ! discrete bayes filter
    !-------------------------
    DO s=1,q
      Hx(:,s)=0.0d0
      DO i=1,p
        Hx(i,s)=Hx(i,s)+H(i,1)*x_da(t,s)+H(i,2)*y_da(t,s)+H(i,3)*z_da(t,s)
      ENDDO
      W(t/obs_interval,s)=1.0d0/SQRT(SUM((Hx(:,s)-Y(t/obs_interval,:))**2)/p)
    ENDDO
    WRITE(*,*)
    WRITE(*,*) "---------Likelihood---------"
    WRITE(*,'(5F25.12)') W(t/obs_interval,:)
    DO s=1,q
      Pr(t/obs_interval,s)=W(t/obs_interval,s)*Pr(t/obs_interval-1,s)
    ENDDO
      Pr(t/obs_interval,:)=Pr(t/obs_interval,:)/SUM(Pr(t/obs_interval,1:q))
    WRITE(*,*)
    WRITE(*,*) "---------Probability---------"
    WRITE(*,'(5F25.12)') Pr(t/obs_interval,:)

    m(t,:)=0
    DO k4=1,Nall
      CALL RANDOM_NUMBER(noise1)
      cpr=0.0d0
      s=0
      DO WHILE (cpr .LE. noise1)
        s=s+1
        cpr=cpr+Pr(t/obs_interval,s)
      ENDDO
      m(t,s)=m(t,s)+1
    ENDDO
    WRITE(*,*)
    WRITE(*,*) "-------optimal ensemble sizes--------"
    WRITE(*,'(5I7)') m(t,:)

!    DEALLOCATE(x_da_k,y_da_k,z_da_k)
!    ALLOCATE(x_da_k(Nall,q),y_da_k(Nall,q),z_da_k(Nall,q))
!    DO s=1,q
!      DO k=1,m(t,s)
!        CALL RANDOM_NUMBER(noise1)
!        CALL RANDOM_NUMBER(noise2)
!        gnoise=SQRT(Pf(1,1))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
!          &*COS(2.0d0*PI*noise2)
!        x_da_k(k,s)=x_da(t,s)+gnoise
!
!        CALL RANDOM_NUMBER(noise1)
!        CALL RANDOM_NUMBER(noise2)
!        gnoise=SQRT(Pf(2,2))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
!          &*COS(2.0d0*PI*noise2)
!        y_da_k(k,s)=y_da(t,s)+gnoise
!
!        CALL RANDOM_NUMBER(noise1)
!        CALL RANDOM_NUMBER(noise2)
!        gnoise=SQRT(Pf(3,3))*SQRT(-2.0d0*LOG(1.0d0-noise1))&
!          &*COS(2.0d0*PI*noise2)
!        z_da_k(k,s)=z_da(t,s)+gnoise
!      ENDDO
!    ENDDO
  END SUBROUTINE dbf
END MODULE pf
