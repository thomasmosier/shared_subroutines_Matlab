!  spatiotemporal neighbors
!
!****************************************************************************
!
!  PROGRAM: spatiotemporal neighbors
!
!  PURPOSE: given n points ordered by increasing date, find the m nearest
!			neighbors to each observation from the past observations
!
!****************************************************************************

!	program spatiotemporalnn

!	implicit none


!	contains

	Subroutine tnmex4(xcoord,ycoord,m,n,istart,nnmat)
	implicit none
	integer n,m,i,j
	integer, intent(in)::istart
	real(kind=8), intent(in):: xcoord(n),ycoord(n)
	real(kind=8) xyc2(n),d2(n)
	integer, intent(out)::nnmat(m,n)
	integer irngt(n)

	irngt=0

	nnmat=1 !if a many 1s shows up, something is wrong or it is before istart
	
!	We precompute this for use below
	xyc2=0.5*(xcoord*xcoord+ycoord*ycoord)

!	OPEN(UNIT=666,FILE='N.ASC')
!	WRITE(666,*) n
!	CLOSE(666)

!	OPEN(UNIT=666,FILE='m.ASC')
!	WRITE(666,*) m
!	CLOSE(666)

!	OPEN(UNIT=666,FILE='xcoord.ASC')
!	WRITE(666,*) xcoord
!	CLOSE(666)

!	OPEN(UNIT=666,FILE='ycoord.ASC')
!	WRITE(666,*) ycoord
!	CLOSE(666)

!	OPEN(UNIT=666,FILE='istart.ASC')
!	WRITE(666,*) istart
!	CLOSE(666)

!	OPEN(UNIT=666,FILE='NNMATin.ASC')
!	WRITE(666,*) nnmat
!	CLOSE(666)


	do i=istart,n
	




!	The interest centers in the correct ordinal distance or correct ranks of some function
!	of distance to the observation i. Clearly, squaring the distances preserves the ranks.
!	The equation below gives squared spatial distance of observations in the past to the observation.
!	Note, the omission of the terms xcoord(i)*xcoord(i) and ycoord(i)*ycoord(i)
!	do not affect the ranks of the other distances relative to observation i as these are scalars
!	which are the same for all observations. One can test this by uncommenting the alternative
!	statement of squared-distance. The use of the straightfoward squared-distance costs
!	2 substractions, 2 multiplications, and one addition. The one used here costs 2 subtractions and
!	2 multiplications and so saves a little time.

	d2(1:i)=dsqrt((xcoord(i)-xcoord(1:i))**2+(ycoord(i)-ycoord(1:i))**2)

!	d2(1:i)=xyc2(1:i)-xcoord(i)*xcoord(1:i)-ycoord(i)*ycoord(1:i)



!This call shows the effect on the times of sorting all i distances (an unconditional ranking)
!	call RNKPAR (d2(1:i), IRNGT, i)



!This call shows the effect of just performing the partial rankings of distances (just need m).
	call UNIPAR (real(d2(1:i)), IRNGT, m)




!Stores the m ranks in the nearest neighbor matrix
	nnmat(:,i)=irngt(1:m)

	end do

!	OPEN(UNIT=666,FILE='NNMATout.ASC')
!	WRITE(666,*) nnmat
!	CLOSE(666)

contains


Subroutine unipar (XDONT, IRNGT, NORD)
!  Ranks partially XDONT by IRNGT, up to order NORD at most,
!  removing duplicate entries
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  quickly as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values smaller than the pivot
!  (ILOWT), and the indices of values larger than the pivot
!  that we might still need later on (IHIGT). It iterates
!  until it can bring the number of values in ILOWT to
!  exactly NORD, and then uses an insertion sort to rank
!  this set, since it is supposedly small. At all times, the
!  NORD first values in ILOWT correspond to distinct values
!  of the input array.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (InOut) :: NORD
! __________________________________________________________
      Real    :: XPIV, XWRK, XWRK1, XMIN, XMAX, XPIV0
!
      Integer, Dimension (SIZE(XDONT)) :: ILOWT, IHIGT
      Integer :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
      Integer :: IDEB, JDEB, IMIL, IFIN, NWRK, ICRS, IDCR, ILOW
      Integer :: JLM2, JLM1, JHM2, JHM1
!
      NDON = SIZE (XDONT)
!
!    First loop is used to fill-in ILOWT, IHIGT at the same time
!
      If (NDON < 2) Then
         If (NORD >= 1) Then
            NORD = 1
            IRNGT (1) = 1
         End If
         Return
      End If
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of low values.
!
     Do ICRS = 2, NDON
        If (XDONT(ICRS) == XDONT(1)) Then
          Cycle
        Else If (XDONT(ICRS) < XDONT(1)) Then
           ILOWT (1) = ICRS
           IHIGT (1) = 1
        Else
           ILOWT (1) = 1
           IHIGT (1) = ICRS
        End If
        Exit
     End Do
!
      If (NDON <= ICRS) Then
         NORD = Min (NORD, 2)
         If (NORD >= 1) IRNGT (1) = ILOWT (1)
         If (NORD >= 2) IRNGT (2) = IHIGT (1)
         Return
      End If
!
      ICRS = ICRS + 1
      JHIG = 1
      If (XDONT(ICRS) < XDONT(IHIGT(1))) Then
         If (XDONT(ICRS) < XDONT(ILOWT(1))) Then
            JHIG = JHIG + 1
            IHIGT (JHIG) = IHIGT (1)
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = ICRS
         Else If (XDONT(ICRS) > XDONT(ILOWT(1))) Then
            JHIG = JHIG + 1
            IHIGT (JHIG) = IHIGT (1)
            IHIGT (1) = ICRS
         End If
      ElseIf (XDONT(ICRS) > XDONT(IHIGT(1))) Then
         JHIG = JHIG + 1
         IHIGT (JHIG) = ICRS
      End If
!
      If (NDON <= ICRS) Then
         NORD = Min (NORD, JHIG+1)
         If (NORD >= 1) IRNGT (1) = ILOWT (1)
         If (NORD >= 2) IRNGT (2) = IHIGT (1)
         If (NORD >= 3) IRNGT (3) = IHIGT (2)
         Return
      End If
!
      If (XDONT(NDON) < XDONT(IHIGT(1))) Then
         If (XDONT(NDON) < XDONT(ILOWT(1))) Then
            Do IDCR = JHIG, 1, -1
              IHIGT (IDCR+1) = IHIGT (IDCR)
            End Do
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = NDON
            JHIG = JHIG + 1
         ElseIf (XDONT(NDON) > XDONT(ILOWT(1))) Then
            Do IDCR = JHIG, 1, -1
              IHIGT (IDCR+1) = IHIGT (IDCR)
            End Do
            IHIGT (1) = NDON
            JHIG = JHIG + 1
         End If
      ElseIf (XDONT(NDON) > XDONT(IHIGT(1))) Then
         JHIG = JHIG + 1
         IHIGT (JHIG) = NDON
      End If
!
      If (NDON <= ICRS+1) Then
         NORD = Min (NORD, JHIG+1)
         If (NORD >= 1) IRNGT (1) = ILOWT (1)
         If (NORD >= 2) IRNGT (2) = IHIGT (1)
         If (NORD >= 3) IRNGT (3) = IHIGT (2)
         If (NORD >= 4) IRNGT (4) = IHIGT (3)
         Return
      End If
!
      JDEB = 0
      IDEB = JDEB + 1
      JLOW = IDEB
      XPIV = XDONT (ILOWT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                   (XDONT(IHIGT(3))-XDONT(ILOWT(IDEB)))
      If (XPIV >= XDONT(IHIGT(1))) Then
         XPIV = XDONT (ILOWT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                      (XDONT(IHIGT(2))-XDONT(ILOWT(IDEB)))
         If (XPIV >= XDONT(IHIGT(1))) &
             XPIV = XDONT (ILOWT(IDEB)) + REAL (2*NORD) / REAL (NDON+NORD) * &
                                          (XDONT(IHIGT(1))-XDONT(ILOWT(IDEB)))
      End If
      XPIV0 = XPIV
!
!  One puts values > pivot in the end and those <= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the IHIGT array as soon as we have more
!  than enough values in ILOWT, i.e. one more than
!  strictly necessary so as to be able to come out of the
!  case where JLOWT would be NORD distinct values followed
!  by values that are exclusively duplicates of these.
!
!
      If (XDONT(NDON) > XPIV) Then
         lowloop1: Do
            ICRS = ICRS + 1
            If (XDONT(ICRS) > XPIV) Then
               If (ICRS >= NDON) Exit
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            Else
               Do ILOW = 1, JLOW
                 If (XDONT(ICRS) == XDONT(ILOWT(ILOW))) Cycle lowloop1
               End Do
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               If (JLOW >= NORD) Exit
            End If
         End Do lowloop1
!
!  One restricts further processing because it is no use
!  to store more high values
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) <= XPIV) Then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               Else If (ICRS >= NDON) Then
                  Exit
               End If
            End Do
         End If
!
!
      Else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         lowloop2: Do ICRS = ICRS + 1, NDON - 1
            If (XDONT(ICRS) > XPIV) Then
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            Else
               Do ILOW = 1, JLOW
                 If (XDONT(ICRS) == XDONT (ILOWT(ILOW))) Cycle lowloop2
               End Do
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               If (JLOW >= NORD) Exit
            End If
         End Do lowloop2
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) <= XPIV) Then
                  If (ICRS >= NDON) Exit
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               End If
            End Do
         End If
      End If
!
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      Do
         if (JLOW == NORD) Exit
         If (JLM2 == JLOW .And. JHM2 == JHIG) Then
!
!   We are oscillating. Perturbate by bringing JLOW closer by one
!   to NORD
!
           If (NORD > JLOW) Then
                XMIN = XDONT (IHIGT(1))
                IHIG = 1
                Do ICRS = 2, JHIG
                   If (XDONT(IHIGT(ICRS)) < XMIN) Then
                      XMIN = XDONT (IHIGT(ICRS))
                      IHIG = ICRS
                   End If
                End Do
!
                JLOW = JLOW + 1
                ILOWT (JLOW) = IHIGT (IHIG)
                IHIG = 0
                Do ICRS = 1, JHIG
                   If (XDONT(IHIGT (ICRS)) /= XMIN) then
                      IHIG = IHIG + 1
                      IHIGT (IHIG ) = IHIGT (ICRS)
                   End If
                End Do
                JHIG = IHIG
             Else
                ILOW = ILOWT (JLOW)
                XMAX = XDONT (ILOW)
                Do ICRS = 1, JLOW
                   If (XDONT(ILOWT(ICRS)) > XMAX) Then
                      IWRK = ILOWT (ICRS)
                      XMAX = XDONT (IWRK)
                      ILOWT (ICRS) = ILOW
                      ILOW = IWRK
                   End If
                End Do
                JLOW = JLOW - 1
             End If
         End If
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
!
!   We try to bring the number of values in the low values set
!   closer to NORD. In order to make better pivot choices, we
!   decrease NORD if we already know that we don't have that
!   many distinct values as a whole.
!
         IF (JLOW+JHIG < NORD) NORD = JLOW+JHIG
         Select Case (NORD-JLOW)
! ______________________________
         Case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            Select Case (JHIG)
!
!   Not enough values in high part either (too many duplicates)
!
            Case (0)
               NORD = JLOW
!
            Case (1)
               JLOW = JLOW + 1
               ILOWT (JLOW) = IHIGT (1)
               NORD = JLOW
!
!   We make a special case when we have so few values in
!   the high values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            Case (2)
               If (XDONT(IHIGT(1)) <= XDONT(IHIGT(2))) Then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
               ElseIf (XDONT(IHIGT(1)) == XDONT(IHIGT(2))) Then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
                  NORD = JLOW
               Else
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
               End If
               Exit
!
            Case (3)
!
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (3)
               If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                  IHIGT (3) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  End If
               End If
               JHIG = 1
               JLOW = JLOW + 1
               ILOWT (JLOW) = IHIGT (1)
               JHIG = JHIG + 1
               IF (XDONT(IHIGT(JHIG)) /= XDONT(ILOWT(JLOW))) Then
                 JLOW = JLOW + 1
                 ILOWT (JLOW) = IHIGT (JHIG)
               End If
               JHIG = JHIG + 1
               IF (XDONT(IHIGT(JHIG)) /= XDONT(ILOWT(JLOW))) Then
                 JLOW = JLOW + 1
                 ILOWT (JLOW) = IHIGT (JHIG)
               End If
               NORD = Min (JLOW, NORD)
               Exit
!
            Case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JHIG
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (IFIN)
               If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                  IHIGT (IFIN) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  End If
               End If
!
               JDEB = JLOW
               NWRK = NORD - JLOW
               IWRK1 = IHIGT (1)
               XPIV = XDONT (IWRK1) + REAL (NWRK) / REAL (NORD+NWRK) * &
                                      (XDONT(IHIGT(IFIN))-XDONT(IWRK1))
!
!  One takes values <= pivot to ILOWT
!  Again, 2 parts, one where we take care of the remaining
!  high values because we might still need them, and the
!  other when we know that we will have more than enough
!  low values in the end.
!
               JHIG = 0
               lowloop3: Do ICRS = 1, IFIN
                  If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                     Do ILOW = 1, JLOW
                        If (XDONT(IHIGT(ICRS)) == XDONT (ILOWT(ILOW))) &
                            Cycle lowloop3
                     End Do
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                     If (JLOW > NORD) Exit
                  Else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                  End If
               End Do lowloop3
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                  End If
               End Do
           End Select
!
! ______________________________
!
         Case (1)
!
!  Only 1 value is missing in low part
!
            XMIN = XDONT (IHIGT(1))
            IHIG = 1
            Do ICRS = 2, JHIG
               If (XDONT(IHIGT(ICRS)) < XMIN) Then
                  XMIN = XDONT (IHIGT(ICRS))
                  IHIG = ICRS
               End If
            End Do
!
            JLOW = JLOW + 1
            ILOWT (JLOW) = IHIGT (IHIG)
            Exit
!
! ______________________________
!
         Case (0)
!
!  Low part is exactly what we want
!
            Exit
!
! ______________________________
!
         Case (-5:-1)
!
!  Only few values too many in low part
!
            IRNGT (1) = ILOWT (1)
            Do ICRS = 2, NORD
               IWRK = ILOWT (ICRS)
               XWRK = XDONT (IWRK)
               Do IDCR = ICRS - 1, 1, - 1
                  If (XWRK < XDONT(IRNGT(IDCR))) Then
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  Else
                     Exit
                  End If
               End Do
               IRNGT (IDCR+1) = IWRK
            End Do
!
            XWRK1 = XDONT (IRNGT(NORD))
            insert1: Do ICRS = NORD + 1, JLOW
               If (XDONT(ILOWT (ICRS)) < XWRK1) Then
                  XWRK = XDONT (ILOWT (ICRS))
                  Do ILOW = 1, NORD - 1
                     If (XWRK <= XDONT(IRNGT(ILOW))) Then
                        If (XWRK == XDONT(IRNGT(ILOW))) Cycle insert1
                        Exit
                     End If
                  End Do
                  Do IDCR = NORD - 1, ILOW, - 1
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  End Do
                  IRNGT (IDCR+1) = ILOWT (ICRS)
                  XWRK1 = XDONT (IRNGT(NORD))
               End If
            End Do insert1
!
            Return
!
! ______________________________
!
         Case (:-6)
!
! last case: too many values in low part
!
            IDEB = JDEB + 1
            IMIL = MIN ((JLOW+IDEB) / 2, NORD)
            IFIN = MIN (JLOW, NORD+1)
!
!  One chooses a pivot from 1st, last, and middle values
!
            If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(IDEB))) Then
               IWRK = ILOWT (IDEB)
               ILOWT (IDEB) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
            End If
            If (XDONT(ILOWT(IMIL)) > XDONT(ILOWT(IFIN))) Then
               IWRK = ILOWT (IFIN)
               ILOWT (IFIN) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
               If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(IDEB))) Then
                  IWRK = ILOWT (IDEB)
                  ILOWT (IDEB) = ILOWT (IMIL)
                  ILOWT (IMIL) = IWRK
               End If
            End If
            If (IFIN <= 3) Exit
!
            XPIV = XDONT (ILOWT(IDEB)) + REAL(NORD)/REAL(JLOW+NORD) * &
                                      (XDONT(ILOWT(IFIN))-XDONT(ILOWT(1)))
            If (JDEB > 0) Then
               If (XPIV <= XPIV0) &
                   XPIV = XPIV0 + REAL(2*NORD-JDEB)/REAL (JLOW+NORD) * &
                                  (XDONT(ILOWT(IFIN))-XPIV0)
            Else
               IDEB = 1
            End If
!
!  One takes values > XPIV to IHIGT
!  However, we do not process the first values if we have been
!  through the case when we did not have enough low values
!
            JHIG = 0
            IFIN = JLOW
            JLOW = JDEB
!
            If (XDONT(ILOWT(IFIN)) > XPIV) Then
               ICRS = JDEB
              lowloop4: Do
                 ICRS = ICRS + 1
                  If (XDONT(ILOWT(ICRS)) > XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                     If (ICRS >= IFIN) Exit
                  Else
                     XWRK1 = XDONT(ILOWT(ICRS))
                     Do ILOW = IDEB, JLOW
                        If (XWRK1 == XDONT(ILOWT(ILOW))) &
                            Cycle lowloop4
                     End Do
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     If (JLOW >= NORD) Exit
                  End If
               End Do lowloop4
!
               If (ICRS < IFIN) Then
                  Do
                     ICRS = ICRS + 1
                     If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ILOWT (ICRS)
                     Else
                        If (ICRS >= IFIN) Exit
                     End If
                  End Do
               End If
           Else
              lowloop5: Do ICRS = IDEB, IFIN
                  If (XDONT(ILOWT(ICRS)) > XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                  Else
                     XWRK1 = XDONT(ILOWT(ICRS))
                     Do ILOW = IDEB, JLOW
                        If (XWRK1 == XDONT(ILOWT(ILOW))) &
                            Cycle lowloop5
                     End Do
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     If (JLOW >= NORD) Exit
                  End If
               End Do lowloop5
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                  End If
               End Do
            End If
!
         End Select
! ______________________________
!
      End Do
!
!  Now, we only need to complete ranking of the 1:NORD set
!  Assuming NORD is small, we use a simple insertion sort
!
      IRNGT (1) = ILOWT (1)
      Do ICRS = 2, NORD
         IWRK = ILOWT (ICRS)
         XWRK = XDONT (IWRK)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK < XDONT(IRNGT(IDCR))) Then
               IRNGT (IDCR+1) = IRNGT (IDCR)
            Else
               Exit
            End If
         End Do
         IRNGT (IDCR+1) = IWRK
      End Do
     Return

end subroutine unipar




	end Subroutine tnmex4


!	end program spatiotemporalnn
