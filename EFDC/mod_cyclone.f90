! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
    ! Parametric cyclone models
    ! 2021-01-20, Nghiem Tien Lam
    module cyclone

    use global
    
    implicit none

    type CycloneTrackPoint
        real(RKD) :: Time  = 0.0        !< time (Julian days)
        real(RKD) :: Lng = 0.0          !< central longitude (degrees)
        real(RKD) :: Lat = 0.0          !< central latitude (degrees)
        real(RKD) :: Xc = 0.0           !< central x position (m)
        real(RKD) :: Yc = 0.0           !< central y position (m)
        real(RKD) :: Pc = 1010.0        !< central pressure (hPa)
        real(RKD) :: Vmax = 0.0         !< maximum wind speed (m/s)
        real(RKD) :: Rmw = 0.0          !< radius of maximum wind speed (km)
        real(RKD) :: Vfx = 0.0          !< moving speed of cyclone in x direction (m/s)
        real(RKD) :: Vfy = 0.0          !< moving speed of cyclone in y direction (m/s)
        real(RKD) :: B = 2.5            !< Holland's shape parameter (-)
    end type
    
    type CycloneTrack
        type(CycloneTrackPoint), allocatable, dimension(:) :: points
        integer :: num_points = 0       !> number of points in this track
        real(RKD) :: start_time         !> time of start simulation (days)
        real(RKD) :: end_time           !> time of end simulation (days)
    contains 

        procedure, pass(self) :: InterpCycloneTrack
    end type
    
    real(RKD) :: rho_a = 1.15           ! Air density (kg/m3)
    real(RKD) :: Pinf = 1010.           ! Ambient pressure (hPa)
    real(RKD) :: thetaMax = 65          ! between (45 - 135 deg.)
    real(RKD) :: TC_KM = 0.7            ! Boundary layer effect
    real(RKD) :: TC_KF = 0.5            ! Forward moving effect
    real(RKD) :: TC_CSWP = 0.88         ! Conversion factor from 1-min. to 10-min sustained wind speed
    real(RKD) :: N_RAMP = 9.0           ! Ramp up factor
    
    real(RKD) :: rad2deg, deg2rad, fcor, sgnW, Rmw, delP, Vfm, Vfa, coeffB, CONVRT2
    type(CycloneTrack), allocatable, dimension(:), target :: cyclone_tracks
    integer :: num_cyclone_tracks = 0   !> number of cyclone tracks
    integer :: ICYCLONE = 0             !> 0=no, 1=McConochie (2004), 2=Willoughby (2006)

    contains

    
    subroutine ReadCyclones
        use convertwgs84
        use fson
        use fson_value_m, only: fson_value_count, fson_value_get
    
        type(fson_value), pointer :: json_data, tracks, track, points, item
        type(CycloneTrackPoint), pointer :: pt(:)
        real(RKD) :: lon(1), lat(1), xutm(1), yutm(1), delT, delX, delY
        integer :: i, j, count
        logical :: file_exists
        
        !ICYCLONE = 0
        
        inquire(FILE="cyclones.jnp", EXIST=file_exists)
        if (.not. file_exists) return
        
        call utmpars
        
        write(*,'(A)') 'READING CYCLONE INPUT FILE IN CYCLONES.JNP' 
        json_data => fson_parse("cyclones.jnp")
        
        !call fson_get(json_data, "cyclone_model", ICYCLONE)
        call fson_get(json_data, "tracks", tracks)
        num_cyclone_tracks = fson_value_count(tracks)
        allocate(cyclone_tracks(num_cyclone_tracks))
        do i = 1, num_cyclone_tracks
            cyclone_tracks(i).start_time = -1.e30
            cyclone_tracks(i).end_time = +1.e30
            track => fson_value_get(tracks, i)
            call fson_get(track, "start", cyclone_tracks(i).start_time)
            call fson_get(track, "end", cyclone_tracks(i).end_time)
            call fson_get(track, "points", points)
            count = fson_value_count(points)
            cyclone_tracks(i).num_points = count
            allocate(cyclone_tracks(i).points(count))
            pt => cyclone_tracks(i).points
            do j = 1, count
                item => fson_value_get(points, j)
                call fson_get(item, "time", pt(j).Time)
                call fson_get(item, "lon", pt(j).Lng)
                call fson_get(item, "lat", pt(j).Lat)
                call fson_get(item, "pc", pt(j).Pc)
                call fson_get(item, "vmax", pt(j).Vmax)
                call fson_get(item, "rmw", pt(j).Rmw)
                call fson_get(item, "b", pt(j).B)
                lon(1) = pt(j).Lng
                lat(1) = pt(j).Lat
                call UTM_WGS84(lon,lat,xutm,yutm)
                pt(j).Xc = xutm(1)
                pt(j).Yc = yutm(1)
                pt(j).Vfx = 0.
                pt(j).Vfy = 0.
                if (j > 1) then
                    delT = 86400.*(pt(j).Time - pt(j-1).Time)
                    delX = pt(j).Xc - pt(j-1).Xc
                    delY = pt(j).Yc - pt(j-1).Yc
                    if (ABS(delT) >= 1.0e-9) then
                        pt(j).Vfx = delX/delT
                        pt(j).Vfy = delY/delT
                        if (j == 2) then
                            pt(j-1).Vfx = pt(j).Vfx
                            pt(j-1).Vfy = pt(j).Vfy
                        endif    
                    else
                        pt(j).Vfx = pt(j-1).Vfx
                        pt(j).Vfy = pt(j-1).Vfy
                    endif
                endif
                if(cyclone_tracks(i).start_time < pt(j).Time) cyclone_tracks(i).start_time = pt(j).Time
                if(cyclone_tracks(i).end_time > pt(j).Time) cyclone_tracks(i).end_time = pt(j).Time
            enddo
        enddo
    end subroutine
    
    integer function InterpCycloneTrack(self, time, point)
    
        implicit none
        
        class(CycloneTrack), intent(in) :: self
        real(kind = RKD), intent(in) :: time
        type(CycloneTrackPoint), intent(inout) :: point
        real(RKD) :: deltaT, factor
        integer :: i

        do i = 1, self.num_points - 1
            if ((time >= self.points(i).Time) .AND. &
                (time <= self.points(i+1).Time)) then
                deltaT = self.points(i+1).Time - self.points(i).Time
                if(ABS(deltaT) > 1.0e-9) then
                    factor = (time - self.points(i).Time)/deltaT
                else
                    factor = 0.5
                endif
                point.Time  = time
                point.Lng = self.points(i).Lng + (self.points(i+1).Lng - self.points(i).Lng)*factor
                point.Lat = self.points(i).Lat + (self.points(i+1).Lat - self.points(i).Lat)*factor
                point.Xc = self.points(i).Xc + (self.points(i+1).Xc - self.points(i).Xc)*factor
                point.Yc = self.points(i).Yc + (self.points(i+1).Yc - self.points(i).Yc)*factor
                point.Pc = self.points(i).Pc + (self.points(i+1).Pc - self.points(i).Pc)*factor
                point.Vmax = self.points(i).Vmax + (self.points(i+1).Vmax - self.points(i).Vmax)*factor
                point.Rmw = self.points(i).Rmw + (self.points(i+1).Rmw - self.points(i).Rmw)*factor
                point.Vfx = self.points(i).Vfx + (self.points(i+1).Vfx - self.points(i).Vfx)*factor
                point.Vfy = self.points(i).Vfy + (self.points(i+1).Vfy - self.points(i).Vfy)*factor
                point.B = self.points(i).B + (self.points(i+1).B - self.points(i).B)*factor
                InterpCycloneTrack = 1
                return
            endif
        enddo
        InterpCycloneTrack = 0
    end function InterpCycloneTrack
    
    subroutine CycloneFields(time)
        real(kind = RKD), intent(in) :: time
        type(CycloneTrackPoint) :: pt
        integer :: ret, i
        real(kind = RKD) :: omega
        
        if(ICYCLONE==0) return
        
        rad2deg = 180. / PI
        deg2rad = PI / 180.

        omega = 2. * PI / (23.934469444*3600.)

        ! *** EE7.3 CONVERT TO 2 METERS FOR ALL CALCULATIONS (0.003 IS OPEN GRASSLAND Z0)
        CONVRT2   = LOG(2.0/0.003)/LOG(10.0/0.003)
        
        IF(TC_KM <= 0.) TC_KM = 0.9
        IF(TC_KF <= 0.) TC_KM = 0.5
        IF(TC_CSWP == 0.) TC_CSWP = 0.88
        
        !write(*,*) 'CycloneFields at @',time
        do i=1, num_cyclone_tracks
            ret = InterpCycloneTrack(cyclone_tracks(i), time, pt)
            if (ret == 1) then
                fcor = 2. * omega * SIN(deg2rad*pt.Lat)
                sgnW = SIGN(1., pt.Lat)
                delP = Pinf - pt.Pc
                Vfm = SQRT(pt.Vfx**2 + pt.Vfy**2)
                Vfa = ATAN2(pt.Vfy, pt.Vfx)
                if (ICYCLONE == 1) then
                    call CycloneFieldHolland1980(pt)
                elseif (ICYCLONE == 2) then
                    call CycloneFieldHubbert1991(pt)
                elseif (ICYCLONE == 3) then
                    call CycloneFieldMcConochie2004(pt)
                elseif (ICYCLONE == 4) then
                    call CycloneFieldWilloughby2006(pt)
                endif
            endif
        enddo        
    end subroutine CycloneFields
    
    subroutine CycloneFieldHolland1980(pt)
        implicit none

        class(CycloneTrackPoint), intent(inout) :: pt
        real(RKD) :: fc, omega, factor
        real(RKD) :: x, y, r, Pa, Wx, Wy, DEL, Vm
        integer :: L,LF,LL,ND
                
        ! Subtract translational speed of storm from maximum wind speed
        ! to avoid distortion in the Holland curve fit.  Added back later
        Vm = pt.Vmax - Vfm
        
        ! Convert wind speed (10 m) to top of atmospheric boundary layer
        ! to avoid distortion in the Holland curve fit.  Added back later
        Vm = Vm / TC_KM

        if(delP < 1.) delP = 1.
        Rmw = pt.Rmw
        coeffB = pt.B 
        if(coeffB <= 0.) coeffB = MIN(MAX(rho_a*EXP(1.0)*Vm**2/delP, 1.0), 2.5)

        !$OMP DO PRIVATE(ND,LF,LL,L,X,Y)
        DO ND=1,NDM  
            LF=2+(ND-1)*LDM  
            LL=MIN(LF+LDM-1,LA)
            
            DO L=LF,LL 
                x = DLON(L) !XCOR(L,5)
                y = DLAT(L) !YCOR(L,5)
                call CyclonePointHolland1980(pt, x, y, r, Pa, Wx, Wy)
                call assimilate_cyclone(L, r, Rmw, Pa, CONVRT2 * WX, CONVRT2 * WY)
            ENDDO
        ENDDO   ! *** END OF DOMAIN LOOP
        !$OMP END DO        
    end subroutine 
    
    subroutine CyclonePointHolland1980(pt, x, y, r, Pa, Wx, Wy)
    
        implicit none

        class(CycloneTrackPoint), intent(inout) :: pt
        !integer, intent(in)  :: L
        real(RKD), intent(in)  :: x, y
        real(RKD), intent(out)  :: r, Pa, WX, WY
        real(RKD) :: dx, dy, rf, Vc2, Vg, V10
        real(RKD) :: phi, theta, gamma, delta

        dx = x - pt.Xc
        dy = y - pt.Yc
        r = SQRT(dx * dx + dy * dy) / 1000.
        theta = ATAN2(dy, dx)
        Pa = pt.Pc
        Vg = 0.
        if ( r > 0.) then
            delta = (Rmw / r) ** coeffB
            Pa = pt.Pc + delP * EXP(-delta)
        
            ! Calculate gradient winds at the point
            rf = 1000. * r * ABS(fcor) / 2.
            Vc2 = 100. * coeffB * delP / rho_a * delta * EXP(-delta)
            Vg = SQRT(Vc2 + rf * rf) - rf
        endif

        phi = sgnW*PI/2. + theta
        
        ! Convert wind velocity from top of atmospheric boundary layer
        ! (which is what the Holland curve fit produces) to wind velocity 
        ! at 10 m above the earth's surface
        ! Also convert from 1 minute averaged winds to 10 minute averaged winds
        V10 = TC_KM * TC_CSWP * Vg
        
        gamma = TC_KF * ABS(V10)/pt.Vmax

        ! The surface wind speed at a 10 m elevation
        Wx = V10 * COS(phi) + gamma*pt.Vfx
        Wy = V10 * SIN(phi) + gamma*pt.Vfy
    end subroutine CyclonePointHolland1980
    
    subroutine CycloneFieldHubbert1991(pt)
        implicit none

        class(CycloneTrackPoint), intent(inout) :: pt
        real(RKD) :: fc, omega
        real(RKD) :: x, y, r, Pa, Wx, Wy, DEL, Vm
        integer :: L,LF,LL,ND
                
        thetaMax = 70.
        Rmw = pt.Rmw
        coeffB = pt.B 
        if(coeffB <= 0.) coeffB = 1.5 + (980. - pt.Pc)/120.
        
        !$OMP DO PRIVATE(ND,LF,LL,L,X,Y)
        DO ND=1,NDM  
            LF=2+(ND-1)*LDM  
            LL=MIN(LF+LDM-1,LA)
            
            DO L=LF,LL 
                x = DLON(L) !XCOR(L,5)
                y = DLAT(L) !YCOR(L,5)
                call CyclonePointHubbert1991(pt, x, y, r, Pa, Wx, Wy)
                call assimilate_cyclone(L, r, Rmw, Pa, CONVRT2 * WX, CONVRT2 * WY)
            ENDDO
        ENDDO   ! *** END OF DOMAIN LOOP
        !$OMP END DO        
    end subroutine 
    
    subroutine CyclonePointHubbert1991(pt, x, y, r, Pa, Wx, Wy)
    
        implicit none

        class(CycloneTrackPoint), intent(inout) :: pt
        !integer, intent(in)  :: L
        real(RKD), intent(in)  :: x, y
        real(RKD), intent(out)  :: r, Pa, WX, WY
        real(RKD) :: dx, dy, rf, Vc2, Vg, W
        real(RKD) :: phi, theta, beta, delta, asym

        beta = 25.  ! inflow angle
        
        dx = x - pt.Xc
        dy = y - pt.Yc
        r = SQRT(dx * dx + dy * dy) / 1000.
        theta = rad2deg*ATAN2(dy, dx)
        Pa = pt.Pc
        Vg = 0.
        if ( r > 0.) then
            delta = (Rmw / r) ** coeffB
            Pa = pt.Pc + delP * EXP(-delta)
            
            ! Calculate gradient winds at the point
            rf = 1000. * r * ABS(fcor) / 2.
            Vc2 = 100. * coeffB * delP / rho_a * delta * EXP(-delta)
            Vg = SQRT(Vc2 + rf * rf) - rf
        endif
        
        phi = deg2rad*(sgnW*(90. + beta) + theta)

        asym = COS(deg2rad*(theta + sgnW*thetaMax)- Vfa)
        W = TC_KM * Vg + TC_KF * asym*Vfm
        
        ! convert to 1-munite sustained wind speed
        W = W / TC_CSWP 
                
        Wx = W * COS(phi)
        Wy = W * SIN(phi)
    end subroutine CyclonePointHubbert1991
    
    subroutine CycloneFieldMcConochie2004(pt)
        implicit none

        class(CycloneTrackPoint), intent(inout) :: pt
        real(RKD) :: x, y, r, Pa, Wx, Wy, DEL
        integer :: L,LF,LL,ND
                
        ! Willoughby & Rahn (2004)
        coeffB = pt.B 
        if(coeffB <= 0.) coeffB = 0.886 + 0.0177 * pt.Vmax - 0.0094 * pt.Lat
        Rmw = pt.Rmw
        if (Rmw <= 0.) Rmw = 46.4 * EXP(-0.0155 * pt.Vmax + 0.0169 * pt.Lat)

        !$OMP DO PRIVATE(ND,LF,LL,L,X,Y)
        DO ND=1,NDM  
            LF=2+(ND-1)*LDM  
            LL=MIN(LF+LDM-1,LA)
            
            DO L=LF,LL 
                x = DLON(L) !XCOR(L,5)
                y = DLAT(L) !YCOR(L,5)
                call CyclonePointMcConochie2004(pt, x, y, r, Pa, Wx, Wy)
                call assimilate_cyclone(L, r, Rmw, Pa, CONVRT2 * WX, CONVRT2 * WY)
            ENDDO
        ENDDO   ! *** END OF DOMAIN LOOP
        !$OMP END DO        
    end subroutine 
    
    subroutine CyclonePointMcConochie2004(pt, x, y, r, Pa, Wx, Wy)
    
        implicit none

        class(CycloneTrackPoint), intent(inout) :: pt
        real(RKD), intent(in)  :: x, y
        real(RKD), intent(out)  :: r, Pa, WX, WY
        real(RKD) :: dx, dy, rf, Vc2, Vg, Vgr, V10, Vabs
        real(RKD) :: theta, gamma, delta, phi, asym, beta

        dx = x - pt.Xc
        dy = y - pt.Yc
        r = SQRT(dx * dx + dy * dy) / 1000.
        theta = rad2deg * ATAN2(dy, dx)
        Pa = pt.Pc
        Vg = 0.
        if ( r > 0.) then
            delta = (Rmw / r) ** coeffB
            Pa = pt.Pc + delP * EXP(-delta)
        
            ! Calculate gradient winds at the point
            rf = 1000. * r * ABS(fcor) / 2.
            Vc2 = 100. * coeffB * delP / rho_a * delta * EXP(-delta)
            Vg = SQRT(Vc2 + rf * rf) - rf
        endif

        ! Inflow angle (Sobey, 1977)
        beta = 25.
        if (r < Rmw) then
            beta = 10. * r / Rmw
        elseif (r < 1.2 * Rmw) then
            beta = 10. + 75. * (r / Rmw - 1.)
        endif
                
        phi = deg2rad*(sgnW*(90. + beta) + theta)
        
        asym = 0.5 * (1. + COS(deg2rad*(theta + sgnW*thetaMax)- Vfa))
        if(Vg /= 0.) then
            asym =  asym * Vg / ABS(Vg)
        endif
        Vgr = Vg + TC_KF * asym * Vfm
        Vabs = ABS(Vgr)
        
        ! Boundary layer coeff. (Harper, 2001)        
        gamma = 0.66
        if (Vabs < 6.) then
            gamma = 0.81
        else if (Vabs < 19.5) then
            gamma = 0.81 - 2.93e-3 * (Vabs - 6.)
        else if (Vabs < 45.) then
            gamma = 0.77 - 4.31e-3 * (Vabs - 19.5)
        endif
        
        ! Convert to 1-minute sustained wind speed:
        V10 = gamma * Vgr / TC_CSWP

        ! The surface wind speed at a 10 m elevation
        Wx = V10 * COS(phi)
        Wy = V10 * SIN(phi)
    end subroutine CyclonePointMcConochie2004

    subroutine CycloneFieldWilloughby2006(pt)
    
        implicit none

        class(CycloneTrackPoint), intent(inout) :: pt
        real(RKD) :: X1, X2, R1, R2, fr, vmax_gl,vmax_ss, nn, A, eq3_rhs, xi, r1_r2
        real(RKD) :: x, y, r, Pa, Wx, Wy
        integer :: L, LF,LL,ND
        
        X2 = 25.
        coeffB = pt.B 
        if(coeffB <= 0.) coeffB = 0.886 + 0.0177 * pt.Vmax - 0.0094 * pt.Lat     
        
        ! Remove forward speed from maximum wind speed
        ! Vmax:	Maximum 10-m 1-min sustained wind for the tropical cyclone	(m/s)
        ! vmax_ss	(Vmax,?sym):	Maximum 10-m 1-min sustained wind for the tropical cyclone with motion asymmetry removed	(m/s)
        vmax_ss = pt.Vmax - 0.5 * Vfm
        if (vmax_ss < 0.) vmax_ss = 0.
        
        ! Converts maximum 10-m 1-minute symmetric sustained wind speed to gradient
        ! wind speed.The conversion factor depends on whether the storm is over land
        ! or water.    
        fr = 0.9

        ! wind_gl (V_G, m/s):	Gradient level winds at grid point
        vmax_gl = vmax_ss / fr
        ! Calculate radius of maximum winds
        Rmw = 46.4 * EXP(-0.0155 * vmax_gl + 0.0169 * pt.Lat)
        ! Calculates X1 using Eqn 10a (Willoughby et al. 2006).
        X1 = 317.1 - 2.026 * vmax_gl + 1.915 * pt.Lat;
        ! Calculates n using Eqn 10b (Willoughby et al. 2006).
        nn = 0.4067 + 0.0144 * vmax_gl - 0.0038 * pt.Lat;
        ! Calculates A using Eqn 10c (Willoughby et al. 2006).
        A = 0.0696 + 0.0049 * vmax_gl - 0.0064 * pt.Lat
        if (A < 0) A = 0.
        ! Calculate right-hand side of Eqn. 3
        eq3_rhs = (nn * ((1 - A) * X1 + 25 * A)) / (nn * ((1 - A) * X1 + 25 * A) + Rmw)
        xi = SolveNewtonRaphson(eq3_rhs)

        ! Calculate radius to start of transition region
        r1_r2 = 15.
        if( Rmw > 20. ) r1_r2 = 25.
        R1 = Rmw - xi * r1_r2     ! R1 (km):	Lower boundary of the transition zone for Willoughby model
        R2 = R1 + r1_r2           ! R2 (km):	Upper boundary of the transition zone for Willoughby model
                
        !$OMP DO PRIVATE(ND,LF,LL,L,X,Y)
        DO ND=1,NDM  
            LF=2+(ND-1)*LDM  
            LL=MIN(LF+LDM-1,LA)
            
            DO L=LF,LL 
                x = DLON(L) !XCOR(L,5)
                y = DLAT(L) !YCOR(L,5)                
                call CyclonePointWilloughby2006(pt, x, y, r, Pa, Wx, Wy, &
                    vmax_gl, R1, R2, X1, X2, A, nn)
                call assimilate_cyclone(L, r, Rmw, Pa, CONVRT2 * WX, CONVRT2 * WY)
            ENDDO
        ENDDO   ! *** END OF DOMAIN LOOP
        !$OMP END DO  
    end subroutine CycloneFieldWilloughby2006

    subroutine CyclonePointWilloughby2006(pt, x, y, r, Pa, Wx, Wy, &
        vmax_gl, R1, R2, X1, X2, A, nn)
        implicit none
        class(CycloneTrackPoint), intent(inout) :: pt
        !integer, intent(in)  :: L
        real(RKD), intent(in)  :: x, y
        real(RKD), intent(inout) :: r, Pa, Wx, Wy
        real(RKD), intent(inout) :: vmax_gl, R1, R2, X1, X2, A, nn
        real(RKD) :: dx, dy, theta, delta, beta, Ws
        real(RKD) :: Vi, Vo, Vg, xi, gwd, phi, fr, cf, ww
        
        dx = x - pt.Xc
        dy = y - pt.Yc
        r = SQRT(dx * dx + dy * dy) / 1000.
        theta = rad2deg * ATAN2(dy, dx)

        Pa = pt.Pc
        if ( r > 0.) then
            delta = (Rmw / r) ** coeffB
            Pa = pt.Pc + delP * EXP(-delta)
        endif            

        ! Calculate gradient winds at the point
        Vi = vmax_gl * (r / Rmw) ** nn                                              !Vi (m/s):	Azimuthal average winds inside R1
        Vo = vmax_gl * ((1. - A) * EXP((Rmw - r) / X1) + A * EXP((Rmw - r) / X2))   !V0 (deg.):	Azimuthal average winds outside R2
        ! Vg (VG(r), m/s):	Azimuthal average winds, varies by radius r
        Vg = 0.
        if (r < R1) then 
            Vg = Vi
        elseif (r > R2) then
            Vg = Vo
        else
            xi = (r - R1) / (R2 - R1)
            ww = WeightingDeriv(xi)
            Vg = Vi * (1. - ww) + Vo * ww
        endif
        if (Vg < 0.) Vg = 0.

        ! Bring back to surface level (surface wind reduction factor) (Knaff et al., 2011)
        ! fr:	Reduction factor for converting between surface and gradient winds
        fr = 0.9
        if (r >= 700.) then
            fr = 0.75
        elseif (r > 100.) then
            fr = 0.90 - 0.15*(r - 100.) / 600.
        endif
        
        Ws = Vg * fr

        ! Calculate inflow angle over water based on radius of location from storm
        ! center in comparison to radius of maximum winds (Phadke et al. 2003)
        beta = 25.
        if (r < Rmw) then
            beta = 10. + (1. + r / Rmw)
        elseif (r < 1.2 * Rmw) then
            beta = 20. + 25. * (r / Rmw - 1.)
        endif

        ! Surface wind direction: Add inflow angle to gradient wind direction
        phi = deg2rad*(sgnW*(90. + beta) + theta)

        ! Add back in storm forward motion component
        ! Calculate u- and v-components of surface wind speed
        Wx = Ws * COS(phi)
        Wy = Ws * SIN(phi)

        ! Add back in component from forward motion of the storm
        cf = (Rmw * r) / (Rmw * Rmw + r * r)

       ! Asymmetric surface sustained wind at grid point: Add tangential and
       ! forward speed components and calculate magnitude of this total wind
       Wx = Wx + cf * pt.Vfx 
       Wy = Wy + cf * pt.Vfy
    end subroutine CyclonePointWilloughby2006

    subroutine assimilate_cyclone(L, r, Rmw, Pa, Wx, Wy)
        implicit none
        integer :: L
        real(RKD), intent(in) :: r, Rmw, Pa, Wx, Wy
        real(RKD) :: factor = 0.
        
        if(N_RAMP > 0.) factor = r/(N_RAMP*rmw)
        factor = factor ** 4
        !if (r <= Rmw) factor = 0.
        WNDVELE(L) = factor * WNDVELE(L) + (1. - factor) * WX
        WNDVELN(L) = factor * WNDVELN(L) + (1. - factor) * WY
        PATMT(L) = factor * PATMT(L) + (1. - factor) * Pa
        ATMP(L) = PATMT(L)*0.0101974*G 
    end subroutine    

    real(RKD) function WeightingFunc(xi)
        implicit none
        real(RKD), intent(in) :: xi
        WeightingFunc = ((((((((70. * xi - 315.) * xi + 540.) * xi - 420.) * xi + 126.) * xi) * xi) * xi) * xi) * xi
    end function

    real(RKD) function WeightingDeriv(xi)
        implicit none
        real(RKD), intent(in) :: xi
        WeightingDeriv = (((((((630. * xi - 2520.) * xi + 3780.) * xi - 2520.) * xi + 630.) * xi) * xi) * xi) * xi
    end function
    
    real(RKD) function SolveNewtonRaphson(rhs)
        implicit none
        real(RKD), intent(in) :: rhs
        real(RKD) :: x0 = 0.5, eps = 10e-4
        real(RKD) :: fx, df, dx, x
        integer :: i = 0, itmax = 100
        
        x = x0
        do i = 1, itmax
            df = WeightingDeriv(x)
            fx = WeightingFunc(x) - rhs
            if (Abs(fx) <= eps) exit
            dx = -fx / df
            x = x + dx
        enddo
        SolveNewtonRaphson = x
    end function

    end module
