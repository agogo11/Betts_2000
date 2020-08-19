PRO runmlm_betts2000

  Path = 'G:\MyIDLpackage\MLM_Betts2000\'
  figname = 'Fig2'

  ;; Initialize Parameters

  ; Reference parameters from Table 1 OF B00
  ; name      value         units         description
  ga_0    =   0.025;      % m/s           aerodynamic conductance
  gv_0    =   1./60;       % m/s           vegetative conductance to water vapor
  kent    =   0.2;        % -             entrainment parameter relating ML-top theta-flux to sfc. buoyancy flux
  Gamma_0 =   6e-4;       % K/Pa          stability |d \theta/d p| above ML top
  Qs_0    =   150.;        % W/m^2         surface net radiation
  Qbrad0  =   3.;          % K/day         ML radiative cooling rate
  Qbevap0 =   0.;          % K/day         ML cooling rate by evaporating rain

  ; renormalize Qbrad AND Qbevap to K/s
  Qbrad0   =   Qbrad0/86400.;
  Qbevap0  =   Qbevap0/86400.;

  ; Other reference parameters (from equations/text OF B00, and new to this code)
  ; name      value         units         description
  Cp      =   1006.;       % J/kg/K        specific heat capacity of dry air (B13-pc)
  Rd      =   287.;        % J/kg/K        gas constant for dry air
  RdCp    =   0.286;      % -             Ratio of Rd/Cp
  Lv      =   2.501e6;    % J/kg          latent heat of vaporization of water at 0 C (B13-pc)
  g       =   9.81;       % m/s^2         gravitational constant
  Th_00   =   303.;        % K             constant potential temperature offset in B00 eqn. (12)
  PH_00   =   60e2;       % Pa            reference P_H offset in B00 eqn. (12)
  P_T0    =   100e2;      % Pa            reference subsaturation just above ML top (B00 eqn. (13))
  eps     =   0.622;      % -             ratio of gas constants / MWs for (water vapor):(dry air)
  Cs      =   1e6;        % J/m^2/K       surface heat capacity (not in B00)
  ps      =   940e2;      % Pa            surface pressure
  p0      =   1e5;        % Pa            reference pressure for definition of theta
  nt      =   15000.;      % -             number of timesteps
  dt      =   200.;        % s             timestep
  afilt   =   0.1;       % -             Asselin filter parameter

  ; Create simulation-specific arrays
  ; (for gamma, P_T, Qs, Qbrad, Qbevap; gv, Qs, ga)

  gv_arr = 1./[60,80,100,150,300,500,700,900]
  nsims = N_ELEMENTS(gv_arr)

  Qs_arr = MAKE_ARRAY(nsims, value = Qs_0,/float)
  ga_arr = MAKE_ARRAY(nsims, value = ga_0,/float)

  CASE figname OF
    'Fig1a': BEGIN
      Gam_arr     =   [4., 5., 6., 7.]*(1e-4);  % K/Pa      --stability array
      P_T_arr     =   [1, 1, 1, 1]*P_T0;     % Pa        --subsaturation array
      Qbrad_arr   =   [1, 1, 1, 1]*Qbrad0;  % K/s       --radiative cooling array
      Qbevap_arr  =   [1, 1, 1, 1]*Qbevap0; % K/s       --evaporative cooling array
      nsens      = 4
    END

    'Fig1b': BEGIN
      Gam_arr     =   [1., 1., 1.]*Gamma_0;  % K/Pa      --stability array
      P_T_arr     =   [60., 100., 140.]*(1e2);     % Pa        --subsaturation array
      Qbrad_arr   =   [1, 1, 1, 1]*Qbrad0;  % K/s       --radiative cooling array
      Qbevap_arr  =   [1, 1, 1, 1]*Qbevap0; % K/s       --evaporative cooling array
      nsens      = 3
    END

    'Fig2': BEGIN
      Gam_arr     =   [1.]*Gamma_0;  % K/Pa      --stability array
      P_T_arr     =   [1.]*P_T0;     % Pa        --subsaturation array
      Qbrad_arr   =   [1.]*Qbrad0;  % K/s       --radiative cooling array
      Qbevap_arr  =   [1.]*Qbevap0; % K/s       --evaporative cooling array
      nsens      = 1
    END
  ENDCASE


  ; initialize plot/output variables to zero FOR new simulation
  P_LCL       = FLTARR(nsens,nsims);
  thetab_out  = FLTARR(nsens,nsims);
  qb_out      = FLTARR(nsens,nsims);
  Tb_out      = FLTARR(nsens,nsims);
  Ts_out      = FLTARR(nsens,nsims);
  LHF_out     = FLTARR(nsens,nsims);
  SHF_out     = FLTARR(nsens,nsims);
  delth_out   = FLTARR(nsens,nsims);
  delq_out    = FLTARR(nsens,nsims);
  omega_out   = FLTARR(nsens,nsims);

  ;; Loop over values OF gamma, P_T, Qbrad, or Qbevap
  FOR isens=0, nsens - 1 DO BEGIN
    ; set control parameters to their appropriate values
    Gamma   =   Gam_arr[isens];
    P_T     =   P_T_arr[isens];
    Qbrad   =   Qbrad_arr[isens];
    Qbevap  =   Qbevap_arr[isens];

    ;; Loop over single factor values
    FOR ivar=0, nsims - 1 DO BEGIN
      gv = gv_arr[ivar];
      Qs = Qs_arr[ivar];
      ga = ga_arr[ivar];

      ;; Variable initialization

      ; Prognostic variables: 3-element arrays are FOR Leapfrog/Asselin timestepping
      ;
      Ts          =   (FLTARR(3) +Th_00)*(ps/p0)^(RdCp);   % surface temperature
      thetab      =   FLTARR(3) +Th_00-1.;                     % ML potential temperature
      Tb          =   thetab[1]*(ps/p0)^(RdCp);            % surface air temperature
      qstar_b     =   mixr_sat(Tb,ps/100.)/1000.;                   % svp of surface air
      RH_b        =   0.7;
      qb          =   FLTARR(3) + RH_b*qstar_b;                % ML mixing ratio -- initialize at 70% RH

      ; prognostic variable tendencies
      Ts_dot      =   0.;
      thetab_dot  =   0.;
      qb_dot      =   0.;

      ; other variables that are used in diagnostic equations only
      ; are declared within the relevant diagnostic equation block

      ;; Integration in time
      FOR it= 0, nt - 1 DO BEGIN
        ; Diagnostic equations (D #*)
        ; (D #0) Diagnose surface air density using ideal gas law
        ; (B13-pc)
        Tb = thetab[1]*(ps/p0)^(RdCp);
        rho0 = ps/(Rd*Tb*(1.+0.608*qb[1]));
        ; ---------- END (D #0) ----------

        ; (D #1) Diagnose P_H based on thetab AND qb using B00 Eqns. (14) AND (22)
        qstar_b = mixr_sat(Tb,ps/100.)/1000.;
        RH_b = qb[1]/qstar_b;
        ; 3 following lines are from B02
        Tsat_b = 55.+2840./(3.5*ALOG(thetab[1])-ALOG(1000.*qb[1]/(0.622+qb[1]))-4.805);
        psat_b = p0*(Tsat_b/thetab[1])^3.4965;
        P_H = ps-psat_b;
        ;A = eps*Lv/(2*Cp*Tb);
        ;P_H = ps*(1-RH_b)/(A+(A-1)*RH_b); % N.B. units of P_H are Pascals

        ; check that P_H is +ve AND does NOT exceed 0.9*p0
        IF P_H LT 0 THEN BEGIN
          P_H = PH_00;
        ENDIF ELSE BEGIN
          IF P_H GT ps THEN P_H = 0.9*ps;
        ENDELSE
        ; ---------- END (D #1) ----------

        ; (D #2) Diagnose thetat based on B00 Eqn. (12)
        thetat = Th_00 + Gamma*(P_H - PH_00);
        ; check FOR static stability at ML top
        IF thetat LT thetab[1] THEN thetat = thetab[1]+0.01;
        ; ---------- END (D #2) ----------

        ; (D #3) Diagnose qt based on B00 Eqn. (14)
        Tt = thetat*((ps-P_H)/p0)^(RdCp);
        qstar_t = mixr_sat(Tt,(ps-P_H)/100.)/1000.;
        A = eps*Lv/(2.*Cp*Tt);
        RH_t = (ps - P_H - A*P_T)/(ps - P_H + (A-1.)*P_T);
        qt = RH_t*qstar_t/(1.+(qstar_t/0.622)*(1.-RH_t));
        ; check that qt is +ve
        IF (qt LT 0) THEN qt=0.;

        ; ---------- END (D #3) ----------

        ; (D #4) Diagnose Omega_t based on ML-top closure FOR SH flux, B00 Eqns.
        ; (8) AND (15)
        SHF = Cp*rho0*ga*(Ts[1]-Tb);
        qstar_s = mixr_sat(Ts[1],ps/100.)/1000.;
        LHF = Lv*rho0*ga*gv/(ga+gv)*(qstar_s-qb[1]);
        ;kFthv = 0.608*Cp*Tb/Lv;
        kFthv = 0.073;

        F0_thv = (SHF+kFthv*LHF)/Cp;
        numer = g*kent*Cp*F0_thv;
        denom = (Cp*(thetat-thetab[1])+kFthv*Lv*(qt-qb[1]));

        ; Check that ML-top omega is NOT diagnosed to be negative FOR
        ; one OR more reasons (-ve sfc. buoyancy flux OR -ve buoyancy
        ; jump at ML top)
        IF (numer LT 0) OR (denom LE 0) THEN BEGIN
          Omega_t = 0.;
        ENDIF ELSE BEGIN
          Omega_t = numer/denom;
        ENDELSE
        ; ---------- END (D #4) ----------

        ;
        ; Calculate tendencies FOR prognostic equations (P #*)
        ;

        ; (P #1) surface energy balance
        Ts_dot = 1./Cs*(Qs - SHF - LHF);
        ; ---------- END (P #1) ----------

        ; (P #2) ML thermal balance
        thetab_dot = (g/(Cp*P_H))*(SHF + (Cp*Omega_t/g)*(thetat-thetab[1]) $
          - (Cp*P_H/g)*(Qbrad+Qbevap));
        ; ---------- END (P #2) ----------

        ; (P #3) ML moisture balance
        qb_dot = (g/P_H)*(LHF/Lv + (Omega_t/g)*(qt-qb[1]) + (Cp*P_H)/(Lv*g)*Qbevap);
        ; ---------- END (P #3) ----------
        ;
        ;
        ;        PRINT, Ts[1],thetab[1],qb[1], p_h/100.
        ;
        ;        IF thetab[1] GT 308. THEN BEGIN
        ;          PRINT, ' '
        ;        ENDIF

        ;
        ; Advance prognostic variables in time
        ;

        Ts[2] = Ts[0]+2.*dt*Ts_dot;
        thetab[2] = thetab[0]+2.*dt*thetab_dot;
        qb[2] = qb[0]+2.*dt*qb_dot;

        Ts[0] = Ts[1]+afilt*(Ts[0]-2.*Ts[1]+Ts[2]);
        thetab[0] = thetab[1]+afilt*(thetab[0]-2.*thetab[1]+thetab[2]);
        qb[0] = qb[1]+afilt*(qb[0]-2*qb[1]+qb[2]);

        Ts[1] = Ts[2];
        thetab[1] = thetab[2];
        qb[1] = qb[2];

      ENDFOR

      ; set P_LCL AND other variables to the equilibrated values
      P_LCL[isens, ivar]         =   P_H/100.;
      thetab_out[isens, ivar]    =   thetab[1];
      qb_out[isens, ivar]        =   qb[1];
      Tb_out[isens, ivar]        =   Tb;
      Ts_out[isens, ivar]       =   Ts[1];
      LHF_out[isens, ivar]      =   LHF;
      SHF_out[isens, ivar]       =   SHF;
      delth_out[isens, ivar]     =   thetat-thetab[1];
      delq_out[isens, ivar]      =   1000.*(qt-qb[1]);
      omega_out[isens, ivar]     =   Omega_t;

      IF (gv EQ gv_0) AND (Qs EQ Qs_0) AND (ga EQ ga_0) THEN BEGIN
        iref = ivar;
        E0 = LHF;
        H0 = SHF;
      ENDIF

      ;PRINT, ' '
    ENDFOR

  ENDFOR

  ;plot figure
  CASE figname OF
    'Fig1a': BEGIN
      cgps_open, Path + Figname + '_Betts2000.eps',/CMYK,$
        /nomatch,Font = 1,/Times, xsize =8/2.54, ysize = 8/2.54, fontsize =14

      pos = cglayout([1,1], OXMargin=[4.5,1.], OYMargin=[4.,1.])

      mycharsize = 1.
      myxrange = [60., 900.]
      myyrange = [50., 260.]

      myxtitle = 'R_v (s m-1)'
      myytitle = 'P_LCL (hPa)'

      mycolors = ['red1','red3','red5','red7']

      location_lg = [600, 110]
      text_lg = ['0.04','0.05','0.06','0.07']
      linestyle_lg = MAKE_ARRAY(nsens, value = 0)

      ;Fig 1
      cgplot,1./gv_arr, P_LCL[0,*],/noerase,/nodata, $
        xtitle = myxtitle, ytitle = myytitle,$
        xrange = myxrange,yrange = myyrange,$
        xtickinterval = 200., ytickinterval = 50,$
        charsize = mycharsize, POSITION = pos[*,0]

      FOR isens = 0, nsens - 1 DO BEGIN
        cgoplot,1./gv_arr, P_LCL[isens,*], thick = 2., color = mycolors[isens]
        ;cgtext, xloc_lg[isens], yloc_lg[isens], text_lg[isens], color = mycolors[isens], charsize = 0.7*mycharsize
      ENDFOR

      cglegend, title = text_lg, linestyle = linestyle_lg, color = mycolors, Location = location_lg, /data,$
        vspace = 0.8, charsize = 0.7*mycharsize, tcolors = mycolors

      cgtext, location_lg[0], location_lg[1] + 8., 'Stability (K hPa-1)', color = 'red7', charsize = 0.7*mycharsize
      cgps_close, /png, density = 1000
    END
    'Fig1b': BEGIN
      cgps_open, Path + Figname + '_Betts2000.eps',/CMYK,$
        /nomatch,Font = 1,/Times, xsize =8/2.54, ysize = 8/2.54, fontsize =14

      pos = cglayout([1,1], OXMargin=[4.5,1.], OYMargin=[4.,1.])

      mycharsize = 1.
      myxrange = [60., 900.]
      myyrange = [50., 250.]

      myxtitle = 'R_v (s m-1)'
      myytitle = 'P_LCL (hPa)'

      mycolors = ['red3','red5','red7']

      location_lg = [600, 110]
      text_lg = ['60','100','140']
      linestyle_lg = MAKE_ARRAY(nsens, value = 0)

      ;Fig 1
      cgplot,1./gv_arr, P_LCL[0,*],/noerase,/nodata, $
        xtitle = myxtitle, ytitle = myytitle,$
        xrange = myxrange,yrange = myyrange,$
        xtickinterval = 200., ytickinterval = 50,$
        charsize = mycharsize, POSITION = pos[*,0]

      FOR isens = 0, nsens - 1 DO BEGIN
        cgoplot,1./gv_arr, P_LCL[isens,*], thick = 2., color = mycolors[isens]
        ;cgtext, xloc_lg[isens], yloc_lg[isens], text_lg[isens], color = mycolors[isens], charsize = 0.7*mycharsize
      ENDFOR

      cglegend, title = text_lg, linestyle = linestyle_lg, color = mycolors, Location = location_lg, /data,$
        vspace = 0.8, charsize = 0.7*mycharsize, tcolors = mycolors

      cgtext, location_lg[0], location_lg[1] + 8., 'ML top P_T (hPa)', color = 'red7', charsize = 0.7*mycharsize
      cgps_close, /png, density = 1000
    END
    'Fig2': BEGIN
      cgps_open, Path + Figname + '_Betts2000.eps',/CMYK,$
        /nomatch,Font = 1,/Times, xsize =9/2.54, ysize = 8/2.54, fontsize =14

      pos = cglayout([1,1], OXMargin=[4.5,4.5], OYMargin=[4.,1.])

      mycharsize = 1.
      myxrange = [50., 260.]
      myyrange = [0., 150.]

      myxtitle = 'P_LCL (hPa)'
      myytitle = 'SH,LH (Wm-2)'

      mycolors = ['blu7','red7','black']

      location_lg = [175, 0.95]
      text_lg = ['LH','SH','EF']
      linestyle_lg = MAKE_ARRAY(nsens, value = 0)

      ;Fig 1
      cgplot,1./gv_arr, P_LCL[0,*],/noerase,/nodata, $
        xtitle = myxtitle, ytitle = myytitle,$
        xrange = myxrange,yrange = myyrange,$
        xtickinterval = 50., ytickinterval = 50,$
        charsize = mycharsize, POSITION = pos[*,0], YStyle=8

      cgoplot, P_LCL[0,*],LHF_out[0,*], thick = 4., color = mycolors[0]
      cgoplot, P_LCL[0,*],SHF_out[0,*], thick = 4., color = mycolors[1]

      cgaxis, YAxis=1, YRange=[0.4, 1.], title = 'EF',/Save, charsize = mycharsize
      cgoplot, P_LCL[0,*],LHF_out[0,*]/(LHF_out[0,*] + SHF_out[0,*]), thick = 4., color = mycolors[2]

      cglegend, title = text_lg, linestyle = linestyle_lg, color = mycolors, Location = location_lg, /data,$
        vspace = 0.8, charsize = 0.7*mycharsize, tcolors = mycolors

      cgps_close, /png, density = 1000
    END
  ENDCASE

  PRINT, ' '


END