;; Create publication quality EoR upper limits with a histogram scheme. Example: see Barry et al 2019b
;; Can input either eppsilon sav files or CHIPS csv files. No further analysis is done...this is
;; just a plotting script. Input normalizations are assumed to be correct.
;;
;; Inputs:
;; input_file_dir: Directory of the location of the eppsilon sav file or the CHIPS csv file
;; sav_file_name: Name of the eppsilon sav file. If not set, script will try to create default name
;; csv_file:
;; pols: Specify polarization names to plot (default is xx and yy). Currently can only do two pols
;; missing_inds: Indices of missing and/or data that should not be plotted. Used for excluding coarse
;;   band contamination lines, and is dependent on resolution. Barry et al 2019b example:
;;   missing_inds=[7,8,9,10,11,21,22,23,24,25,36,37,38,39,40,50,51,52,53]
;; outdir: Output directory (default current directory)
;; low_y_range: Lowest y-axis value in mK^2 (default 10^0)
;; high_y_range: Highest y-axis value in mK^2 (default 10^8)
;; low_x_range: Lowest x-axis value in h Mpc^-1 (default 0.14)
;; high_x_range: Highest x-axis value h Mpc^-1 (default 1.9)
;; theory_rgb: RGB array for the theory (default [197,120,62], brown)
;; upper_lim_rgb: RGB array for the upper limits (default [73, 142, 217], light blue)
;; sim: Option to only plot the measured power, useful for simulations without noise 
;; pdf: Create a pdf output (default not set)

; ******** Internal inverf function
function inverf, erfx

  on_error, 2

  ; ERROR FUNCTION CAN ONLY RETURN 0 <= ERRORF(X) <= 1, SO
  ; THE INVERSE ERROR FUNCTION OF ANYTHING GREATER THAN ONE IS
  ; MEANINGLESS...
  if (total((erfx le 1d0) AND (erfx ge 0)) ne N_elements(erfx)) then $
    message, 'The range of the error function is 0 < erf(x) <= 1!'

  ; IS ERFX WITHIN MACHINE PRECISION OF UNITY...
  type = size(erfx,/TYPE)
  epsn = (machar(DOUBLE=(type eq 5L))).epsneg

  ; CHECK THE MACHINE PRECISION HERE...
  p = 1d0 - (erfx - epsn*( (erfx - 0.5*epsn) eq fix(1.,TYPE=type) ))
  p = 0.5d0*p
  t = sqrt(-2d0*alog(p))
  num =  (0.010328d0*t + 0.802853d0)*t + 2.515517d0
  den = ((0.001308d0*t + 0.189269d0)*t + 1.432788d0)*t + 1d0
  return, 0.70710678118654752440d0 * ((t - num/den) > 0)

end
; ******** End internal inverf function


; ******** Main histogram_limits_pub plotting script
pro histogram_limits_pub, input_file_dir=input_file_dir, sav_file_name = sav_file_name, csv_file=csv_file, pols=pols, $
  missing_inds=missing_inds, outdir = outdir, low_y_range=low_y_range, high_y_range=high_y_range, $
  low_x_range=low_x_range, high_x_range=high_x_range, theory_rgb=theory_rgb, upper_lim_rgb=upper_lim_rgb, $
  sim=sim, pdf=pdf


  if ~keyword_set(sav_file_name) AND ~keyword_set(csv_file) then $
    message, 'Please set sav_file_name (eppsilon output) or csv_file (CHIPS output).'

  if ~keyword_set(pols) then pols = ['xx','yy']
  if ~keyword_set(input_file_dir) then input_file_dir = './'
  if ~keyword_set(outdir) then outdir = './'
  if ~keyword_set(output_file) then output_file='histogram_upper_limit.pdf'

  ;Sav file name default
  if ~keyword_set(sav_file_name) and ~keyword_set(csv_file) then begin
    int_name = 'btl_noalltv_noocc4'
    chans='ch9-126_'
    cubes='res_'
    basefile='/Combined_obs_' + int_name + '_cubeXX__even_odd_joint_fullimg_' + chans + cubes
    kperp_end='80'
    kperp_start='18'
    horizon='120'
    endfile='_averemove_swbh_dencorr_' + 'no_'+horizon+'deg_wedge_kperplambda'+kperp_start+'-'+kperp_end+'_kpar0.15-200_1dkpower.idlsave'
  endif else if keyword_set(sav_file_name) and ~keyword_set(csv_file) then begin
    if strmatch(sav_file_name,'*'+pols[0]+'*') then begin
      split_strings = strsplit(sav_file_name, /EXTRACT, pols[0])
      basefile = split_strings[0]
      endfile = split_strings[1]
    endif else if strmatch(sav_file_name,'*'+pols[1]+'*') then begin
      split_strings = strsplit(sav_file_name, /EXTRACT, pols[1])
      basefile = split_strings[0]
      endfile = split_strings[1]
    endif else message, 'sav_file_name does not have a polarization identifier'
  endif

  ;Plotting range defaults
  if ~keyword_set(low_y_range) then low_y_range=10^0.
  if ~keyword_set(high_y_range) then high_y_range=10^8.
  if ~keyword_set(low_x_range) then low_x_range=.14
  if ~keyword_set(high_x_range) then high_x_range=1.9

  limit_percent = 0.97725 ;2 sigma
  hubble_param=0.71

  ;Brad Greig's fiducial theory and confidence limits for redshift 7
  ;GalaxyParam_LF+NF+Tau
  ;Fourier mode (k [Mpc^-1]), 68th, 95th and 99th percentiles on the power spectrum in units of mK^2
  brad_k = [0.02094395,0.03183814,0.0477084,0.06481693,0.0846983,0.1122589,0.1512655,0.2037167,0.2749778,0.3717481,0.5018322,0.6776517,0.9149306,1.235268,1.6676,2.207298,2.79783]
  brad_68 = [2.06836030818,9.33650238046,13.5354885655,27.9094808815,31.8379169811,35.6549002097,$
    37.7766299375,34.7228577191,22.5256589937,19.1373452916,19.0715957566,18.8232315706,18.1639861695,$
    18.6493671092,21.4087292929,28.7580668946,40.4441293602]
  brad_95 = [3.00633765723,12.8642551767,18.6682120582,38.322192403,43.5651253669,49.1075159329,$
    51.9213041607,48.790639303,31.3894142399,26.2910138229,26.8292876974,26.5396841581,26.0301170622,$
    27.1653881665,31.5002616162,42.7424864461,60.1632270707]
  brad_99 = [3.277767,13.62931,19.91,40.66932,46.2304,52.1699,55.14123,52.86583,34.2097,28.67557,29.69817,29.57222,$
    29.46968,30.93632,35.9361,48.86067,68.80418]

  ;k (Mpc^-1) and fiducial
  brad_k_2 = [2.094395e-02,3.183814e-02,4.770840e-02,6.481693e-02,8.469830e-02,1.122589e-01, 1.512655e-01,2.037167e-01,2.749778e-01,3.717481e-01,$
    5.018322e-01,6.776517e-01,9.149306e-01,1.235268e+00,1.667600e+00,2.207298e+00,2.797836e+00]
  brad_fiducial = [3.148450e-01,1.371862e+00,2.112182e+00,4.716543e+00,6.822401e+00,1.036046e+01,1.477729e+01,9.501991e+00,8.793401e+00,7.287165e+00,$
    6.430653e+00,5.709567e+00,5.390083e+00,5.408442e+00,6.280324e+00,8.959497e+00,1.353618e+01]

  brad_k_2 = brad_k_2/hubble_param
  brad_k = brad_k/hubble_param


  ;Color defaults and setting
  if ~keyword_Set(theory_rgb) then theory_rgb = [197,120,62]
  if ~keyword_set (upper_lim_rgb) then upper_lim_rgb = [73, 142, 217]
  ;73, 142, 217 little brighter blue, 54,107,163 little darker blue
  TVLCT, theory_rgb[0], theory_rgb[1], theory_rgb[2], 10
  theory_color = 10B
  TVLCT, upper_lim_rgb[0], upper_lim_rgb[1], upper_lim_rgb[2], 12
  upper_lim_color = 12B

  if keyword_set(pdf) then thickness = 5 else thickness = 3


  ; ******* Plotting for loop through polarizations
  for j=0,N_elements(pols)-1 do begin

    ;Restore eppsilon sav file
    if keyword_set(sav_file_name) then begin
      file = input_file_dir + basefile + pols[j] + endfile
      restore,file

      ;Construct k center from k edges
      n_k_edges=n_elements(k_edges)
      k=(k_edges[1:(n_k_edges-1)]+k_edges[0:(n_k_edges-2)])/2.
      n_k=n_elements(k)

      ;Construct power in mK^2 units
      delta=power*(k^3.)/(2.*!pi^2.)

      ;Construct sigma in mK^2 units.
      dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)

      ;Construct k in h Mpc^-1
      k=k/hubble_param
    endif

    if keyword_set(csv_file) then begin
      ;expected file header:
      ;k [h Mpc^-1],P XX [Delta^2],P YY [Delta^2],xerrhigh [h Mpc^-1],xerrlow [h Mpc^-1],thermal P XX [Delta^2],thermal P YY [Delta^2]
      ;Code does not currently include x error bars
      csv_struct = READ_CSV(csv_file, HEADER=header, types='Double')
      
      k_edges = csv_struct.field1
      n_k_edges=n_elements(k)
      k=(k_edges[1:(n_k_edges-1)]+k_edges[0:(n_k_edges-2)])/2.
      n_k=n_elements(k) 
     
      if j EQ 0 then begin
        delta = csv_struct.field2
        dsigma = csv_struct.field6 / 2. ;Defined by default as 2 sigma
      endif else begin
        delta = csv_struct.field3
        dsigma = csv_struct.field7 / 2. ;Defined by default as 2 sigma
      endelse
    endif

    ;Some versions of CSV reader do not understand NAN or Infinity
    nans = where(delta EQ 'NaN',n_count)
    if n_count GT 0 then delta[nans]=0 
    nans = where(delta EQ '-NaN',n_count)
    if n_count GT 0 then delta[nans]=0 
    nans = where(dsigma EQ 'Infinity',n_count)
    if n_count GT 0 then dsigma[nans]=0 
    dsigma = Double(dsigma)
    delta = Double(delta)

    ;Erf function not defined for 0, will error
    zeros = where(dsigma EQ 0,n_count)
    if n_count GT 0 then dsigma[zeros] = !Values.F_INFINITY
    zeros = where(delta EQ 0,n_count)
    if n_count GT 0 then delta[zeros] = !Values.F_INFINITY
    i_nan = where(~finite(delta),n_count,complement=i_defined)
   
    delta_def = delta[i_defined]
    dsigma_def = dsigma[i_defined]

    ;Construct upper limits using two sigma prior
    limits_def=dsigma_def*(inverf(limit_percent-(1.-limit_percent)*erf((delta_def)/dsigma_def/sqrt(2)))*sqrt(2))+(delta_def)
    limits = DBLARR(N_elements(delta))
    limits[*] = !Values.F_NAN
    limits[i_defined] = limits_def

    ;Print lowest limit to screen
    lim=min(limits,ind,/NAN)
    header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind])+' h Mpc^-1 for '+pols[j]
    print,header

    ;Make missing inds NAN for plotting purposes
    if keyword_set(missing_inds) then begin
      delta[missing_inds] = !Values.F_NAN
      dsigma[missing_inds] = !Values.F_NAN
      limits[missing_inds] = !Values.F_NAN
    endif

    if (j EQ 0) and keyword_set(pdf) then cgPS_Open,outdir + output_file,/quiet,/nomatch

    ;Plot positions
    position1=[0.1, 0.6, 0.5, 0.9]
    position2=[0.1, 0.45, 0.5, 0.6]
    position4=[0.5, 0.6, 0.9, 0.9]
    position5=[0.5, 0.45, 0.9, 0.6]

    if (j EQ 0) then begin
      position_use=position1
      xtitle_use='k (h Mpc$\up-1$)'
      ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
      XTICKFORMAT_use="(I0)"
      YTICKFORMAT_use="(I0)"
    endif

    if (j EQ 1) then begin
      position_use=position4
      xtitle_use='k (h Mpc$\up-1$)'
      ytitle_use=''
      XTICKFORMAT_use="(I0)"
      YTICKFORMAT_use="(A1)"
    endif

    
    if (j EQ 0) then begin
      ;Create plot shell without data. Add pol name and legend
      cgplot, k,limits,/xlog,/ylog,psym=10, xrange=[low_x_range,high_x_range], yrange=[low_y_range,high_y_range],ytitle=ytitle_use, $
        xtitle=xtitle_use, charsize =.8, color=upper_lim_color,thick=thickness,position=position_use,/noerase,/nodata
        
      if N_elements(pols) EQ 2 then begin    
        XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
        al_legend, ['fiducial theory','95% confidence'],color=[theory_color,theory_color], pos=[.5,.89], charsize=.8, thick=thickness,$
          linestyle=[0,2],linsize=.7, box=0, /normal
      endif else if N_elements(pols) EQ 1 then begin
        if pols EQ 'xx' then pol_name='E-W' else if pols EQ 'yy' then pol_name='N-S'
        XYOuts, position_use[0]+0.2, position_use[3]+.02, pol_name+', z=7', /Normal, Alignment=0.5, Charsize=.8
        if ~keyword_set(sim) then begin
          al_legend, ['measured power','2$\sigma$ upper limit','noise level','fiducial theory','95% confidence'],$
            color=['black',upper_lim_color,upper_lim_color,theory_color,theory_color], pos=[.1,.89], $
            charsize=.7, thick=thickness,linestyle=[0,0,2,0,2],linsize=.8, box=0, /normal
        endif else begin
          al_legend, ['measured power','fiducial theory','95% confidence'],color=[upper_lim_color,theory_color,theory_color],$
            pos=[.1,.89], charsize=.7, thick=thickness,linestyle=[0,0,2],linsize=.8, box=0, /normal      
        endelse  
      endif
    endif
    if (j EQ 1) then begin
      ;Create plot shell without data. Add pol name and legend
      cgplot, k,limits,/xlog,/ylog,psym=10, xrange=[low_x_range,high_x_range], yrange=[low_y_range,high_y_range],ytitle=ytitle_use, $
        xtitle=xtitle_use, charsize =.8, color=upper_lim_color,thick=thickness,position=position_use,YTICKFORMAT=YTICKFORMAT_use,/noerase,/nodata
      XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
      if ~keyword_set(sim) then begin
        al_legend, ['measured power','2$\sigma$ upper limit','noise level'],color=['black',upper_lim_color,upper_lim_color], pos=[.1,.89], $
          charsize=.7, thick=thickness,linestyle=[0,0,2],linsize=.8, box=0, /normal
      endif else begin
        al_legend, ['measured power'],color=[upper_lim_color], pos=[.1,.89], $
          charsize=.7, thick=thickness,linestyle=[0],linsize=.8, box=0, /normal      
      endelse
    endif

    ;Construct error bars. Force them to end at plot boundaries for aesthetics
    inds1 = where(k GT low_x_range)
    inds2 = where(k[inds1] LT high_x_range)
    k_in_plot = k[inds1[inds2]]
    delta_k=(k_in_plot[2]-k_in_plot[1])/2.
    error_low = (dsigma[inds1[inds2]]*2.)
    lows = where(delta[inds1[inds2]] LT 0, n_count)
    if n_count GT 0 then error_low[lows]=abs(delta[inds1[inds2[lows]]])
    lows = where((abs(delta[inds1[inds2]]) - error_low) LT low_y_range,n_count)
    if n_count GT 0 then error_low[lows]=abs(delta[inds1[inds2[lows]]])-low_y_range


    if ~keyword_set(sim) then begin
      ;Plot error bars
      for k_i=0, n_elements(k_in_plot)-2 do begin
        delta_k_high=(k_in_plot[k_i+1]-k_in_plot[k_i])/2.
        if k_i NE 0 then delta_k_low=(-k_in_plot[k_i-1]+k_in_plot[k_i])/2. else delta_k_low = delta_k_high

        cgColorFill, [k_in_plot[k_i]-delta_k_low, k_in_plot[k_i]+delta_k_high, k_in_plot[k_i]+delta_k_high,k_in_plot[k_i]-delta_k_low], $
          [limits[inds1[inds2[k_i]]], limits[inds1[inds2[k_i]]], abs(delta[inds1[inds2[k_i]]])-error_low[k_i],abs(delta[inds1[inds2[k_i]]])-error_low[k_i]], $
          Color='grey',/checkforfinite
      endfor

      ;Plot limits, thermal noise, and power
      cgoplot, [k[0]-(k[1]-k[0]),k,k[n_k-1]+(k[n_k-1]-k[n_k-2])],[!Values.F_NAN,limits,!Values.F_NAN],$
        /xlog,/ylog,psym=10,color=upper_lim_color,thick=thickness,position=position_use
      cgoplot, [k[0]-(k[1]-k[0]),k,k[n_k-1]+(k[n_k-1]-k[n_k-2])],[!Values.F_NAN,dsigma,!Values.F_NAN],$
        /xlog,/ylog,psym=10,linestyle=2, color=upper_lim_color,thick=thickness,position=position_use
      cgoplot, [k[0]-(k[1]-k[0]),k,k[n_k-1]+(k[n_k-1]-k[n_k-2])],[!Values.F_NAN,delta,!Values.F_NAN],$
        /xlog,/ylog,color='black',psym=10,thick=thickness-1,position=position_use

    endif else begin
      ;Just plot the power for simulation comparisons
        cgoplot, [k[0]-(k[1]-k[0]),k,k[n_k-1]+(k[n_k-1]-k[n_k-2])],[!Values.F_NAN,delta,!Values.F_NAN],$
          /xlog,/ylog,color=upper_lim_color,psym=10,thick=thickness-1,position=position_use
    endelse


    ;Plot theory
    cgoplot, brad_k_2,brad_fiducial,/xlog,/ylog,color=theory_color,thick=thickness-1,position=position_use
    cgoplot, brad_k,brad_95,/xlog,/ylog,color=theory_color,thick=thickness-1,position=position_use,linestyle=2


  endfor
  ; ******* End plotting for loop through polarizations

  if keyword_set(pdf) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage 
end
; ********

