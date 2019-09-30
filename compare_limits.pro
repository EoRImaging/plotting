;; Script to make comparison plots of the various EoR upper limits. To include/exclude certain
;; papers, uncomment the associated *papers_array, name_array, symbol_array, and symbol_name_array
;; (and optionally linestyle_array). Change n_papers to match.
;;
;; Options:
;; plot_2d: Make one 2D plot
;; plot_3d: Make a series of rotating 3D plots that can then be made into an animation
;; outdir: The output directory of the plots
;; colortable_num: Color table number for the continuous redshift color scheme
;;                 See https://harrisgeospatial.com/docs/LoadingDefaultColorTables.html
;;
;; Last update by NBarry on 24/09/19

pro compare_limits, plot_2d = plot_2d, plot_3d = plot_3d, outdir = outdir, colortable_num = colortable_num

  ;defaults
  if ~keyword_set(plot_3d) then plot_2d=1
  if ~keyword_set(outdir) then outdir = './'
  if ~keyword_set(colortable_num) then colortable_num = 34
  hubble_param=.7
  max_delta = 1E6
  min_delta = 1E3

  ;Pointers organized by paper: mK value, center of k bin/range in h Mpc^-1, and redshift
  ;*****************
  ;Number of papers to compare in the plot
  n_papers=7
  papers_array = PTRARR(n_papers,/allocate)
  name_array = STRARR(n_papers)
  symbol_array = FLTARR(n_papers)
  symbol_name_array = STRARR(n_papers)
  linestyle_array = FLTARR(n_papers)

  ;*****Paciga et al., 2013. DOI: 10.1093/mnras/stt753
  paciga_2013 = transpose([6.15E4,.5,8.6])
  ;
  *papers_array[0] = paciga_2013
  name_array[0] = 'Paciga, 2013'
  symbol_array[0] = 14
  symbol_name_array[0] = 'Filled Diamond'

  ;*****Dillon et al., 2014. DOI: 10.1103/PhysRevD.89.023002
  dillon_2014 = FLTARR(11,3)
  dillon_2014[*,0] = [2.6E7,1.16E6,8.64E5,6.7E5,1.3E6,1.28E7,5.26E7,5.67E8,4.58E6,2.93E8,6.92E8]
  dillon_2014[*,1] = [.058,.06,.063,.065,.0678,.0712,.0715,.149,.15,.15,.089]
  dillon_2014[*,2] = [11.68,10.868,10.153,9.518,8.444,7.985,7.57,7.1896,6.84,6.52,6.23]
  ;
  *papers_array[1] = dillon_2014
  name_array[1] = 'Dillon, 2014'
  symbol_array[1] = 17
  symbol_name_array[1] = 'FilledUpTriangle'

  ;*****Dillon et al., 2015. DOI: 10.1103/PhysRevD.91.123011
  dillon_2015 = FLTARR(3,3)
  dillon_2015[*,0] = [3.8E4,3.69E4,4.67E4]
  dillon_2015[*,1] = [.18,.18,.16]
  dillon_2015[*,2] = [6.4,6.8,7.25]
  ;
  *papers_array[2] = dillon_2015
  name_array[2] = 'Dillon, 2015'
  symbol_array[2] = 18
  symbol_name_array[2] = 'FilledDownTriangle'

  ;*****Beardsley et al., 2016. DOI: 10.3847/1538-4357/833/1/102
  beardsley_2016 = FLTARR(9,3)
  beardsley_2016[*,0] = [3.67E4,2.70E4,3.56E4,3.02E4,4.70E4,3.22E4,3.2E4,2.6E4,2.5E4]
  beardsley_2016[*,1] = [.231,.27,.24,.24,.20,.24,.16,.14,.14]
  beardsley_2016[*,2] = [7.1,7.1,6.8,6.8,6.5,6.5,7.1,6.8,6.5]
  ;
  beardsley_2016_short = FLTARR(6,3)
  beardsley_2016_short[*,0] = [2.70E4,3.02E4,3.22E4,3.2E4,2.6E4,2.5E4]
  beardsley_2016_short[*,1] = [.27,.24,.24,.16,.14,.14]
  beardsley_2016_short[*,2] = [7.1,6.8,6.5,7.1,6.8,6.5]
  ;
  beardsley_2016_extra_short = FLTARR(4,3)
  beardsley_2016_extra_short[*,0] = [2.70E4,3.02E4,3.2E4,2.5E4]
  beardsley_2016_extra_short[*,1] = [.27,.24,.16,.14]
  beardsley_2016_extra_short[*,2] = [7.1,6.8,7.1,6.5]
  ;
  *papers_array[3] = beardsley_2016
  name_array[3] = 'Beardsley, 2016'
  symbol_array[3] = 24;44
  symbol_name_array[3] = 'Filled Bowtie';'Filled Right Half Circle'

  ;*****Patil et al., 2017. DOI: 10.3847/1538-4357/aa63e7
  patil_2017 = FLTARR(15,3)
  patil_2017[*,0] = [131.5^2.,242.1^2.,220.9^2.,337.4^2.,407.7^2.,86.4^2.,144.2^2.,184.7^2.,$
    296.1^2.,342.0^2.,79.6^2.,108.8^2.,148.6^2.,224.0^2.,366.1^2.]
  patil_2017[*,1] = [.053,.067,.083,.103,.128,.053,.067,.083,.103,.128,.053,.067,.083,.103,.128]
  patil_2017[*,2] = [8.3,8.3,8.3,8.3,8.3,9.15,9.15,9.15,9.15,9.15,10.1,10.1,10.1,10.1,10.1]
  ;
  patil_2017_short = FLTARR(6,3)
  patil_2017_short[*,0] = [337.4^2.,407.7^2.,296.1^2.,342.0^2.,224.0^2.,366.1^2.]
  patil_2017_short[*,1] = [.103,.128,.103,.128,.103,.128]
  patil_2017_short[*,2] = [8.3,8.3,9.15,9.15,10.1,10.1]
  ;
  ;**Change the Patil array to include small k (patil_2017) or not (patil_2017_short)
  *papers_array[4] = patil_2017
  name_array[4] = 'Patil, 2017'
  symbol_array[4] = 40
  symbol_name_array[4] = 'Filled Lower Half Circle'

  ;*****Kolopanis et al, 2019. Accepted
  paper_collab = FLTARR(6,3)
  paper_collab[*,0] = [8.09e4,1.28e5,3.61e4,1.46e5,3.8e6,2.41e6]
  paper_collab[*,1] = [0.39,0.53,0.31,0.36,0.334,0.33]
  paper_collab[*,2] = [7.49,8.13,8.37,8.68,9.93,10.88]
  ;
  *papers_array[5] = paper_collab
  name_array[5] = 'Kolopanis, 2019'
  symbol_array[5] = 28
  symbol_name_array[5] = 'Filled Laying Bar'

  ;*****Barry et al., 2019. Accepted (N-S only)
  barry_2019 = FLTARR(6,3)
  barry_2019[*,1]=[0.17426223,0.20330593,0.23234963,0.26139334,0.29043704,0.31948075]
  barry_2019[*,0]=[9176.4701,3273.5603,5681.5812,9564.6415,15549.955,20607.678]
  barry_2019[*,2]=[7,7,7,7,7,7]
  ;
  barry_long_2019 = FLTARR(18,3)
  barry_long_2019[*,1]=[0.174,0.203,0.232,0.261,0.290,0.319,0.349,0.523,0.552,0.581,0.610,0.639,$
    0.668,0.697,0.726,0.755,0.929,0.958]
  barry_long_2019[*,0]=[9.22E3,3.89E3,6.32E3,9.61E3,1.56E4,2.07E4,8.32E4,7.09E4,3.2E4,4.43E4,6.38E4,$
    1.11E5,2.62E5,2.28E5,4.25E4,1.22E5,1.35E6,7.01E4]
  barry_long_2019[*,2]=[7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7]
  ;Same limits, but with infinities on either side of coarse channels and beginning and end. Forces the histogram step
  ;to be of the correct x-width without plotting extra data.
  barry_long_step_2019 = FLTARR(45,3)
  barry_long_step_2019[*,1]=[0.145,0.174,0.203,0.232,0.261,0.290,0.319,0.349,0.3776,0.4937,0.523,0.552,0.581,0.610,0.639,$
    0.668,0.697,0.726,0.755,0.784,0.900,0.929,0.958,0.987,1.017,1.046,1.075,1.104,1.133,1.162,1.191,1.220,1.249,1.278,$
    1.307,1.336,1.365,1.394,1.423,1.452,1.481,1.510,1.539,1.568,1.597]
  barry_long_step_2019[*,0]=[0,9.22E3,3.89E3,6.32E3,9.61E3,1.56E4,2.07E4,8.32E4,0,0,7.09E4,3.2E4,4.43E4,6.38E4,$
    1.11E5,2.62E5,2.28E5,4.25E4,1.22E5,0,0,1.35E6,7.01E4, 5.91E4, 3.24E5,5.32E5,2.42E5,7.41E4,8.57E4,1.13E5,4.65E5,$
    0,0,0,0,0,2.32E6,7.22E5,5.48E5,6.58E5,4.67E5,6.16E5,1.48E6,2.15E6,1.73E6]
  inds = where(barry_long_step_2019[*,0] EQ 0,n_count)
  if n_count GT 1 then barry_long_step_2019[inds,0] = !Values.F_INFINITY
  barry_long_step_2019[*,2]=[7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7]
  ;
  *papers_array[6] = barry_long_step_2019
  name_array[6] = 'Barry, 2019'
  symbol_array[6] = 0 ;16
  symbol_name_array[6] = 'Step' ;'Filled Circle'
  linestyle_array[6] = 0

  ;***** Metens et al., 2019. In prep
  ;mertens_2019 = transpose([72.4^2,0.075,10.1])
  ;
  ;*papers_array[8] = mertens_2019
  ;name_array[8] = 'Mertens, 2019 (in prep)'
  ;symbol_array[8] = 15
  ;symbol_name_array[8] = 'Filled Square'
  ;
  ;*****************

  ;Fiducial EoR
  ;*****************
  ;Not currently included
  ;GalaxyParam_LF+NF+Tau by Brad Greg
  ;Fourier mode (k [Mpc^-1]) and 95th percentiles on the power spectrum in units of mK^2
  brad_k = [0.02094395,0.03183814,0.0477084,0.06481693,0.0846983,0.1122589,0.1512655,0.2037167,$
    0.2749778,0.3717481,0.5018322,0.6776517,0.9149306,1.235268,1.6676,2.207298,2.79783]
  brad_k = brad_k/hubble_param
  brad_95 = [3.00633765723,12.8642551767,18.6682120582,38.322192403,43.5651253669,49.1075159329,$
    51.9213041607,48.790639303,31.3894142399,26.2910138229,26.8292876974,26.5396841581,26.0301170622,$
    27.1653881665,31.5002616162,42.7424864461,60.1632270707]
  n_k=n_elements(brad_k)
  brad_k=(brad_k[1:(n_k-1)]+brad_k[0:(n_k-2)])/2. ;recenter for histogram look

  ;k (Mpc^-1) and fiducial power spectrum in units of mK^2
  k_centers = [2.094395e-02,3.183814e-02,4.770840e-02,6.481693e-02,8.469830e-02,1.122589e-01, $
    1.512655e-01,2.037167e-01,2.749778e-01,3.717481e-01,5.018322e-01,6.776517e-01,9.149306e-01,$
    1.235268e+00,1.667600e+00,2.207298e+00,2.797836e+00]
  k_centers = k_centers/hubble_param
  delta_eor = [3.148450e-01,1.371862e+00,2.112182e+00,4.716543e+00,6.822401e+00,1.036046e+01,$
    1.477729e+01,9.501991e+00,8.793401e+00,7.287165e+00,6.430653e+00,5.709567e+00,5.390083e+00,$
    5.408442e+00,6.280324e+00,8.959497e+00,1.353618e+01]
  ;*****************


  for papers_i=0, n_papers-1 do begin
    if papers_i EQ 0 then begin
      redshift_arr = (*papers_array[papers_i])[*,2]
      delta_arr = (*papers_array[papers_i])[*,0]
    endif else begin
      redshift_arr = [redshift_arr, (*papers_array[papers_i])[*,2]]
      delta_arr = [delta_arr, (*papers_array[papers_i])[*,0]]
    endelse
  endfor
  inds = where(delta_arr GT max_delta,n_count)
  if n_count GT 0 then redshift_arr[inds] = !Values.F_INFINITY
  max_redshift= max(redshift_arr,/NAN)
  min_redshift=min(redshift_arr,/NAN)
  ncolors=(round(max_redshift*10.) - round(min_redshift*10.))+1
  cgLOADCT, colortable_num;, /reverse

  xtitle_use='k (h Mpc$\up-1$)'
  ytitle_use='$\Delta$$\up2$ (mK$\up2$)'

  ;***************** 2D plotting
  if keyword_set(plot_2d) then begin

      cgPS_Open, outdir + 'limits_comparison_2D.pdf',/quiet,/nomatch

      ;Main plot setup -- No data included
      cgplot, (*papers_array[0])[0,1],(*papers_array[0])[0,0], /xlog,/ylog, xrange=[.04,2.0], yrange=[min_delta,max_delta],$
        aspect=.5, YGridStyle=1,YTicklen=1.0, xtitle=xtitle_use, ytitle = ytitle_use,/NoData, charsize=1
        

      ;Find list of papers for histogram plotting option.
      ;step_inds is -1 if not histogram plotting option
      step_inds = where(symbol_name_array EQ 'Step',n_step_count)

      ;Loop over papers
      for papers_i=0, n_papers-1 do begin
        arr = *papers_array[papers_i]
        redshift_colors = BytScl(round((arr[*,2]-min_redshift)/(max_redshift-min_redshift)*256.),min=0,max=256,top=256)

        ;Symbol plotting option
        tmp = where(papers_i EQ step_inds,n_count)
        if n_count EQ 0 then begin
          for arr_i=0,N_elements(arr[*,0])-1 do begin

            cgoplot, arr[arr_i,1],arr[arr_i,0],psym=symbol_array[papers_i],thick=3, color=redshift_colors[arr_i], /xlog,/ylog, $
              symsize=2.0,charsize =1.
          endfor

          ;Histogram plotting option
        endif else begin

          for red_i=0, N_elements(uniq(arr[*,2]))-1 do begin
            redshift_inds = where(arr[*,2] EQ arr[(uniq(arr[*,2]))[red_i],2])

            cgoplot, arr[redshift_inds,1],arr[redshift_inds,0],psym=10,thick=8, color=redshift_colors[redshift_inds[0]],$
              /xlog,/ylog, charsize =1., linestyle=linestyle_array[papers_i]

          endfor

        endelse
      endfor

      ;Redshift legend on the right-hand side
      cgcolorbar,range=[min_redshift,max_redshift], position=[0.98, 0.27, 0.999, 0.75],/vertical,charsize=1, title = 'redshift'

      two_rows = 1
      ;three_rows=1

      for papers_i=0,n_papers-1 do begin

        locations = FLTARR(2)
        if keyword_set(two_rows) then begin
          ;2 legend rows
          locations[0] = 0.1 + (papers_i/2)*((0.95 - 0.1)/round(n_papers/2.))
          locations[1] = 0.16 - (papers_i mod 2)*.05
        endif else if keyword_set(three_rows) then begin
          ;3 legend rows
          locations[0] = 0.1 + (papers_i/3)*((0.95 - 0.1)/round(n_papers/3.))
          locations[1] = 0.16 - (papers_i mod 3)*.05
        endif

        if symbol_name_array[papers_i] NE 'Step' then $
          al_legend, name_array[papers_i], color='grey', pos=locations, /normal, psym=symbol_array[papers_i], $
          charsize=1, linsize=0, box=0, symsize=1.5 $
        else $
          al_legend, name_array[papers_i],color='grey',pos=locations,linsize=.3,thick=4,charsize=1,$
          linestyle=linestyle_array[papers_i],box=0,/normal
      endfor

      cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
      stop

    ;***************** 3D plotting
  endif else if keyword_set(plot_3d) then begin

    n_rot=160
    for rotz_i=0,n_rot-1 do begin

      cgPS_Open, outdir + "limits_comparison_3D_"+string(strtrim(rotz_i,2),FORMAT='(I03)')+".png",$
        FONT=1, Charsize=3.0,/quiet,/nomatch

      ; Set the 3D coordinate space with axes.
      ;x = k
      ;y = z
      ;z = mk^2
      rotz_j = rotz_i
      if rotz_i GT n_rot/2 - 1 then rotz_j = n_rot - rotz_i
      cgSurf, DIST(5), /NODATA, /SAVE, XRANGE=[.04,1], $
        YRANGE=[6,11], ZRANGE=[1000, 1e6], XSTYLE=1, $
        YSTYLE=1, ZSTYLE=1, CHARSIZE=2.0, $
        POSITION=[0.1, 0.1, 0.95, 0.95, 0.1, 0.95], $
        XTICKLEN=1, YTICKLEN=1, XGRIDSTYLE=1, YGRIDSTYLE=1, $
        xtitle=xtitle_use, ytitle='redshift',ztitle=ytitle_use,title=title, $
        /zlog, /xlog, rotz=30. - float(rotz_j)/4.
      cgAXIS, XAXIS=1, /T3D, CHARSIZE=2.0
      cgAXIS, YAXIS=1, /T3D, CHARSIZE=2.0

      for papers_i=0, n_papers-1 do begin
        arr = *papers_array[papers_i]
        z = arr[*,0]
        x = arr[*,1]
        y = arr[*,2]
        FOR j=0,N_elements(z)-1 DO begin
          if z[j] GT 1.0E6 then continue
          cgPlotS, x[j], y[j], z[j], COLOR=color_array[papers_i], PSYM=sym_array[papers_i], SYMSIZE=2.5, /T3D
        ENDFOR
        FOR j=0,N_elements(z)-1 DO begin
          if z[j] GT 1.0E6 then continue
          cgPlotS, [x[j], x[j]], [y[j], y[j]], [1000, z[j]], COLOR=color_array[papers_i], /T3D
        ENDFOR


      endfor

      cglegend, Title=name_array[[0,1]], color=color_array[[0,1]], location=[.02,.06], psym=sym_num[[0,1]], $
        charsize=1.25,Length=0,symsize=1.5
      cglegend, Title=name_array[[2,3]], color=color_array[[2,3]], location=[.20,.06], psym=sym_num[[2,3]], $
        charsize=1.25,Length=0,symsize=1.5
      cglegend, Title=name_array[[4,5]], color=color_array[[4,5]], location=[.41,.06], psym=sym_num[[4,5]], $
        charsize=1.25,Length=0,symsize=1.5
      cglegend, Title=name_array[[6,7]], color=color_array[[6,7]], location=[.60,.06], psym=sym_num[[6,7]], $
        charsize=1.25,Length=0,symsize=1.5
      cglegend, Title=name_array[[8]], color=color_array[[8]], location=[.80,.06], psym=sym_num[[8]], charsize=1.25,$
        Length=0,symsize=1.5

      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    endfor
  endif

end
