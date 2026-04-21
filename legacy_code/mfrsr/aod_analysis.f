      program aod_analysis
c
c This program operates on one day of visible MFRSR data to produce
c a daily time series of cloudmasked aerosol optical depth.  
c
c It applies a Vo and associated error interpolated to the day
c being processed. The interpolation is done with equations that 
c describe the variation of Vo and error over a period of an 
c MFRSR's deployment.  
c
c The equations from which Vo and error are derived are derived 
c by running the program mfr_6ch_calibrate.f on a monthly or bi-monthly
c basis and deriving a mean Vo and error for that period.  
c After this is done for the duration of an MFRSR deployment, or some
c period, e.g., 2 years, the sample of mean Vo's are subjected to a fit.
c The fitted equation has  constant, linear (slope), sine, and cosine 
c terms.  Periodic terms were included because of a repeatable 
c cyclical behavior observed in the 2-month Vo time series' of all
c MFRSRs tested that resembled a seasonal temperature variation. 
c Errors in the mean Vo's are fit to a linear function. In both,
c the time variable used is in the yyyy.ffff format.  For example,
c 1999.5000 would be June 30, 1999.  In applying the Vo equation, 
c that date must first be converted to an angle by multiplying
c by 2*pi.  
c 
c In this program, those equations are used to interpolate 
c a representative Vo and error to the day being processed.  
c Such equations are are derived for Vo's of each MFRSR channel. 
c The interpolations are made in sub get_izero.
c
c Total optical depth is computed using the Vo and sigma 
c generated in get_izero. It then corrects the Vo for the actual
c earth-sun distance for the day being processed. Rayleigh  
c scattering tau is computed dynamically for each time of day 
c that is processed by accessing the measured station pressure at 
c the SURFRAD station.  Ozone absorption is accessed from a file
c that is written by a Perl script that gets daily ozone over the 
c station being processed from NASA TOMS OMI or OMPS. 
c Daily total ozone is listed sequentially in a file for each 
c station for the duration of the network. Those files are named:
c
c bon_ozone.dat
c fpk_ozone.dat
c gwn_ozone.dat
c tbl_ozone.dat
c dra_ozone.dat
c psu_ozone.dat
c sxf_ozone.dat
c
c In this program total ozone for the particular day being processed
c is retrieved from those files through interpolation, because, not 
c all days may be represented.
c 
c---------------------------------------------------------------
c
c AOD is computed by subtracting tau_Rayleigh and tau_ozone from
c the TOD. After the AOD is computed, all values greater than 3.0
c are removed.  Also, the Angstrom exponent is computed for each
c data time using the 500nm and 870nm channel AODs.    
c
c Additional cloud screening is done by applying Alexandrov's
c method (Alexandrov et al. 2004, GRL).  Points that pass that
c test are marked with an icloud value of 0.  Points that do not
c pass are marked with an icloud value of 1.  The icloud value 
c is written on each line of the output file.   
c
c Results are written to an output file named with the convention:
c [sta]_yyyymmdd.aod, where sta is the three-letter identifier of the 
c station, yyyy is the 4-digit year, mm is the month number, 
c and dd is the two-digit day of month.
c
c The aerosol optical depth product files are placed in the directory
c structure /q/surfrad/aod/[sta]/[yyyy], where sta is the 3-letter
c identifier of the station (bon for Bondville, dra for Desert Rock,
c etc), and yyyy is the 4-digit year.  For example, the aod data for
c Table Mountain on 17-apr-2001 will be in the file:
c /q/surfrad/aod/tbl/2001/tbl20010417.aod
c
c In a modification made on 29-Nov-2011 the BND and FPE SURFRAD MFRSRs
c were added as an option for the USDA MFRSRs that were historically
c used at those sites.  The output files from these MFRSR data still 
c use the bon and fpk station identifiers.
c
c On 19-Feb-2015, all problem data identifiers were moved from this program
c to mod files unique to each station, and were named [sta]_aod_errors.dat. 
c Code to read those files and apply the specified data deletions
c was added.
c
c 1625 channel processing:
c
c In 2022 the way of computing 1625 channel optical depths was changed.
c A very high resolution absorption spectrum was used to compute optical
c depths at 1625 for co2, ch4 and h2o more accurately. Now, instead
c of scaling the 5-cm h2o od by the square root of the actual value, it
c is linearly scaled to the measured value. Water vapor optical depth 
c is computed from the 1200 and following 0000 UTC soundings interpolated
c to the SURFRAD sites and then interpolated linearly between those times
c to time being processed.  Also, a dependence of co2, ch4 and h2o optical
c depths on air mass is described by a quadratic fit as a function to air
c mass for each trace gas species. Also, the trace gas optical depths of co2
c and ch4 are now scaled by the ratio of the measured surface pressure to 1013 mb.
c The value of co2 used in RT modeling iof the 1625 filters has been 
c updated to 410 ppm.
c
      parameter (nsta=10)
      parameter (nchannel=5)
      parameter (mfrdata=1500)
      parameter (nvals=5000)
      parameter (nbaro=3000)
c
      common/channels/wave_l(6),rayleigh(nchannel),mfrsr_sn,od_co2,
     1od_ch4,od_h2o_5cm,pwv_snd(4),i4time_pwv(4)
      common/unit/lun_mfr,lun_mfrhead,lun_aod,lun_o3,lun_p,lun_cof,
     1lun_hom,lun_err,lun_snd
      common/directory/in_directory,in_directory_mfr,out_directory,
     1 in_directory_pres,in_directory_o3,home_directory,snd_directory,
     2 lstring_in,lstring_mfr,lstring_out,lstring_p,lstring_o3,len_hom,
     3 len_snd
      common/files/isite,extensions(nsta),names(nsta),lstring(nsta),
     1i4time_date,ista_to_gmt(nsta),pressure_nominal(nsta),xlat,xlon
      common/baro/real_i4tim70_loc(nbaro),p(nbaro),n_pres
      common/calibration/i_zero(nchannel),error(nchannel),
     1  constant_coef(nchannel),slope_coef(nchannel),
     2  sine_coef(nchannel),cos_coef(nchannel),
     3  err_slope_i_zero(nchannel),err_intercept_i_zero(nchannel)
      common/ch1625/iserial(20),number_batch(20),batch(20),channel_1625
c define the fit coefficients (intercept, x and x^2) of 1625 od as a 
c function of air mass for co2, ch4 and h2o
      common/coefficients/ch4int,ch4x,ch4x2,co2int,co2x,co2x2,
     1h2oint,h2ox,h2ox2
c
      logical*1 channel_1625
c
      integer*4 ldate_mfr(mfrdata),ltime_mfr(mfrdata)
      integer*4 i4time70_local(mfrdata)
      integer*4 i4time_date
      integer*4 number_batch
c
      real*4 wave(596),o3coeff(596)
      real*4 i_zero
      real*4 intercept_i_zero
      real*4 direct(mfrdata,nchannel)
      real*4 cosz_mfr(mfrdata)
      real*4 ln_direct(mfrdata,nchannel)
      real*4 path_len_mfr(mfrdata) ! path lengtgh for all MFRSR data lines
      real*4 path(mfrdata) ! path length for all accepted MFRSR data
      real*4 error_aod(mfrdata,nchannel)
      real*4 aod(mfrdata,nchannel)
      real*4 alpha(mfrdata)
      real*4 alpha_clear(mfrdata)
      real*4 ozone_coef(nchannel)
      real*4 t_o(nchannel)
      real*4 aod_clear(mfrdata,nchannel)
      real*4 error_aod_clear(mfrdata,nchannel)
      real*4 pres_int(mfrdata),pres_int_clear(mfrdata)
c
      real*4 snd_time(4),pv(4)
c
c Alexandrov cloud screen algorithm arrays
      real*4 tau_prime(mfrdata)
      real*4 e_prime(mfrdata)
      real*4 aodmax(mfrdata)
      real*4 aodmin(mfrdata)
      real*4 time(mfrdata)
      real*4 aod_max_interp(mfrdata)
      real*4 aod_min_interp(mfrdata)
c
      integer*4 icloud(mfrdata)
      integer*4 ltime(mfrdata),ltime_clear(mfrdata)
      integer*4 i4time70_mfr_local(mfrdata)
c
      character*1 batch
      character*1 ans
      character*1 answer
      character*3 ch_ext
      character*24 names
      character*24 a_date
      character*24 atime_1min
      character*24 atime_psu
      character*24 atime_fpk
      character*3 extensions
      character*80 in_directory,in_directory_mfr,out_directory
      character*80 in_directory_pres,in_directory_o3
      character*80 home_directory
      character*80 snd_directory
      character*4 ch_yyyy
      character*2 ch_yy,ch_mm,ch_dd
      character*80 full_filename_mfr
      character*80 full_filename_o3coef
      character*80 ozone_file
      character*80 aod_filename
      character*8 mfr_serial
      character*9 fnam_date
      character*9 fnam_tbl_fire_90610
      character*9 fnam_tbl_fire_90710
      character*9 fnam_tbl_fire_102120
c
c declarations for the error listings files and associated deletion code
c
      character*18 error_file
      character*5 ch_begin_date
      character*5 ch_end_date
      character*1 ch_interval
      character*9 ch_begin_date_full
      character*9 ch_end_date_full
c
c declarations for header records of the output file
c
      character*6 ch_wave_l(5)
      character*5 ch_error(5)
      character*3 ch_iwave
      character*4 ch_iwave4
      character*80 site_name
c
c depug file for 1625 channel trace gas optical depths
c open a temporary file for recording 1625 optical depths
      lun_od = 24
      open(unit=lun_od,file='1625_od',status='unknown')
      write(lun_od,51)
  51  format(1x,'localtime path_len od_co2 od_ch4 od_h2o Rayleigh')
c
c Define the station and date to process before defining directories
c
d     print *,'call initial'
      call initial(istatus_init)
c
c
c Define the directories where to find the input data and write the output
c
d     print *,'call define directories'
      call define_directories
c
c Get the latitude and longitude of the selected station for the SZA calc.
c
      ch_ext = extensions(isite)
      if(extensions(isite).eq.'bnd') then
         ch_ext = 'bon'
      else if(extensions(isite).eq.'fpe') then
         ch_ext = 'fpk'
      endif
c
      open(unit=lun_hom,file=home_directory(1:len_hom)//
     1ch_ext//'head',status='old',err=35)
c
      read(lun_hom,31) site_name
  31  format(a)
      read(lun_hom,*) xlat,xlon,elev
      close(unit=lun_hom)
      go to 37 !successful reading the lat,lon and elevation from the header file
c
  35     print *,'ERROR'
         print *,'couldnt open the header file for this station to get'
         print *,'the latitude, longitude and elevation'
         stop
c 
  37  continue
c
c
c
c Read the ozone absorption coefficient file
c
      full_filename_o3coef = in_directory(1:lstring_in)//'ozone.coefs'
      open(unit=lun_cof,file=full_filename_o3coef,status='old',readonly,
     1     err=255)
      do i = 1,596
         read(lun_cof,*) wave(i),o3coeff(i)
      enddo
      close(unit=lun_cof)
c
c
c Get the station pressure time series for the chosen date
c
      call i4tim70_fnam(i4time_date,fnam_date,istatus)
c
c Read the pressure for the SURFRAD file corresponding to i4time_date
c and the next day
c
      ch_ext = extensions(isite)
      if(extensions(isite).eq.'bnd') then
         ch_ext = 'bon'
      else if(extensions(isite).eq.'fpe') then
         ch_ext = 'fpk'
      endif
c
      call read_pres(fnam_date,ch_ext,ista_to_gmt(isite),
     1 istatus_pres)
      if(istatus_pres.ne.0) then
         print *,'pressure could not be accessed for this date,'
         print *,'default to the mean station pressure for Rayleigh'
         print *,'scattering tau'
         n_pres = 1
         p(1) = pressure_nominal(isite)
      endif
c
c
      call i4tim70_atime(i4time_date,a_date,istatus)
c
c pad the ascii date with a 0 in front if necessary
      print *,'the input time and date is ',a_date
      if(a_date(1:1).eq.' ')then
         a_date(1:1)='0'
      endif
c
c
c Start processing
c
c   determine the year of the requested date (i4time70)
      call i4tim70_int(i4time_date,nyear,nmonth,nday,nhour,nmin,nsec,
     1istatus)
      print 10,nmonth,nday,nyear
  10  format(1x,'Processing ',i2.2,'/',i2.2,'/',i2.2)
      if(nyear.gt.93)then
         nyyyy = 1900 + nyear
      else
         nyyyy = 2000 + nyear
      endif
c   convert the date and time integers to character variables
      write(ch_yyyy,15) nyyyy
  15  format(i4.4)
      write(ch_yy,20)nyear
      write(ch_mm,20)nmonth
      write(ch_dd,20)nday
  20  format(i2.2)
c
c
c First process the MFRSR data
c
c   get the MFRSR head serial number, channel central wavelengths, and
c   calibration linear equation coefficients for the chosen date and station
c
      call mfrsr_head(i4time_date,mfr_serial,lstring_serial)
c
c If this MFRSR has a 1625 channel, the following arrays will be used
c to interpolate precipitable water vapor to various times of day
c from soundings analyzed in sub. get_pwv
c
c   convert the sounding i4times to real and load the pv array
c   for use in the interpolation routine
c
      do i = 1,4
         snd_time(i) = float(i4time_pwv(i))
         pv(i) = pwv_snd(i)
      enddo
c
c
c Compute the ozone absorption coefficients for the central wavelengths 
c of the first 5 spectral channels of the MFRSR being used
c
      do i = 1,nchannel
         ozone_coef(i) = trpt(wave,o3coeff,596,wave_l(i),-9.999)
         if(ozone_coef(i).lt.0.0)ozone_coef(i)=0.0
         print *,'ozone absorption coef for ',wave_l(i),' = ',
     1           ozone_coef(i)
      enddo
c
c
c   construct the appropriate mfrsr data filename and open
c
      full_filename_mfr = in_directory_mfr(1:lstring_mfr)//
     1ch_yyyy//'/'//extensions(isite)//ch_yyyy//ch_mm//ch_dd//'_'//
     2mfr_serial(1:lstring_serial)//'.ccc'
      print *,'the mfrsr file is ',full_filename_mfr
      print *,' '
c
      open(unit=lun_mfr,file=full_filename_mfr,status='old',err=200)
c
      call read_mfr(full_filename_mfr,ldate_mfr,ltime_mfr,cosz_mfr,
     1 direct,path_len_mfr,i4time70_local,n_mfr_recs,istatus_read_mfr)
c
      if(istatus_read_mfr.ne.0) then
         print *,'Error reading MFRSR file ',full_filename_mfr
         print *,'...there may be a problem'
         print *,' '
      endif
c
c Skip day if less than 10 MFRSR records read
c
      if(n_mfr_recs.lt.10) then
         close(unit=lun_mfr)
         print *,' '
         print *,'less than 10 MFRSR records found, skip this day'
         go to 350 ! skip day
      endif
      close(unit=lun_mfr)
c
c Enough data, ready the MFRSR data for processing
c
      print *,' '
      print *,n_mfr_recs,' MFRSR records read'
      print *,' '
      n_aod = 0
      do i = 1,n_mfr_recs
c
c  reject all data with path length greater than 7 atmospheres
         if(path_len_mfr(i).gt.7.0)then
           continue! skip
         else ! path length acceptable
           n_aod=n_aod+1
           ltime(n_aod) = ltime_mfr(i)
           i4time70_mfr_local(n_aod) = i4time70_local(i)
           do j=1,nchannel
             if(direct(i,j).le.0.0)then
c     print *,'channel = ',j
c     print *,'ln_direct set to 0 because direct (',direct(i,j),'< 0'
c     print *,' '
               ln_direct(n_aod,j)=-9.999
             else
               ln_direct(n_aod,j)=alog(direct(i,j))
             endif
           enddo
           path(n_aod)=path_len_mfr(i)
         endif
      enddo
c
      go to 300 ! compute aerosol optical depths
c
 200  print *,' '
      print *,'ERROR--Could not open MFRSR data file for ',a_date
      print *,'stop,  filename: ',full_filename_mfr
      close(unit=lun_mfr)
      stop
c
 255  print *,' '
      print *,'ERROR-could not open the ozone.coefs file, stop'
      stop
c
 300  continue 
c
c
c  Get the total ozone for the station and date requested
c
c new total ozone access code (as of 29-apr-2020)
c
         call interp_ozone(ch_ext,total_ozone)
c
         if(total_ozone.lt.0.0) then ! use default o3 of 300
            total_ozone = 300.
         endif
         print *,' '
         print *,'Total ozone used for this day is ',total_ozone
         print *,' '
c
c   Compute Tau_ozone for the five channels
         do j = 1,nchannel
            t_o(j) = ozone_coef(j)*(total_ozone/1000.)
            print *,' '
            print *,'the optical depth due to ozone is ',t_o(j),' for ',
     1              wave_l(j),'nm'
         enddo
c
      p_o = 1013.25 ! define sea level pressure
c
c
c  Compute aerosol optical depths
c
c     loop through all of the accepted mfrsr data
c
      do i = 1,n_aod
c
c  Get the appropriate ln(Vo) for the day being processed
c
d        print *,' '
d        print *,'calling get_izero'
c
         call get_izero(i4time70_mfr_local(i),i)
c
c
c  Compute tau for Rayleigh scattering for the station perssure 
c  corresponding to the current time
c
c     get the station pressure
c
         if(n_pres.le.1)then ! no pressure data, use the default (mean)
c
            p_s = pressure_nominal(isite)
d           print *,'default pressure used : ',p_s
c
         else ! interpolate the station pressure for the current time
c
d           print *,'1st element of pres time and pres arrays = ',
d    1              real_i4tim70_loc(1),p(1)
d           print *,'last element of pres time and pres arrays = ',
d    1              real_i4tim70_loc(n_pres),p(n_pres)
d           print *,'i4time to interp to= ',float(i4time70_mfr_local(i))
c
            p_s = trpt(real_i4tim70_loc,p,n_pres,
     1                 float(i4time70_mfr_local(i)),-9999.9)
c
c    if pressure is missing, set it to the nomimal station pressure
            if(p_s.lt.0.)then
               p_s = pressure_nominal(isite)
            endif
c
d           print *,'interpolated pressure = ',p_s
c
         endif
c
c   load the pressure array that is interpolated to MFRSR times for
c   the output file
c
         pres_int(i) = p_s
c
c  Compute tau Rayleigh scattering (molecular)
c
         if(i.eq.1) print *,' '
c
         do j = 1,nchannel
            exponent = -4.15 + (0.2*(wave_l(j)/1000.))
            exp_term = (wave_l(j)/1000.)**exponent
c
            rayleigh(j) = (0.0088*exp_term*p_s)/p_o
c
            if(i.eq.1) then
              print *,'Tau_Rayleigh for this station is : ',rayleigh(j),
     1            ' for ',wave_l(j),'nm'
            endif
c
         enddo
c
c Compute the optical depths for the 1625 channel as a function of air mass
c and scale co2 and ch4 ODs to the measured station pressure
c
         if(channel_1625.eq..TRUE.)then
            p_scale = p_s/p_o
            od_co2 = (co2int + (co2x*path(i)) + (co2x2*path(i)**2))
     1*p_scale
            od_ch4 = (ch4int + (ch4x*path(i)) + (ch4x2*path(i)**2))
     1*p_scale
            od_h2o_5cm = (h2oint + (h2ox*path(i)) + (h2ox2*path(i)**2))
         endif
c
c  print the ODs for CO2, CH4, and H2O(5cm)--should be 0 if not 1625 nm
         if(i.eq.1) then
            print *,'Tau for CO2 is ',od_co2
            print *,'Tau for CH4 is ',od_ch4
            print *,'Tau for H2O(5cm) for MFRSR ',
     1mfr_serial(1:lstring_serial),' is ',od_h2o_5cm
            print *,' '
         endif
c
c
c  Compute total optical depth for each MFRSR channel
c
         do j = 1,nchannel
c
            if(ln_direct(i,j).gt.0.0)then
c
               tod = (i_zero(j)-ln_direct(i,j))/path(i)
c
c   Subtract effects of molecular scattering and ozone absorption
c   to get aerosol optical depth
c
c   if channel not 3, just subtract Rayleigh and ozone ODs
c
               if (j.ne.3) then
c
                  aod(i,j) = tod - rayleigh(j) - t_o(j)
c
               else ! channel is 3
c
                  if(wave_l(j).lt.1600.) then ! ch 3 is 614nm, not 1625
c
                     aod(i,j) = tod - rayleigh(j) - t_o(j)
c
                  else ! channel 3 is 1625 nm, account for CO2, CH4, and H2O absorption
c
c       first interpolate integrated water vaopr between sounding times to current time
c
                     pv_last = 1.5 ! initialize a default PWV
                     igmt_hhmm = ltime(i) + (ista_to_gmt(isite)*100)
                     ihour_now = igmt_hhmm/100
                     imin_now = mod(igmt_hhmm,100)
                     i_gmt_time_now = i4time_date + ihour_now*3600 +
     1                                imin_now*60
                     current_gmt_time = float(i_gmt_time_now)
c
c              interpolate pwv between 12, 00, and the following 12 GMT to the current time
                     pv_now=trpt(snd_time,pv,4,current_gmt_time,
     1                           -999.0)
c
c              keep updating pv_last for PWV to use if interpolatoin impossible
                     if(pv_now.ne.-999.0)then
                        pv_last = pv_now
                     endif
c
d                    print *,'interpolated pv_now = ', pv_now
c
c               compute optical depth for H2O
c
                     if(pv_now.ne.-999.0)then
c
c
c The H2O od calculation is made for 5 cm integrated water vapor and must 
c be scaled to the actual value (pv_now)
c
c new method scales h2o linearly instead of by the square root
                        od_h2o = od_h2o_5cm*(pv_now/5.0) !scale current PWV from OD at 5cm
c
                     else ! interpolation impossible, use the last PWV value computed
c
c new method scales h2o linearly instead of by the square root
                        od_h2o = od_h2o_5cm*(pv_last/5.0)
c
                     endif
c
c  debug prints for 1625 trace gas optical depths
           write(lun_od,52)ltime(i),path(i),od_co2,od_ch4,od_h2o,
     1rayleigh(j)
  52       format(4x,i4.4,1x,3x,f7.3,4(1x,f7.5))
c
c
                     aod(i,j) = tod - rayleigh(j) - t_o(j) - od_co2 - 
     1                          od_ch4 - od_h2o
d                    print *,'od_h2o = ',od_h2o
                  endif
c
               endif
c
               if(aod(i,j).lt.0.)then
                 print *,'*** negative aod computed for hour',ltime(i)
                 print *,'channel = ',j,'  od_h2o = ',od_h2o
                 print *,'tod, ray, t_o = ',tod,rayleigh(j),t_o(j)
                 print *,'aod = ',aod(i,j)
               endif
c
c       compute the AOD error for the present pathlength 
c
               error_aod(i,j)= error(j)/path(i)
c
            else
c
               aod(i,j)=-9.999
               tod=-9.999
               error_aod(i,j)=-9.9999
c
            endif ! end AOD calculation for one time
c
d           print *,'tod= ',tod,'  Rayleigh = ',rayleigh(j),' TO3 = '
d    1             ,t_o(j)
d           print *,'aod = ',aod(i,j)
c 
         enddo ! end channel (j) loop
c
      enddo ! end AOD computation time (i) loop through the day
c
      close (unit=lun_od)
c--------------- end AOD calculations ----------------------------------
c
c-----------begin pre-specified AOD error deletion section--------------
c
c Now delete AOD data listed in the station error files that have been
c pre-identified as bad
c
c   Open the error listing file for this station
c
      error_file = ch_ext//'_aod_errors.dat'
      open(unit=lun_err,file=in_directory(1:lstring_in)//error_file,
     1status='old',err=21)
c
      print *,' '
      print *,'AOD error file ',in_directory(1:lstring_in)//error_file
      print *,'opened'
      go to 23
c
  21  print *,' '
      print *,'No error file found for ',ch_ext,' continue processing'
      go to 102
c
c    Read the AOD error listing and apply deletions
c
  23  read(lun_err,24,end=102) num_site,ch_begin_date,ch_end_date,
     1ihhmm_start,ihhmm_end,ch_interval,ichannel_start,ichannel_last
  24  format(i2,1x,a5,1x,a5,1x,i4,1x,i4,1x,a1,1x,i1,1x,i1)
c
c the following code does not work with bon and bnd of fpk and fpe
c so it is commented out
c     if(num_site.ne.isite)then
c        print *,' '
c        print *,'Site number in errors file (',num_site,')does not'
c        print *,'match that of the site being analyzed (',isite,')'
c        stop
c     endif
c
      ch_begin_date_full = ch_begin_date//'0000'
      ch_end_date_full = ch_end_date//'0000'
c
      call fnam_i4tim70(ch_begin_date_full,i4time_begin,istatus)
      call fnam_i4tim70(ch_end_date_full,i4time_end,istatus)
c
      if(i4time_date.lt.i4time_begin)go to 23
c
      if(i4time_date.gt.i4time_end) then
         go to 23
      else ! delete specified AOD data
c
c--------------- AOD errors deletion code ----------------------
c
      do i = 1,n_aod
         do j = 1,nchannel
            call i4tim70_int(i4time70_mfr_local(i),idum1,idum2,idum3,
     1                       i_hr,i_min,idum4,istatus)
            ihhmm_test = i_hr*100 + i_min
            if(ch_interval.eq.'b')then ! bad data between two times
               if((ihhmm_test.ge.ihhmm_start).and.
     1            (ihhmm_test.le.ihhmm_end))then
                  if((j.ge.ichannel_start).and.(j.le.ichannel_last))then
                     print *,'deleting specified bad data'
                     aod(i,j) = -9.999
                     error_aod(i,j) = -9.9999
                  endif
               endif
            else if(ch_interval.eq.'o')then ! bad data outside of two times
               if((ihhmm_test.le.ihhmm_start).or.
     1            (ihhmm_test.ge.ihhmm_end))then
                  if((j.ge.ichannel_start).and.(j.le.ichannel_last))then
                     print *,'deleting specified bad data'
                     aod(i,j) = -9.999
                     error_aod(i,j) = -9.9999
                  endif
               endif
            else
               print *,' '
               print *,'Error...at line:'
               print 24,num_site,ch_begin_date,ch_end_date,ihhmm_start,
     1             ihhmm_end,ch_interval,ichannel_start,ichannel_last
               print *,'ch_interval parameter ',ch_interval,' not'
               print *,'recognized, ...stop'
               stop
            endif
         enddo ! end channel (j) loop of deletion code
c
      enddo ! end time (i) loop of deletion code
c
      endif
      go to 23 ! return to read another line from the AOD error file
c
 102  continue ! finished deleting pre-specified bad AOD data
      close (unit=lun_err)
c
c--------------end AOD error deletion section---------------------------
c
c
c
c Finished computing AOD, and deleting bad data, now compute the Angstrom exponent
c
c       compute the alpha denominator (lambda expressed in micrometers)
c
         alpha_denomimator = alog(wave_l(5)/1000.)-alog(wave_l(2)/1000.)
c
c       compute alpha -[dln(aod)/dln(lambda)] based on the 868 and 500 nm
      do i = 1,n_aod
c
         if((aod(i,2).gt.0.0).and.
     1      (aod(i,5).gt.0.0)) then
           alpha_numerator = alog(aod(i,5)) - alog(aod(i,2))
           if(alpha_denomimator.ne.0.)then
              alpha(i) = -(alpha_numerator/alpha_denomimator)
           else
              alpha(i) = -9.999
           endif
         else
            alpha(i) = -9.999
         endif
         if((alpha(i).lt.-9.990).or.(alpha(i).gt.9.990)) then
            alpha(i) = -9.999
         endif
c
      enddo ! end alpha calculation loop for this day
c
c
c-----------------------------------------------------------------------
c
c Now screen AOD data for cloud contamination using the Alexandrov method
c
c   Preliminary cloud screening--no AOD value greater than 3 written to the output file
c
c   870 nm channel used for this screening; if that channel fails, then no
c   output written for that time
c
c   The window size is set based on the frequency of the MFRSR data. The 
c   data frequency is determined by calling subroutine "find_res"
c
c--------------------------------
c      Idenfify 6 & 7-sep-2010 as a special cases (4-mile canyon fire near Tbl Mt.)
c      to have no cloud screening done
c           fnam_tbl_fire_90610 = '102490000'
c           call fnam_i4tim70(fnam_tbl_fire_90610,i4time_90610,istatus)
c           fnam_tbl_fire_90710 = '102500000'
c           call fnam_i4tim70(fnam_tbl_fire_90710,i4time_90710,istatus)
c
c      Identify the date of the calwood fire near Table Mt.
c          fnam_tbl_fire_102120 = '202950000'
c          call fnam_i4tim70(fnam_tbl_fire_102120,i4time_102120,istatus)
c--------------------------------
c
c Increase the acceptable AOD upper limit from 2 to 4 for TBL 6-7 Sept 2010
c
      if(((isite.eq.4).and.(i4time_date.eq.i4time_90610)).or.
     1   ((isite.eq.4).and.(i4time_date.eq.i4time_90710))) then
c     2   ((isite.eq.4).and.(i4time_date.eq.i4time_102120))) then
         aod_limit = 4.0
      else
         aod_limit = 2.0
      endif
c
c override the above limits
c     aod_limit = 5.0
c
      n_clear = 0
      do i = 1,n_aod
         if((aod(i,2).le.aod_limit).and.(aod(i,2).ne.-9.999)) then
            n_clear = n_clear +1
            ltime_clear(n_clear) = ltime(i)
            pres_int_clear(n_clear) = pres_int(i)
            alpha_clear(n_clear) = alpha(i)
            do j = 1,nchannel
               aod_clear(n_clear,j) = aod(i,j)
               error_aod_clear(n_clear,j) = error_aod(i,j)
            enddo
         endif
      enddo
c
      if(n_clear.eq.0) then
         print *,' '
         print *,'No AODs less than ',aod_limit,' found'
         print *,'no output file written for this date'
         stop
      endif
c
c
c---------------------------------------------------------------------
c  Apply the Alexandrov et al. cloud screening method
c
c   Initially mark all times as cloudy
c
      do i = 1,n_clear
         icloud(i) = 1
      enddo
c
c  check if all 870 data are missing
c  if so, set all to good and skip Alexandrov screening
c
      ikount = 0
      do i = 1,n_clear
        if(aod_clear(i,5).ne.-9.999) then 
           ikount = ikount + 1
        endif
      enddo
      if(ikount.eq.0) then ! skip Alexandrov screening
         print *,' '
         print *,'All 870 data missing, Alexandrov screening skipped'
         print *,'Set all data to clear and exit Alexandrov'
         print *,' '
         do i = 1,n_clear
            icloud(i) = 0
         enddo
         go to 600 ! skip the Alexandrov cloud screening
      endif

c
c  find the frequency of the MFRSR data in minutes to set the window size
c
      call find_res(i4time70_local,n_mfr_recs,nres)
      print *,' '
      print *,'Frequency of this days MFRSR data is ',nres,' minutes'
c
c
c Skip cloud screening for Table Mt on 6 & 7-sep-2010, 4-mile canyon fire--clear all day
c
c     if((isite.eq.4).and.(i4time_date.eq.i4time_90610)) then
c
c        do i = 1,n_clear
c           icloud(i) = 0
c        enddo
c        go to 600 ! skip the Alexandrov cloud screening
c
c     endif
c
c     if((isite.eq.4).and.(i4time_date.eq.i4time_90710)) then
c
c        do i = 1,n_clear
c           icloud(i) = 0
c        enddo
c        go to 600 ! skip the Alexandrov cloud screening
c
c     endif
c
c     if((isite.eq.4).and.(i4time_date.eq.i4time_102120)) then
cc
c        do i = 1,n_clear
c           icloud(i) = 0
c        enddo
c        go to 600 ! skip the Alexandrov cloud screening
cc
c     endif
c
c   Initial variables
      loop_stop = 0
      n_screened = 0
c
c
c   Define the cloud-screening window size (n_test), aims at a 15-min window
c 
      n_test = 15/nres
c
c   if the MFRSR data frequency is 5 min or greater, set the test window to 4
      if(nres.ge.5)then
         n_test = 4
      endif
      print *,' '
      print *,'Alexandrov test window set to: ',n_test,' data points'
      print *,' '
c
c  define the date after which the change to 1-min mfrsr data at SRRB sites occurred
      atime_1min = '29-feb-2008 00:00:00.0'
      call atime_i4tim70(atime_1min,i4time_1min,istatus)
c
c
c  The initial cloud screen loop
c
      do k = 1,n_clear
c
         iend = k + (n_test-1) ! define the end of the running window
c
c  check if at the end of the AOD time series, if so, set flag
c  to exit loop later
c
         if(iend.ge.n_clear) then
            iend = n_clear
            loop_stop = 1
         endif
c
c     Average the 870-nm AODs in the current window
c
         tau_total = 0.0
         ncount = 0
         do i = k, iend
            if(aod_clear(i,5).ne.-9.999) then
               tau_total = tau_total + aod_clear(i,5) ! sum the 870nm aod
               ncount = ncount + 1
            endif
         enddo
         if(ncount.gt.0) then
            tau_ave = tau_total/float(ncount)
         else
c
c     870 nm data must be missing, can't do Alexandrov screening
c     therefore set Alexandrov cloud screen index for mid-window 
c     to 9 and escape
c
            index_bad = (k+iend)/2
c           icloud(index_bad) = 9
c ***changed to good on 3-13-2023 to skip alex screening if 870 missing
            icloud(index_bad) = 0
            go to 355 !escape the current running window, leaving
c                      its central point marked as bad (9) 
         endif
c
c     Compute Alexandrovs tau_prime (normalized aods)
c
         do i = k,iend
            if(aod_clear(i,5).ne.-9.999) then
               tau_prime(i) = aod_clear(i,5) - tau_ave + 0.2
            else
               tau_prime(i) = -9.999
            endif
c
c     To avoid problems with taking logs of negative numbers,
c     in Alexandrovs method, set all negative tau_prime times to cloudy
c
            if (tau_prime(i).le.0.0) then
               tau_prime(i) = -9.999
d              print *,'tau_prime negative, set to ',tau_prime(i)
            endif
c
         enddo
c
c     Compute epsilon prime and use it to screen for clouds
c
         tau_total = 0.0
         tau_exp_total = 0.0
         ncount = 0
         do i = k, iend
            if(tau_prime(i).ne.-9.999) then
               tau_total = tau_total + tau_prime(i)
               tau_exp_total = tau_exp_total + alog(tau_prime(i))
               ncount = ncount + 1
            endif
         enddo
c
         if(ncount.le.0) go to 355 ! escape this window, leave center pt. cloudy
c
c   Compute window averages of tau_prime and tau_prime_log
c
         tau_prime_ave = tau_total/float(ncount)
         tau_prime_log_ave = tau_exp_total/float(ncount)
d        print *,'tau_exp_total= ',tau_exp_total,'  tau_prime_ave = ',
d    1tau_prime_ave,'  ncount = ',ncount
c
c     Define the midpoint of the window (index) and mark only it as 
c     cloudy or clear based on the diagnostics
c
         index = (k+iend)/2
         e_prime(index) = 1 - (exp(tau_prime_log_ave)/tau_prime_ave)
d        print *,'e_prime = ',e_prime(index)
c
         if(tau_prime(index).eq.-9.999) then
            go to 355 ! if the window center point has a negative 
c                       tau_prime, leave the point marked as cloudy
         else
            if(e_prime(index).ge.0.0002) then
               icloud(index) = 1
            else
               icloud(index) = 0
               n_screened = n_screened + 1
            endif
         endif
c
c    If this is the end of the time series, then get out
c
 355     if(loop_stop.eq.1) then
            go to 400 ! escape out of cloud screen loop
         endif
c
      enddo ! end Alexandrov initial cloud screen loop
c
c
c Second step in Alexandrov cloud screen--min/max test
c
c   Construct the AOD min/max envelopes for initial cloud-screened AODs
c
 400  loop_stop = 0
      icount = 0
c
c  Determine the min and max AOD within the moving windows of  
c  each point that passed the first test
c
      do k = 1,n_clear
         iend = k + (n_test-1)
         if (iend.ge.n_clear) then
            iend = n_clear
            loop_stop = 1
         endif
c
         aod_min = 10.0
         aod_max = 0.0
c
         do i = k,iend
            if((aod_clear(i,5).ne.-9.999).and.
     1         (icloud(i).eq.0)) then
               aod_max = amax1(aod_max,aod_clear(i,5))
               aod_min = amin1(aod_min,aod_clear(i,5))
            endif
         enddo
c
         index = (k+iend)/2
         if((aod_min.ne.10.0).and.
     1      (aod_max.ne.0.0)) then
            icount = icount + 1
c          Alexandrov scales the mins and maxes by 1.2
            aodmax(icount) = aod_max*1.2
            aodmin(icount) = aod_min/1.2
            ihour = ltime_clear(index)/100
            min = mod(ltime_clear(index),100)
            time(icount) = float(ihour) + float(min)/60.
         endif
         if(loop_stop.eq.1) go to 500
      enddo
c
c
c  Interpolate the AOD mins and maxes to all n_clear points
c
 500  nmax = icount
      if (nmax.le.1) go to 600
c
      do k = 1,n_clear
         ihour = ltime_clear(k)/100
         min = mod(ltime_clear(k),100)
         current_time = float(ihour) + float(min)/60.
         aod_max_interp(k) = trpt(time,aodmax,nmax,current_time,-9.999)
         aod_min_interp(k) = trpt(time,aodmin,nmax,current_time,-9.999)
d        print 501,current_time,aod_min_interp(k),aod_max_interp(k)
d501     format('time = ',f5.2,' aod_min = ',f7.2,'  aod max = ',f7.2)
      enddo
c
c All AODs that did not pass the first test and lie within the
c min/max envelope are now labeled as "cloud-screened" 
c
      do k = 1,n_clear

         if((aod_max_interp(k).ne.-999.).and.
     1      (aod_min_interp(k).ne.-999.)) then
            if((icloud(k).eq.1).and.
     1         (aod_clear(k,5).lt.aod_max_interp(k)).and.
     2         (aod_clear(k,5).gt.aod_min_interp(k))) then 
               icloud(k) = 0
               n_screened = n_screened + 1 !add to the total screened
            endif
         endif
      enddo
c
 600  continue ! escape point for min/max interpolation code
c
c End Alexandrov cloud screening
c
c
c  Compute daily averages of cloud screened AOD for each channel
c
d     print *,' '
d     print *,'Before averaging alexandrov n_clear = ',n_clear
d     print *,' '
      aod_ave_415 = 0.0
      aod_ave_500 = 0.0
      aod_ave_614 = 0.0
      aod_ave_670 = 0.0
      aod_ave_868 = 0.0
      icount_415 = 0
      icount_500 = 0
      icount_614 = 0
      icount_670 = 0
      icount_868 = 0
      do l = 1,n_clear
d        print *,'icloud(l) = ',icloud(l)
         if(icloud(l).eq.0) then
            if(aod_clear(l,1).gt.0.0)then
               aod_ave_415 = aod_ave_415 + aod_clear(l,1)
               icount_415 = icount_415 + 1
            endif
d           print *,'clear aod 500 = ',aod_clear(l,2)
            if(aod_clear(l,2).gt.0.0)then
               aod_ave_500 = aod_ave_500 + aod_clear(l,2)
               icount_500 = icount_500 + 1
            endif
            if(aod_clear(l,3).gt.0.0)then
               aod_ave_614 = aod_ave_614 + aod_clear(l,3)
               icount_614 = icount_614 + 1
            endif
            if(aod_clear(l,4).gt.0.0)then
               aod_ave_670 = aod_ave_670 + aod_clear(l,4)
               icount_670 = icount_670 + 1
            endif
            if(aod_clear(l,5).gt.0.0)then
               aod_ave_868 = aod_ave_868 + aod_clear(l,5)
               icount_868 = icount_868 + 1
            endif
         endif
      enddo
d     print *,' '
d     print *,'after tally, icount_500 = ',icount_500
d     print *,'aod_ave_500 total = ',aod_ave_500
d     print *,' '
c
c Special case for Table Mt. 7 Sept. 2010
c
      if((i4time_date.eq.i4time_90710).and.(isite.eq.4)) then
c
      aod_ave_415 = 0.0
      aod_ave_500 = 0.0
      aod_ave_614 = 0.0
      aod_ave_670 = 0.0
      aod_ave_868 = 0.0
      icount_415 = 0
      icount_500 = 0
      icount_614 = 0
      icount_670 = 0
      icount_868 = 0
      do l = 1,n_clear
         if((icloud(l).eq.0).and.(ltime_clear(l).lt.1518)) then
            if(aod_clear(l,1).gt.0.0)then
               aod_ave_415 = aod_ave_415 + aod_clear(l,1)
               icount_415 = icount_415 + 1
            endif
            if(aod_clear(l,2).gt.0.0)then
               aod_ave_500 = aod_ave_500 + aod_clear(l,2)
               icount_500 = icount_500 + 1
            endif
            if(aod_clear(l,3).gt.0.0)then
               aod_ave_614 = aod_ave_614 + aod_clear(l,3)
               icount_614 = icount_614 + 1
            endif
            if(aod_clear(l,4).gt.0.0)then
               aod_ave_670 = aod_ave_670 + aod_clear(l,4)
               icount_670 = icount_670 + 1
            endif
            if(aod_clear(l,5).gt.0.0)then
               aod_ave_868 = aod_ave_868 + aod_clear(l,5)
               icount_868 = icount_868 + 1
            endif
         else
            icloud(l) = 1 ! after 1518 clouds moved in
         endif
      enddo
      endif
c 
c
c    set the daily average AOD missing value
      aod_miss = -9.999
c
      if(icount_415.gt.0)then
         aod_ave_415 = aod_ave_415/float(icount_415)
      else
         aod_ave_415 = aod_miss
      endif
      if(icount_500.gt.0)then
         aod_ave_500 = aod_ave_500/float(icount_500)
      else
         aod_ave_500 = aod_miss
      endif
      if(icount_614.gt.0)then
         aod_ave_614 = aod_ave_614/float(icount_614)
      else
         aod_ave_614 = aod_miss
      endif
      if(icount_670.gt.0)then
         aod_ave_670 = aod_ave_670/float(icount_670)
      else
         aod_ave_670 = aod_miss
      endif
      if(icount_868.gt.0)then
         aod_ave_868 = aod_ave_868/float(icount_868)
      else
         aod_ave_868 = aod_miss
      endif
c
c
c Do not write an output file if no AOD could be computed
c
      if(aod_ave_500.lt.-9.990)then
        print 205
 205    format(/,"No viable AOD computed, skip writing the output file")
        print *,' '
        print *,'aod_ave_500 (post alexandrov) = ',aod_ave_500
        print *,'icount_500 = ',icount_500
        go to 305 ! no aod computed, skip writing an output file
      endif
c
c
c special case for Table Mt on 6 & 7-sep-2010, 4-mile canyon fire--clear all day
c
c     if((isite.eq.4).and.(i4time_date.eq.i4time_90610)) then
c        n_screened = n_clear
c     endif
c
c     if((isite.eq.4).and.(i4time_date.eq.i4time_90710)) then
c        n_screened = icount_500
c     endif
c
c  Write the results to a file
c
c     Construct the output filename
c
      aod_filename = out_directory(1:lstring_out)//
     1ch_ext//'_'//ch_yyyy//ch_mm//ch_dd//'.aod'
c
      print *,' '
      print *,'the output AOD file name is ',aod_filename
c
c  Open the output file
c
      open(unit=lun_aod,file=aod_filename,status='unknown',err=250)
      print *,' '
      print *,aod_filename,' opened'
c
c  Construct  and write the header records
c
      write(lun_aod,210) names(isite)
 210  format(a<lstring(isite)>,1x,'SURFRAD aerosol optical depth')
c
c
c set the daily average AOD to missing 
      aod_miss = -9.999
c
      write(lun_aod,225) a_date(1:11),nday,nmonth,nyyyy,n_clear
 225  format(a11,1x,i2.2,1x,i2.2,1x,i4.4,1x,i3,1x,'lines of data')
c
      write(lun_aod,227) (wave_l(k),k=1,nchannel)
 227  format(<nchannel>(f6.1,1x),' channel central wavelengths (nm)')
c
      write(lun_aod,229) aod_ave_415,aod_ave_500,aod_ave_614,
     1 aod_ave_670,aod_ave_868,n_screened
 229  format(5(f6.3,1x),'Daily average AODs -- ',i3,' sample size')
c
      write(lun_aod,228) int(total_ozone)
 228  format(i3.3,1x,'Dobson units of ozone')
c
c  Construct the column header record
c
c     Convert wavelength channels to character variables
c
      do i = 1,nchannel
         i_wave = nint(wave_l(i))
         if(i_wave.le.999)then
            write(ch_iwave,240)i_wave
 240        format(i3.3)
            ch_wave_l(i) = 'OD'//ch_iwave
            ch_error(i) = ch_iwave//'E'
         else ! accomodate wavelengths beyond 999 nm
            write(ch_iwave4,241)i_wave
 241        format(i4.4)
            ch_wave_l(i) = 'OD'//ch_iwave4
            ch_error(i) = ch_iwave4//'E'
         endif
      enddo
c
c  Write the column header record
c
      write(lun_aod,242) (ch_wave_l(k),k=1,nchannel),
     1                   (ch_error(k),k=1,nchannel)
 242  format('ltime  0=good  ',<nchannel>(a,1x),<nchannel-1>(a,3x),
     1a4,3x,'p_mb',2x,'Ang_exp')
c
c  Write the data
c
      do i = 1,n_clear
         write(lun_aod,230) ltime_clear(i),icloud(i),
     1     (aod_clear(i,j),j=1,nchannel),
     2     (error_aod_clear(i,j),j=1,nchannel),pres_int_clear(i),
     3     alpha_clear(i)
 230     format(i4.4,5x,i1,3x,<nchannel>(1x,f6.3),<nchannel>(1x,f7.4),
     1          1x,f6.1,1x,f6.3)
d       print 230,ltime_clear(i),icloud(i),
d    1     (aod_clear(i,j),j=1,nchannel),
d    2     (error_aod_clear(i,j),j=1,nchannel),pres_int_clear(i),
d    3     alpha_clear(i)
      enddo
c
      close(unit=lun_aod)
c
      stop !Normal Stop
c
c
 250  print *,' '
      print *,'ERROR--Could not open a new output file'
      stop
c
 305  print *,' '
      print *,'No clear points found, skip output file'
      stop
c
 350  print *,' '
      print *,'only ',n_mfr_recs,' in input file'
      print *,'no output AOD file created, skip this day'
      stop
c
      end
c-----------------------------------------------------------------------
      subroutine coeff_a(index)
c
c this subroutine contains absorption fit coefficients as a function 
c of air mass for ch4, co2 and h20 for the "a" batch of 1625 filters
c
      common/coefficients/ch4int,ch4x,ch4x2,co2int,co2x,co2x2,
     1h2oint,h2ox,h2ox2
c
      real*4 ch4_int(4),ch4_x(4),ch4_x2(4)
      real*4 co2_int(4),co2_x(4),co2_x2(4)
      real*4 h2o_int(4),h2o_x(4),h2o_x2(4)
c
      data (ch4_int(i),i=1,4)/0.00228,
     c                        0.00184,
     c                        0.00200,
     c                        0.00211/
      data (ch4_x(i),i=1,4)/-0.00006,
     c                      -0.00004,
     c                      -0.00005,
     c                      -0.00006/
      data (ch4_x2(i),i=1,4)/0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000/
c
      data (co2_int(i),i=1,4)/0.00508,
     c                        0.00693,
     c                        0.00601,
     c                        0.00572/
      data (co2_x(i),i=1,4)/-0.00044,
     c                      -0.00065,
     c                      -0.00054,
     c                      -0.00051/
      data (co2_x2(i),i=1,4)/0.00003,
     c                       0.00005,
     c                       0.00004,
     c                       0.00004/
c
      data (h2o_int(i),i=1,4)/0.00424,
     c                        0.00445,
     c                        0.00433,
     c                        0.00432/
      data (h2o_x(i),i=1,4)/-0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011/
      data (h2o_x2(i),i=1,4)/0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000/
c
      ch4int = ch4_int(index)
      ch4x= ch4_x(index)
      ch4x2 = ch4_x2(index)
c
      co2int = co2_int(index)
      co2x= co2_x(index)
      co2x2 = co2_x2(index)
c
      h2oint = h2o_int(index)
      h2ox= h2o_x(index)
      h2ox2 = h2o_x2(index)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine coeff_b(index)
c
c this subroutine contains absorption fit coefficients as a function
c of air mass for ch4, co2 and h20 for the "b" batch of 1625 filters
c
      common/coefficients/ch4int,ch4x,ch4x2,co2int,co2x,co2x2,
     1h2oint,h2ox,h2ox2
c
      real*4 ch4_int(10),ch4_x(10),ch4_x2(10)
      real*4 co2_int(10),co2_x(10),co2_x2(10)
      real*4 h2o_int(10),h2o_x(10),h2o_x2(10)
c
      data (ch4_int(i),i=1,10)/0.00228,
     c                        0.00221,
     c                        0.00229,
     c                        0.00214,
     c                        0.00235,
     c                        0.00211,
     c                        0.00232,
     c                        0.00229,
     c                        0.00240,
     c                        0.00219/
      data (ch4_x(i),i=1,10)/-0.00006,
     c                      -0.00006,
     c                      -0.00006,
     c                      -0.00005,
     c                      -0.00007,
     c                      -0.00006,
     c                      -0.00006,
     c                      -0.00006,
     c                      -0.00007,
     c                      -0.00005/
      data (ch4_x2(i),i=1,10)/0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000/
c
      data (co2_int(i),i=1,10)/0.00527,
     c                        0.00523,
     c                        0.00516,
     c                        0.00544,
     c                        0.00494,
     c                        0.00566,
     c                        0.00516,
     c                        0.00504,
     c                        0.00486,
     c                        0.00527/
      data (co2_x(i),i=1,10)/-0.00046,
     c                      -0.00045,
     c                      -0.00044,
     c                      -0.00047,
     c                      -0.00042,
     c                      -0.00050,
     c                      -0.00045,
     c                      -0.00043,
     c                      -0.00041,
     c                      -0.00045/
      data (co2_x2(i),i=1,10)/0.00003,
     c                       0.00003,
     c                       0.00003,
     c                       0.00003,
     c                       0.00003,
     c                       0.00003,
     c                       0.00003,
     c                       0.00003,
     c                       0.00003,
     c                       0.00003/
c
      data (h2o_int(i),i=1,10)/0.00429,
     c                        0.00426,
     c                        0.00428,
     c                        0.00427,
     c                        0.00424,
     c                        0.00432,
     c                        0.00429,
     c                        0.00424,
     c                        0.00424,
     c                        0.00426/
      data (h2o_x(i),i=1,10)/-0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011/
      data (h2o_x2(i),i=1,10)/0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000/
c
      ch4int = ch4_int(index)
      ch4x= ch4_x(index)
      ch4x2 = ch4_x2(index)
c
      co2int = co2_int(index)
      co2x= co2_x(index)
      co2x2 = co2_x2(index)
c
      h2oint = h2o_int(index)
      h2ox= h2o_x(index)
      h2ox2 = h2o_x2(index)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine coeff_c(index)
c
c this subroutine contains absorption fit coefficients as a function
c of air mass for ch4, co2 and h20 for the "c" batch of 1625 filters
c
      common/coefficients/ch4int,ch4x,ch4x2,co2int,co2x,co2x2,
     1h2oint,h2ox,h2ox2
c
      real*4 ch4_int(6),ch4_x(6),ch4_x2(6)
      real*4 co2_int(6),co2_x(6),co2_x2(6)
      real*4 h2o_int(6),h2o_x(6),h2o_x2(6)
c
      data (ch4_int(i),i=1,6)/0.00278,
     c                        0.00338,
     c                        0.00321,
     c                        0.00315,
     c                        0.00321,
     c                        0.00322/
      data (ch4_x(i),i=1,6)/-0.00009,
     c                      -0.00013,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00012/
      data (ch4_x2(i),i=1,6)/0.00000,
     c                       0.00001,
     c                       0.00001,
     c                       0.00001,
     c                       0.00001,
     c                       0.00001/
c
      data (co2_int(i),i=1,6)/0.00505,
     c                        0.00396,
     c                        0.00419,
     c                        0.00429,
     c                        0.00419,
     c                        0.00418/
      data (co2_x(i),i=1,6)/-0.00045,
     c                      -0.00033,
     c                      -0.00035,
     c                      -0.00036,
     c                      -0.00035,
     c                      -0.00035/
      data (co2_x2(i),i=1,6)/0.00003,
     c                       0.00002,
     c                       0.00003,
     c                       0.00003,
     c                       0.00002,
     c                       0.00002/
c
      data (h2o_int(i),i=1,6)/0.00450,
     c                        0.00458,
     c                        0.00455,
     c                        0.00452,
     c                        0.00455,
     c                        0.00453/
      data (h2o_x(i),i=1,6)/-0.00012,
     c                      -0.00014,
     c                      -0.00013,
     c                      -0.00013,
     c                      -0.00013,
     c                      -0.00013/
      data (h2o_x2(i),i=1,6)/0.00000,
     c                       0.00001,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000/
c
      ch4int = ch4_int(index)
      ch4x= ch4_x(index)
      ch4x2 = ch4_x2(index)
c
      co2int = co2_int(index)
      co2x= co2_x(index)
      co2x2 = co2_x2(index)
c
      h2oint = h2o_int(index)
      h2ox= h2o_x(index)
      h2ox2 = h2o_x2(index)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine coeff_d(index)
c
c this subroutine contains absorption fit coefficients as a function
c of air mass for ch4, co2 and h20 for the "d" batch of 1625 filters
c
      common/coefficients/ch4int,ch4x,ch4x2,co2int,co2x,co2x2,
     1h2oint,h2ox,h2ox2
c
      real*4 ch4_int(12),ch4_x(12),ch4_x2(12)
      real*4 co2_int(12),co2_x(12),co2_x2(12)
      real*4 h2o_int(12),h2o_x(12),h2o_x2(12)
c
      data (ch4_int(i),i=1,12)/0.00206,
     c                        0.00190,
     c                        0.00181,
     c                        0.00191,
     c                        0.00191,
     c                        0.00181,
     c                        0.00199,
     c                        0.00186,
     c                        0.00200,
     c                        0.00204,
     c                        0.00206,
     c                        0.00211/
      data (ch4_x(i),i=1,12)/-0.00006,
     c                      -0.00004,
     c                      -0.00004,
     c                      -0.00005,
     c                      -0.00005,
     c                      -0.00004,
     c                      -0.00005,
     c                      -0.00005,
     c                      -0.00005,
     c                      -0.00005,
     c                      -0.00006,
     c                      -0.00005/
      data (ch4_x2(i),i=1,12)/0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000/
c
      data (co2_int(i),i=1,12)/0.00575,
     c                        0.00659,
     c                        0.00689,
     c                        0.00704,
     c                        0.00625,
     c                        0.00771,
     c                        0.00594,
     c                        0.00743,
     c                        0.00598,
     c                        0.00589,
     c                        0.00581,
     c                        0.00560/
      data (co2_x(i),i=1,12)/-0.00051,
     c                      -0.00061,
     c                      -0.00064,
     c                      -0.00067,
     c                      -0.00056,
     c                      -0.00076,
     c                      -0.00053,
     c                      -0.00072,
     c                      -0.00054,
     c                      -0.00052,
     c                      -0.00052,
     c                      -0.00049/
      data (co2_x2(i),i=1,12)/0.00004,
     c                       0.00004,
     c                       0.00005,
     c                       0.00005,
     c                       0.00004,
     c                       0.00005,
     c                       0.00004,
     c                       0.00005,
     c                       0.00004,
     c                       0.00004,
     c                       0.00004,
     c                       0.00003/
c
      data (h2o_int(i),i=1,12)/0.00434,
     c                        0.00443,
     c                        0.00447,
     c                        0.00448,
     c                        0.00440,
     c                        0.00454,
     c                        0.00433,
     c                        0.00451,
     c                        0.00436,
     c                        0.00435,
     c                        0.00434,
     c                        0.00432/
      data (h2o_x(i),i=1,12)/-0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00012,
     c                      -0.00011,
     c                      -0.00012,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011,
     c                      -0.00011/
      data (h2o_x2(i),i=1,12)/0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000,
     c                       0.00000/
c
      ch4int = ch4_int(index)
      ch4x= ch4_x(index)
      ch4x2 = ch4_x2(index)
c
      co2int = co2_int(index)
      co2x= co2_x(index)
      co2x2 = co2_x2(index)
c
      h2oint = h2o_int(index)
      h2ox= h2o_x(index)
      h2ox2 = h2o_x2(index)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine find_res(i4time70_local,n_mfr_recs,nres)
c
c This routine determines the frequency of the MFRSR data for the day
c being processed in minutes. Output is "nres," the resolution in min.
c
      parameter(mfrdata=1500)
c
      integer*4 ifrequency(15)
      integer*4 i4time70_local(mfrdata)
c
c test for 1 through 15-min. frequencies
      do i = 1,15
         ifrequency(i) = 0
      enddo
c
c construct a histogram of frequencies
      do i = 1,(n_mfr_recs-1)
         dif = float((i4time70_local(i+1))-(i4time70_local(i)))/60.
         idif = nint(dif)
d        print *,'idif = ',idif,'minutes'
         if((idif.gt.0).and.(idif.le.15)) then
            ifrequency(idif) = ifrequency(idif) + 1
         endif
      enddo
c
d     print *,' '
d     do i = 1,15
d        print *,'freq ',i,' minutes occurred ',ifrequency(i),' times'
d     enddo
c
c determine the most frequent frequency in the histogram
      most_frequent = 0
      do i=1,15
         if(ifrequency(i).gt.most_frequent) then
            index = i
            most_frequent = ifrequency(i)
         endif
      enddo
c
c default to 1-minute frequency if the resolution could not be
c determined
      if(most_frequent.lt.1) then
         index=1
         print *,' '
         print *,'***Could not determine the temporal resolution'
         print *,'   defaulted to 1-minute'
         print *,' '
      else
         print *,' ' 
         print *,'Most frequent resolution found is ',index
         print *,'minutes'
         print *,' '
      endif
c
c compute the temporal resolution in minutes
      nres = index
c
      return
      end
c----------------------------------------------------------------------
        subroutine get_izero(i4time_local,icount)
c
c This subroutine restores the normalized i-zero value, and corrects it
c for the earth-sun distance,that is appropriate to the day of year being 
c processed.
c
c
      parameter (nchannel=5)
c
      common/calibration/i_zero(nchannel),error(nchannel),
     1  constant_coef(nchannel),slope_coef(nchannel),
     2  sine_coef(nchannel),cos_coef(nchannel),
     3  err_slope_i_zero(nchannel),err_intercept_i_zero(nchannel)
c
      common/unit/lun_mfr,lun_mfrhead,lun_aod,lun_o3,lun_p,lun_cof,
     1lun_hom,lun_err,lun_snd
c
      integer*4 njul_days(12)
      integer*4 date_fmt_changed
      integer*4 i4time_local
c
      real*4 i_zero
      real*4 intercept_i_zero
      real*4 i_zero_corr
      real*4 i_zero_norm(nchannel)
c
      data njul_days/0,31,59,90,120,151,181,212,243,273,304,334/
c
d     if (icount.eq.1)then
d        print *,'...in get_i_zero'
d        do k = 1,nchannel
d           print *,' '
d           print *,'constant_coef(',k,') = ',constant_coef(k)
d           print *,'slope_coef(',k,') = ',slope_coef(k)
d           print *,'sine_coef(',k,') = ',sine_coef(k)
d           print *,'cos_coef(',k,') = ',cos_coef(k)
d        enddo
d     endif
c
      pi = 3.1415927
      two_pi = 2.*pi
c
c First compute the day of the year
      call i4tim70_int(i4time_local,nyear,nmonth,nday,nhour,nmin,
     1                         nsec,istatus)
d     if (icount.eq.1)then
d     print *,'in get_izero'
d     print *,'date/time from i4tim70_int is'
d     print *,nyear,nmonth,nday,nhour,nmin,nsec
d     endif
c
      njulian = njul_days(nmonth) + nday
c
c Now convert the date to yyyy.ffff format
      if(nyear.le.93) then
         nyyyy = 2000 + nyear
      else
         nyyyy = 1900 + nyear
      endif
c
      if(mod(nyyyy,4).eq.0) then
         fract_year = float(njulian)/366.
      else
         fract_year = float(njulian)/365.
      endif
c
      yyyy_frac = float(nyyyy)+fract_year
      if (icount.eq.1)then
         print *,' '
         print *,'fractional date = ',yyyy_frac
      endif
c
      angular_date = yyyy_frac*two_pi

c
c Compute a normalized i-zero and error using the calibration history 
c for each channel (normalized means appropriate to a circular orbit)
c
      do k = 1,nchannel
         if((sine_coef(k).eq.0.).and.(cos_coef(k).eq.0.)) then
            if (icount.eq.1)then
               print *,'Linear fit to i_zero used'
            endif
            i_zero_norm(k) = constant_coef(k) + slope_coef(k)*yyyy_frac
         else
            if (icount.eq.1)then
               print *,'cyclic fit to i_zero used'
            endif
            i_zero_norm(k) = constant_coef(k) +
     1                    slope_coef(k)*angular_date +
     2                    sine_coef(k)*sin(angular_date) +
     3                    cos_coef(k)*cos(angular_date)
         endif
c
         if (icount.eq.1)then
         print *,'normalized i_zero(',k,') = ',i_zero_norm(k)
         endif
         error(k) = err_slope_i_zero(k)*yyyy_frac+
     1              err_intercept_i_zero(k)
      end do
c
d     if (icount.eq.1)then
d        do k = 1,nchannel
d           print *,' '
d           print *,'angular date = ',angular_date
d           print *,'constant_coef(',k,') = ',constant_coef(k)
d           print *,'slope_coef(',k,')*angular_date = ',slope_coef(k)*
d    1angular_date
d           print *,'sine_coef(',k,')*sin = ',sine_coef(k)*
d    1sin(angular_date)
d           print *,'cos_coef(',k,')*cos = ',cos_coef(k)*
d    1cos(angular_date)
d           print *,'normalized i_zero(',k,') = ',i_zero_norm(k)
d           print *,'error in i_zero(',k,') = ',error(k)
d        enddo
d     endif
c
c
c Compute the correction factor for the correct earth-sun distance for 
c the day being processed
c
      doy = float(njulian - 1)! day-1 used in gamma calculation below 
c
c Compute the correction factor for the current earth-sun distance
c
      gamma = (2.*pi*doy)/365.
c E-zero is the correction factor to normalize the I zero values to one
c astronomical unit
      ezero = 1.000110+(0.034221*(COS(gamma)))+
     1           (0.001280*(SIN(gamma)))+(0.000719*(COS(2.*gamma)))
     2          +(0.000077*(SIN(2.*gamma)))
      if(icount.eq.1) print *,'AU correction factor, ezero = ',ezero
c
c correct the normalized i_zero to the current earth-sun distance for all 5
c channels
      do k = 1,nchannel
         i_zero_corr = ezero*i_zero_norm(k)
c
         if(icount.eq.1)print *,'i_zero_corr = ',i_zero_corr
c
c compute the log of the i_zero 
         i_zero(k) = alog(i_zero_corr)
c
         if(icount.eq.1)print *,'log i_zero_corr = ',i_zero(k)
      end do
c
c
      return
      end
c-----------------------------------------------------------------------
      subroutine interp_ozone(ch_ext,total_ozone)
c
c This subroutine attempts ito acquire total ozone for the day
c being analyzed
c
      parameter (nsta=10)
      parameter (nchannel=5)
      parameter (mfrdata=1500)
      parameter (nvals=5000)
      parameter (nbaro=3000)
      parameter (no3=15000)
c
      common/unit/lun_mfr,lun_mfrhead,lun_aod,lun_o3,lun_p,lun_cof,
     1lun_hom,lun_err,lun_snd
      common/directory/in_directory,in_directory_mfr,out_directory,
     1 in_directory_pres,in_directory_o3,home_directory,snd_directory,
     2 lstring_in,lstring_mfr,lstring_out,lstring_p,lstring_o3,len_hom,
     3 len_snd
      common/files/isite,extensions(nsta),names(nsta),lstring(nsta),
     1i4time_date,ista_to_gmt(nsta),pressure_nominal(nsta),xlat,xlon
c
      character*24 names
      character*3 extensions
      character*3 sta
      character*3 ch_ext
      character*80 in_directory,in_directory_mfr,out_directory
      character*80 in_directory_pres,in_directory_o3
      character*80 home_directory
      character*80 snd_directory
      character*80 ozone_file
c
      real*4 x4time_o3(no3),t_o3(no3)
c
c  define the station-specific total ozone file name
c    e.g., /home/surfrad/aod/ozone/tbl_ozone.dat
c
      ozone_file = in_directory_o3(1:lstring_o3)//
     1 ch_ext//'_ozone.dat'
c
c  Open the total ozone file for this station
c
      open(unit=lun_o3,file=ozone_file,status='old',
     1     err=50)
      icount = 0
      do i =1,no3
         read(lun_o3,*,end=10) iyyyy,mm,idd,i4time_line,i_ozone
c  23     format(i4,', ',i2,', ',i2,', ',i11,', ',i3.3)
         icount = icount + 1
         x4time_o3(icount) = float(i4time_line)
         t_o3(icount) = float(i_ozone)
      enddo
      print *,' '
      print *,'ERROR--too many lines in ',ozone_file
      stop
c
  10  print *,' '
      print *,icount,' lines read from ',ozone_file
      close (unit=lun_o3)
c
c Interpolate total ozone for the day being analyzed
c
      call trpt_o3(x4time_o3,t_o3,icount,float(i4time_date),
     1iflag,index_i,index_i_plus_1,total_ozone)
c
      if(iflag.eq.0) then !interpolation not possible
         total_ozone = -9999.9
         print *,' '
         print *,'Total ozone interpolation not possible'
         print *,'Default value of 300 will be used'
      else
         print *,' '
         print *,'Interpolated total ozone for this day is ',total_ozone
      endif
c
      return
c
  50  print *,' '
      print *,'Could not open ', ch_ext//'_ozone.dat'
      stop
      end
c-----------------------------------------------------------------------
      subroutine trpt_o3(p,a,n,zen,int_flg,index_i,index_i_plus_1,
     1o3_interp)
c
c     linear interpolation of total ozone between two days
c        returns -9999.9 if interpolation not possible
c
c     p- real i4time array for days
c     a- array of daily total ozone corresponding to the i4time array
c     n- size of the input arrays p and a
c     zen real i4time to interpolate to
c     o3_interp interpolated o3 value
c     int_flg - interpolation flag
c             =0, then interpolation not possible
c             =1, then exact match, interpolation not done 
c             =2, then interpolation done between index_i
c                                      and index_i_plus_1
c
      dimension p(n),a(n)
c
      int_flg = 0 ! assume interpolation not possible
      index_i = 0
      index_i_plus_1 = 0
c
      if(zen.gt.p(n))go to 5
c
c
      nm1=n-1
      do 10 i=1,nm1
      if(zen.lt.p(i)) go to 5
      if(zen.gt.p(i+1)) go to 10
c   found the time interval to interpolate across
c      first check if one of the two is the exact i4time needed
        if(zen.eq.p(i))then
           o3_interp = a(i)
           int_flg = 1 !exact match to desired i4time
           return
        else if(zen.eq.p(i+1)) then
           o3_interp = a(i+1)
           int_flg = 1 !exact match to desired i4time
           return
        else
           continue
        endif
c
c   interploation possible between airmasses p(i) and p(i+1)
         x1=p(i)
         x2=p(i+1)
         x=zen
         o3_interp=(a(i)*(x-x2)-a(i+1)*(x-x1))/(x1-x2)
         int_flg = 2 !interpolated to desired airmass
         index_i = i
         index_i_plus_1 = i+1
         return
   10 continue
c
    5 o3_interp = -9999.9
      int_flg = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine average_o3(ave_o3,sta)
c
c This subroutine attempts to compute the average total ozone of the
c days before and after a day when the TOMS total ozone for this site
c is missing
c
      parameter (nsta=10)
      parameter (nchannel=5)
      parameter (mfrdata=1500)
      parameter (nvals=5000)
      parameter (nbaro=3000)
c
      common/unit/lun_mfr,lun_mfrhead,lun_aod,lun_o3,lun_p,lun_cof,
     1lun_hom,lun_err,lun_snd
      common/directory/in_directory,in_directory_mfr,out_directory,
     1 in_directory_pres,in_directory_o3,home_directory,snd_directory,
     2 lstring_in,lstring_mfr,lstring_out,lstring_p,lstring_o3,len_hom,
     3 len_snd
      common/files/isite,extensions(nsta),names(nsta),lstring(nsta),
     1i4time_date,ista_to_gmt(nsta),pressure_nominal(nsta),xlat,xlon
c
      character*24 names
      character*24 atime_before,atime_after
      character*3 extensions
      character*3 sta
      character*80 in_directory,in_directory_mfr,out_directory
      character*80 in_directory_pres,in_directory_o3
      character*80 home_directory
      character*80 ozone_file_before
      character*80 ozone_file_after
      character*80 snd_directory
c
c define the day before amd after in i4time, base 1970
c
      iday_before = i4time_date - 86400
      iday_after = i4time_date + 86400
c
      call i4tim70_atime(iday_before,atime_before,istatus)
      if(atime_before(1:1).eq.' ')atime_before(1:1)='0'
c
c  Construct the ozone file name for the day before
c    e.g., /home/surfrad/aod/ozone/tbl/tbl_12-jan-2001.o3
c
      ozone_file_before = in_directory_o3(1:lstring_o3)//
     1 sta//'/'//sta//'_'//atime_before(1:11)//'.o3'
c
c  Open the total ozone file from the day before
c
      open(unit=lun_o3,file=ozone_file_before,status='unknown',
     1     err=50)
      read(lun_o3,*,end=50)total_ozone
      close (unit=lun_o3)
c
      if(total_ozone.lt.0.)then
         go to 50
      else
         total_ozone_before = total_ozone
      endif
c
c Now try the day after
      call i4tim70_atime(iday_after,atime_after,istatus)
      if(atime_after(1:1).eq.' ')atime_after(1:1)='0'
c
c  Construct the ozone file name
c    e.g., /home/surfrad/aod/ozone/tbl/tbl_12-jan-2001.o3
c
      ozone_file_after = in_directory_o3(1:lstring_o3)//
     1 sta//'/'//sta//'_'//atime_after(1:11)//'.o3'
c 
c  Open the total ozone file from the day after
c 
      open(unit=lun_o3,file=ozone_file_after,status='unknown',
     1     err=50)
      read(lun_o3,*,end=50)total_ozone
      close (unit=lun_o3)
c
      if(total_ozone.lt.0.)then
         go to 50
      else
         total_ozone_after = total_ozone
      endif
c
c Successful accessing the total ozone for the day before and after
c Now compute the average
c
      ave_o3 = (total_ozone_before + total_ozone_after)/2.
c
      return ! normal exit
c
  50  ave_o3 = -9999.9
      close (unit=lun_o3)
      return
      end
c-----------------------------------------------------------------------
         subroutine read_mfr(filename,ldate_mfr,ltime_mfr,cosz_mfr,
     1     direct,path_len_mfr,i4time70_local,n_mfr_recs,
     2     istatus_read_mfr)
c
c This routine reads MFRSR data from the file filename that is
c passed as an input augument.  These files are organized in gmt
c and the times must be corrected to local time in order to line 
c up with the clear-id files. Also, the time in these ccc files 
c (in the first column) is in decimal days since Jan. 1, 1900.  
c To avoid dealing with that, the fractional day is stripped off
c to calculate the time (gmt), and the year, month, and day are 
c input to the subroutine via input arguments.  They were used to 
c construct the ccc filename.  
c
        parameter (nsta=10,mfrdata=1500)
        parameter (nchannel=5)
c
      common/unit/lun_mfr,lun_mfrhead,lun_aod,lun_o3,lun_p,lun_cof,
     1lun_hom,lun_err,lun_snd
      common/files/isite,extensions(nsta),names(nsta),lstring(nsta),
     1i4time_date,ista_to_gmt(nsta),pressure_nominal(nsta),xlat,xlon
c
        character*24 names
        character*3 extensions
        character*24 atime_mfr
        character*80 filename
c
        character*9 ch_start_date_bon
        character*9 ch_end_date_bon
        character*24 atime_gmt
        character*9 ch_yydddhhmm
c
        integer*4 ista_to_gmt
        integer*4 ldate_mfr(mfrdata),ltime_mfr(mfrdata)
        integer*4 i4time70_local(mfrdata)
        integer*4 i4time_date
c
        real*4 direct(mfrdata,nchannel)
        real*4 cosz_mfr(mfrdata)
        real*4 path_len_mfr(mfrdata) !computed path for each line of data read
        real*4 datetime
c
        real*8 day_fract_day
        real*8 d_iday
        real*8 d_fract_day
c
cd       print *,' '
cd       print *,'---in read_mfr ---'
c
c Identify the special period from 5 Feb. through 14 June 1998 
c for Bondville when only morning MFRSR data were good
c
        ch_start_date_bon = "980360000" ! Feb 5, 1998
        ch_end_date_bon = "981652358" ! June 14, 1998
c
c   convert these dates to 14time, base 1970
        call fnam_i4tim70(ch_start_date_bon,i4time70_start_bon,istatus)
        call fnam_i4tim70(ch_end_date_bon,i4time70_end_bon,istatus)
c
c   end special period definition for Bondville
c
c
        istatus_read_mfr = 1
c
        ihr_to_gmt = ista_to_gmt(isite)
cd       print *,'the value of ista_to_gmt has been set to ',ihr_to_gmt
cd       print *,'for station ',extensions(isite)
cd       print *,' '
c
        print *,' '
        print *,'Reading the mfrsr data'
c
        n_mfr_recs = 0
c
        do i = 1,mfrdata
c
           read(lun_mfr,*,end=120,err=132)day_fract_day,cosz,tot1,
     1      tot2,tot3,tot4,tot5,tot6,tot7,dif1,dif2,dif3,dif4,dif5,
     2      dif6,dif7,dir1,dir2,dir3,dir4,dir5,dir6,dir7
c
c  defined above:
c       real*8 day_fract_day
c       real*8 d_iday
c       real*8 d_fract_day
c  have to read the day_fract_day as real*8
c  below, that real*8 variable is broken apart into an int*4 day (iday)
c  and a real*4 fractional day (fract_day)
c
           iday = jidint(day_fract_day) !arg real*8, result(iday) int*4
           d_iday = dflotj(iday)        !arg int*4, result(d_iday) real*8
           d_fract_day = day_fract_day-d_iday !real*8 subtraction
           fract_day = sngl(d_fract_day) !arg real*8, result(fract_day) real*4
c
cd          print *,iday,fract_day,cosz,tot1,
cd    1      tot2,tot3,tot4,tot5,tot6,tot7,dif1,dif2,dif3,dif4,dif5,
cd    2      dif6,dif7,dir1,dir2,dir3,dir4,dir5,dir6,dir7
c
c    convert the input gmt date/time to i4time
c
cd          print *,' '
cd          print *,'input to date_to_i4time:'
cd          print *,'iday = ',iday
cd          print *,'fract_day = ',fract_day
c
           call date_to_i4time(iday,fract_day,i4time70_gmt)
cd          print *,' mfr i4time70_gmt = (from date_to_i4time)',
cd    1      i4time70_gmt
cd          call i4tim70_atime(i4time70_gmt,atime_mfr,istatus)
cd          if(atime_mfr(1:1).eq.' ')atime_mfr(1:1)='0'
cd          print *,'corresponding atime = ',atime_mfr
c
c
c Do not trust the SZA in the MFRSR file
c Compute an accurate solar zenith angle (subroutine needs UTC date/time)
c
           call i4tim70_int(i4time70_gmt,iyy_gmt,mon_gmt,iday_gmt,
     1 ihh_gmt,mm_gmt,isec_gmt,istatus)
c     
c  get the 4-digit year
           call i4tim70_atime(i4time70_gmt,atime_gmt,istatus)
           read(atime_gmt(8:11),15) yyyy
  15       format(f5.0)
c  compute the decimal hour
           dec_hr = float(ihh_gmt) + (float(mm_gmt)/60.)
c  get the day-of-year
           call i4tim70_fnam(i4time70_gmt,ch_yydddhhmm,istatus)
           read(ch_yydddhhmm(3:5),20)day_of_year
  20       format(f4.0)         
c
c   compute the solar zenith angle(sza) and path length (am in argument list)
c
           call sunpos(yyyy,day_of_year,dec_hr,xlat,xlon,az,el,sza,am,
     1 soldst)
c
d           if(i.le.10) then
d              print *,yyyy,day_of_year,dec_hr,xlat,xlon,az,el,sza,am
d           endif
c
c      convert the gmt i4time to local i4time 
cd          print *,'ihr_to_gmt = ',ihr_to_gmt
           i4time70_local_new = i4time70_gmt-(ihr_to_gmt*3600)
cd          print *,'adjusted i4time to local = ',
cd    1            i4time70_local_new
c
           call i4tim70_int(i4time70_local_new,local_year,
     1       local_month,local_day,local_hour,local_min,local_sec,
     2       istatus)
cd          print *,'adjusted to local time:'
cd          print *,local_year,local_month,local_day,local_hour,local_min,
cd    1           local_sec
c
           if(local_year.gt.93)then
              localyyyy = 1900 + local_year
           else
              localyyyy = 2000 + local_year
           endif
c
c   Test for the special period at Bondville - 5 Feb - 14 June 1998
c
           if(isite.eq.1) then
              if((i4time70_local_new.ge.i4time70_start_bon).and.
     1           (i4time70_local_new.le.i4time70_end_bon)) then
                 if(local_month.eq.5) then
                    go to 132 !skip all of May 1998 at Bondville
                 else if(local_hour.ge.12) then
                    go to 132 ! skip all afternoon data
                 endif
              endif
           endif
c
c   end special check for Bondville 1998 data
c 
           n_mfr_recs = n_mfr_recs + 1
c
           i4time70_local(n_mfr_recs) = i4time70_local_new
c
           ldate_mfr(n_mfr_recs) = 
     1          (localyyyy*10000)+(local_month*100)+local_day
           ltime_mfr(n_mfr_recs) = (local_hour*100)+local_min
           cosz_mfr(n_mfr_recs) = cosd(sza)
           path_len_mfr(n_mfr_recs) = am
           direct(n_mfr_recs,1) = dir2 ! 415 nm channel
           direct(n_mfr_recs,2) = dir3 ! 500 nm channel
           direct(n_mfr_recs,3) = dir4 ! 615 nm channel
           direct(n_mfr_recs,4) = dir5 ! 670 nm channel
           direct(n_mfr_recs,5) = dir6 ! 868 nm channel
c
cd          print *,' '
cd          print *,'ldate_mfr(n_mfr_recs) = ',ldate_mfr(n_mfr_recs) 
cd          print *,'ltime_mfr(n_mfr_recs) = ',ltime_mfr(n_mfr_recs) 
c
c   Escape point for bad data line
c
 132        continue
c
        enddo
c
 120    print *, 'MFRSR file read, ',n_mfr_recs,' good records found'
c
        print *,' '
        istatus_read_mfr = 0
c
d       print *,'ldate_mfr  ltime  computed cosZ   pathlen   direct3'
d       do i = 1,n_mfr_recs
d          print 125,ldate_mfr(i),ltime_mfr(i),cosz_mfr(i),
d    1     path_len_mfr(i),(direct(i,j),j=1,nchannel)
 125       format(1x,i8.8,1x,i4.4,1x,f7.5,1x,f6.3,<nchannel>(1x,f8.3))
d       enddo
c
        return !end reading MFRSR data
c
        end
c-----------------------------------------------------------------------
      subroutine sunpos(year,doy,dechr,lat,long,az,el,sza,am,
     1 soldst)
c                                                                     
c   Input to this algorithm is UTC
c
c   This subroutine calculates the local azimuth and elevation of the 
c   sun at a specific location and UTC time using an approximation to   
c   equations used to generate tables in The Astronomical Almanac.  
c   refraction correction is added so sun position is apparent one. 
c                                                                    
c   The Astronomical Almanac, U.S. Gov't Printing Office, Washington, 
c     D.C. (1985).                                                    
c                                                                    
c   input parameters                                                  
c     year=year, e.g., 1986                                           
c     doy=day of year, e.g., feb 1=32                                 
c     dechr=hours plus fraction in UT, e.g., 8:30 am eastern daylight  
c       time is equal to 8.5 + 5(5 hours west of Greenwich) -1(for    
c       daylight savings time correction)                             
c     lat=latitude in decimal degrees (north is positive)                    
c     long=longitude in decimal degrees (east is positive)                    
c                                                                     
c   output parameters                                                 
c     az=sun azimuth angle (measured east from north, 0 to 360 deg)   
c     el=sun elevation angle (deg)                                   
c     sza=solar zenith angle (deg)
c     plus others, but note the units indicated before return         
c                                                                    
c  work with real variables and define some constants, including
c  one to change between degs and radians
c
c
      implicit real (a-z)
c      
      data twopi,pi,rad/6.2831853,3.1415927,.017453293/
c
c compute the current "real" julian date (actually add 2,400,000 for jd)
c
      delta=year-1949.
      leap=aint(delta/4.)
      jd=32916.5+delta*365.+leap+doy+dechr/24.
c
c   1st no. is mid. 0 jan 1949 minus 2.4e6; leap=leap days since 1949
c  the last yr of century is not leap yr unless divisible by 400
      if(amod(year,100.).eq.0.0.and.amod(year,400.).ne.0.0)jd=jd-1.
c
c   calculate ecliptic coordinates
      time=jd-51545.0
c   51545.0 + 2.4e6 = noon 1 jan 2000
c
c   force mean longitude between 0 and 360 degs
c
      mnlong=280.460+.9856474*time
      mnlong=mod(mnlong,360.)
      if(mnlong.lt.0.)mnlong=mnlong+360.
c
c   mean anomaly in radians between 0 and 2*pi
c
      mnanom=357.528+.9856003*time
      mnanom=mod(mnanom,360.)
      if(mnanom.lt.0.)mnanom=mnanom+360.
      mnanom=mnanom*rad
c
c   compute the ecliptic longitude and obliquity of ecliptic in radians
c
      eclong=mnlong+1.915*sin(mnanom)+.020*sin(2.*mnanom)
      eclong=mod(eclong,360.)
      if(eclong.lt.0.)eclong=eclong+360.
      oblqec=23.439-.0000004*time
      eclong=eclong*rad
      oblqec=oblqec*rad
c
c   calculate right ascension and declination
c
      num=cos(oblqec)*sin(eclong)
      den=cos(eclong)
      ra=atan(num/den)
c
c   force ra between 0 and 2*pi
c
      if(den.lt.0)then
          ra=ra+pi
      elseif(num.lt.0)then
          ra=ra+twopi
      endif
c
c   dec in radians
      dec=asin(sin(oblqec)*sin(eclong))
c
c   calculate Greenwich mean sidereal time in hours
c
      gmst=6.697375+.0657098242*time+dechr 
c
c   hour not changed to sidereal time since 'time' includes
c   the fractional day 
c
      gmst=mod(gmst,24.)
      if(gmst.lt.0.)gmst=gmst+24.
c
c   calculate local mean sidereal time in radians 
c
      lmst=gmst+long/15.
      lmst=mod(lmst,24.)
      if(lmst.lt.0.)lmst=lmst+24.
      lmst=lmst*15.*rad
c
c   calculate hour angle in radians between -pi and pi
c
      ha=lmst-ra
      if(ha.lt.-pi)ha=ha+twopi
      if(ha.gt.pi)ha=ha-twopi
c
c   change latitude to radians
c
      lat=lat*rad
c
c   calculate azimuth and elevation
c
      el=asin(sin(dec)*sin(lat)+cos(dec)*cos(lat)*cos(ha))
      az=asin(-cos(dec)*sin(ha)/cos(el))
c
c   this puts azimuth between 0 and 2*pi radians
c
      if(sin(dec)-sin(el)*sin(lat).ge.0.)then
         if(sin(az).lt.0.)az=az+twopi
      else
         az=pi-az
      endif
c
cc   if az=90 degs, elcritical=asin(sin(dec)/sin(lat))
cc    elc=asin(sin(dec)/sin(lat))
cc    if(el.ge.elc)az=pi-az
cc    if(el.le.elc.and.ha.gt.0.)az=twopi+az
c
c   Calculate refraction correction for US stan. atmosphere
c   need to have el in degs before calculating correction
c
      el=el/rad
c
      if(el.ge.19.225) then 
         refrac=.00452*3.51823/tan(el*rad)
      else if (el.gt.-.766.and.el.lt.19.225) then
         refrac=3.51823*(.1594+.0196*el+.00002*el**2)/
     1   (1.+.505*el+.0845*el**2)
      else if (el.le.-.766) then
         refrac=0.0
      endif
c
c   note that 3.51823=1013.25 mb/288 C
c
c   elevation in degs
c
      el=el+refrac

c compute the solar zenith angle
c
      sza = 90.0 - el
c
c compute the air mass
c
      am=armass(el)
c
c   calculate distance to sun in A.U. & diameter in degs
c
      soldst=1.00014-.01671*cos(mnanom)-.00014*cos(2.*mnanom)
      soldia=.5332/soldst
c
c   convert az and lat to degs before returning
c
      az=az/rad
      lat=lat/rad
      ha=ha/rad
      dec=dec/rad
c
c   mnlong in degs, gmst in hours, jd in days if 2.4e6 added;
c   mnanom,eclong,oblqec,ra,and lmst in radians
c
      return
      end
c-----------------------------------------------------------------------
c
c   This function calculates air mass using kasten's
c   approximation to bemporad's tables
c
c
	function armass(el)
	z=(90.-el)*3.141592654/180.
	armass=1./(cos(z)+.50572*(6.07995+el)**-1.6364)
	return
	end
c-----------------------------------------------------------------------
      subroutine date_to_i4time(iday,fract_day,i4time70_gmt)
c
c This routine converts the inane MRFSR time to i4time (1970 base)
c note: MFRSR time is given as the number of days since Jan. 1 1900
c       plus the fractional day.
c
c  subtract the number of days between 1900 and Jan. 1, 1970 from
c  iday  
c
      integer*4 isec_in_years,isec_hours,isec_min
      integer*4 i4time_gmt
c
c     idays_since_1960 = iday - 21915
      idays_since_1970 = iday - 25568
c
c  extract the time from the fractional gmt day
cd     print *,'fract_day = ',fract_day
      ihour_gmt = int(fract_day*24.)
cd     print *,'int(ihour_gmt) = ',ihour_gmt
cd     print *,'float(hour_gmt) = ',float(ihour_gmt)
      frac_hour = (fract_day*24.) - float(ihour_gmt)
cd     print *,'frac_hour = ',frac_hour
      min = nint(frac_hour*60.)
cd     print *,'min = ',min
c
c  compute the number of seconds since Jan. 1, 1960 to the input date
c     isec_in_days = idays_since_1960*24*3600
      isec_in_days = idays_since_1970*24*3600
cd     print *,'isec_in_days = ',isec_in_days
      isec_hours = ihour_gmt*3600
      isec_min = min*60
c
      i4time70_gmt = isec_in_days + isec_hours + isec_min
c
      return
      end
c---------------------------------------------------------------------------
      subroutine initial(istatus_init)
c
      parameter (nsta=10)
      parameter (nchannel=5)
      parameter (nvals=5000)
c
c  Determine the site, dates of the period to be processed
c  and read the appropriate ozone file
c
      common/calibration/i_zero(nchannel),error(nchannel),
     1  constant_coef(nchannel),slope_coef(nchannel),
     2  sine_coef(nchannel),cos_coef(nchannel),
     3  err_slope_i_zero(nchannel),err_intercept_i_zero(nchannel)
c
      common/channels/wave_l(6),rayleigh(nchannel),mfrsr_sn,od_co2,
     1od_ch4,od_h2o_5cm,pwv_snd(4),i4time_pwv(4)
      common/unit/lun_mfr,lun_mfrhead,lun_aod,lun_o3,lun_p,lun_cof,
     1lun_hom,lun_err,lun_snd
      common/files/isite,extensions(nsta),names(nsta),lstring(nsta),
     1i4time_date,ista_to_gmt(nsta),pressure_nominal(nsta),xlat,xlon
      common/o3data/ozone(nvals),real_o3_i4time(nvals),n_rec
      common/ch1625/iserial(20),number_batch(20),batch(20),channel_1625
c
      logical*1 channel_1625
c
      real*4 i_zero
      real*4 intercept_i_zero
c
      character*1 ans
      character*24 names
      character*3 extensions
      character*24 atime_date
      character*1 batch
c
      integer*4 i4time_date
      integer*4 number_batch
c
      data (names(i),i=1,nsta)/'Bondville',
     1                         'Fort Peck',
     2                         'Goodwin Creek',
     3                         'Table Mountain',
     4                         'Desert Rock',
     5                         'Penn State',
     6                         'Sioux Falls',
     7                         'Canaan Valley',
     8                         'BND SURFRAD',
     9                         'FPE SURFRAD'/
c
      data (lstring(i),i=1,nsta)/9,9,13,14,11,10,11,13,11,11/
c
      data (extensions(i),i=1,nsta)/'bon',
     1                              'fpk',
     2                              'gwn',
     3                              'tbl',
     4                              'dra',
     5                              'psu',
     6                              'sxf',
     7                              'can',
     8                              'bnd',
     9                              'fpe'/
c
      data (pressure_nominal(i),i=1,nsta)/992.,
     1                                    942.,
     2                                    1007.,
     3                                    830.,
     4                                    900.,
     5                                    973.,
     6                                    960.,
     7                                    900.,
     8                                    992.,
     9                                    942./
c
      data (ista_to_gmt(i),i=1,nsta)/6,7,6,7,8,5,6,5,6,7/
c
c optical depths for 1625 channel absorption as function of MFRSR serial #,
c filter batch (a, b, c, or d), and sequence number within the batch 
c
      data (iserial(i),i=1,20)/632,633,648,649,650,651,659,660,661,
     c662,663,664,665,669,670,671,672,673,2*999/
c
c 1625 filter batch (a,b,c,or d)
      data (batch(i),i=1,20)/"a","a","b","b","b","b","b","b","c","c",
     1"c","c","c","b","b","b","d","d",2*"x"/
c
c 1625 number in batch (e.g., a1, a4, b1, b5,...)
      data (number_batch(i),i=1,20)/1,4,1,5,7,9,8,10,2,3,4,5,6,2,3,4,
     c10,11,2*999/
c
c old 1625 optical depths before air mass depencence ---------------
c     data (co2_392ppm(i),i=1,20)/0.0065,0.0072,0.0067,0.0064,0.0066,
c    c0.0063,0.0065,0.0067,0.0053,0.0056,0.0057,0.0056,0.0055,0.0067,
c    c0.0066,0.0069,0.0073,0.0072,2*0.0/
c
c     data (ch4_1_7ppm(i),i=1,20)/0.0036,0.0034,0.0036,0.0036,0.0036,
c    c0.0037,0.0036,0.0035,0.0047,0.0045,0.0044,0.0045,0.0045,0.0035,
c    c0.0036,0.0034,0.0033,0.0033,2*0.0/
c
c     data (h2o_5cm(i),i=1,20)/0.0055,0.0056,0.0055,0.0055,0.0055,
c    c0.0055,0.0055,0.0055,0.0058,0.0058,0.0058,0.0058,0.0058,0.0055,
c    c0.0055,0.0055,0.0056,0.0056,2*0.0/
c-------------------------------------------------------------------
c
      istatus_init = 1 !initialize subroutine status as bad
c
c initialize the logical variable channel_1625 to false
      channel_1625 = .FALSE. 
c
c  assign logical unit numbers
c
      lun_mfr = 10 !unit for MFRSR 500 nm data
      lun_mfrhead = 11
      lun_aod = 13
      lun_o3 = 14
      lun_p = 15
      lun_cof = 17
      lun_hom = 20
      lun_err = 21
      lun_snd = 22
c
c Choose a site
c
  5     print *,' '
        do i = 1,nsta
           write(6,8)names(i),i
  8        format(1x,a,t30,'[',i2,']')
        enddo
        print 9,nsta
  9     format(1x,'enter the number of the desired site (1-',i2,')')
        accept 10, isite
  10    format(i2)
d       print *,'site = ',isite
        if((isite.lt.0).or.(isite.gt.nsta)) then
           print *,'choices are only 1 through ',nsta,', try again'
           go to 5
        endif
c
c
c Choose a date to process
c
  20     print *,' '
         print *,' Enter a date to process, dd-mmm-yyyy: '
         print *,' '
         print *,'              e.g., 06-nov-1999'
         read(5,22)ina,atime_date
  22     format(q,a)
d        print *,'length of input string = ',ina
         if(ina .eq. 11) then
           continue
         else if((ina.eq.10).and.(atime_date(2:2).eq.'-').and.
     1    (atime_date(6:6).eq.'-')) then
           ina=11
           atime_date(1:11) = '0'//atime_date(1:10)
         else
           print *,'*ERROR* Need all 11 characters for dd-mmm-yyyy' 
           print *,'stop'
           stop
         endif
c
c Convert atime to yyyy, mm, dd, etc
         call atime_i4tim70(atime_date,i4time_date,istatus)
         call i4tim70_int(i4time_date,nyear,nmonth,
     1                 nday,nhour,nmin,nsec,istatus)
c
d        print *,'in initial'
d        print *,'atime_date = ',atime_date
d        print *,'nyear nmonth nday nhour nmin nsec'
d        print *,nyear,nmonth,nday,nhour,nmin,nsec
c
c
      print 30,atime_date
  30  format(/,1x,'the chosen time is ',a24,
     1//,1x,'  continue ? [<cr> or n]')
      accept 35, ans
  35  format(a1)
      if(ans.eq.'n')then
         go to 20
      else
         istatus_init = 0
d     print *,'leaving initial'
d     print *,' '
         return
      endif
c
      end
c-----------------------------------------------------------------------
      subroutine read_pres(fnam,sta,ihr_to_gmt,istatus_pres)
c
c This routine reads SURFRAD data from the file corresponding to the 
c date fnam and for the following date.  This is done because SURFRAD 
c data are in GMT days, and AOD processing is done in local standard 
c time.  The resultant two-day arrays of i4time (local, base 1970) and 
c station pressure are stored in common block baro.  Because the 
c i4time is used for interpolation, it is real in the common block.
c
      parameter (nsta=10,nbaro=3000,xmiss=-9999.9)
c
      common/unit/lun_mfr,lun_mfrhead,lun_aod,lun_o3,lun_p,lun_cof,
     1lun_hom,lun_err,lun_snd
      common/directory/in_directory,in_directory_mfr,out_directory,
     1 in_directory_pres,in_directory_o3,home_directory,snd_directory,
     2 lstring_in,lstring_mfr,lstring_out,lstring_p,lstring_o3,len_hom,
     3 len_snd
      common/baro/real_i4tim70_loc(nbaro),p(nbaro),n_pres
c
      integer i4time_p_local
c
      real*4 p
      real*4 real_i4tim70_loc
c
      character*9 fnam,fnam_day2
      character*33 filename
      character*80 station_name
      character*4 ch_yyyy
      character*2 ch_yy
      character*3 ch_jul
      character*3 sta
      character*4 ch_hhmm
      character*9 fnam_p
      character*80 in_directory,in_directory_mfr,out_directory
      character*80 in_directory_pres,in_directory_o3
      character*80 home_directory
      character*80 snd_directory
c
d     print *,' '
d     print *,'---in read_pres---'
      istatus_pres = 1
c
d     print *,'fnam day 1 = ',fnam
      call fnam_i4tim70(fnam,i4time_day1,istatus)
      i4time_day2 = i4time_day1 + (3600*24)
      call i4tim70_fnam(i4time_day2,fnam_day2,istatus)
d     print *,'fnam day 2 = ',fnam_day2
c
c read the SURFRAD data file for day 1
c    name the potential files
      filename = in_directory_pres(1:lstring_p)//sta//fnam(1:5)//
     1'.dat'
c
c  open the file
      open(unit=lun_p,file=filename,status='old',err=200) 
      print *,'file opened, reading SURFRAD data from ',filename
c
      n_pres = 0
c
c read header records
c
      read(lun_p,102)ina,station_name
 102  format(1x,q,a)
d     print *,' '
d     print *,'station name in header is ',station_name(1:ina)
c
      read(lun_p,*)xlat,xlon,elevation
      print *,xlat,xlon,elevation
d     print *,'In read_pressure '
d     print *,'latitude = ',xlat
d     print *,'longitude = ',xlon
d     print *,'elevation = ',elevation
d     print *,' '
c
c
      do i =1,nbaro
c
       read(lun_p,*,end=120)iyear,jday,month,iday,ihour,min,
     1 dt,zen,dpsp,iqc_dpsp,upsp,iqc_upsp,xnip,iqc_xnip,
     1 diffuse,iqc_diffuse,dpir,iqc_dpir,dpirc,iqc_dpirc,
     2 dpird,iqc_dpird,upir,iqc_upir,upirc,iqc_upirc,
     2 upird,iqc_upird,uvb,iqc_uvb,par,iqc_par,
     3 xnetpsp,iqc_xnetpsp,xnetpir,iqc_xnetpir,totalnet,
     3 iqc_totalnet,tc,iqc_tc,rh,iqc_rh,speed,iqc_speed,
     4 dir,iqc_dir,pres,iqc_pres
c
c change all values with bad qc flags to missing
c
         if(iqc_pres.eq.0) then
c
            n_pres = n_pres + 1
c
            p(n_pres) = pres
cd           print *,'pres = ',p(n_pres)
c
c construct an i4time for each line
            write(ch_yyyy,50)iyear
  50        format(i4.4)
            write(ch_jul,51)jday
  51        format(i3.3)
            write(ch_hhmm,52)ihour,min
  52        format(i2.2,i2.2)
            fnam_p = ch_yyyy(3:4)//ch_jul//ch_hhmm
            call fnam_i4tim70(fnam_p,i4time,istatus)
c        adjust the i4time to local
            i4time_p_local = i4time - (ihr_to_gmt*3600)
            real_i4tim70_loc(n_pres) = float(i4time_p_local)
cd           print *,'pres data i4time_local= ',real_i4tim70_loc(n_pres)
         endif
c
      enddo
      go to 120
c
 200  print *,filename,' not found, not used'
c
c
 120  close(unit=lun_p)
c
c read the data for the following day
c
      filename = in_directory_pres(1:lstring_p)//sta//
     1fnam_day2(1:5)//'.dat'
c
c  open the file
      open(unit=lun_p,file=filename,status='old',err=400) 
      print *,'file opened, reading SURFRAD data from ',filename
c
c read header records
c
 201  continue
      read(lun_p,102)ina,station_name
c     print *,' '
c     print *,'station name in header is ',station_name(1:ina)
c
      read(lun_p,*)xlat,xlon,elevation
d     print *,' '
d     print *,'latitude = ',xlat
d     print *,'longitude = ',xlon
d     print *,'elevation = ',elevation
d     print *,' '
c
c
      do i =1,nbaro
c
       read(lun_p,*,end=220)iyear,jday,month,iday,ihour,min,
     1 dt,zen,dpsp,iqc_dpsp,upsp,iqc_upsp,xnip,iqc_xnip,
     1 diffuse,iqc_diffuse,dpir,iqc_dpir,dpirc,iqc_dpirc,
     2 dpird,iqc_dpird,upir,iqc_upir,upirc,iqc_upirc,
     2 upird,iqc_upird,uvb,iqc_uvb,par,iqc_par,
     3 xnetpsp,iqc_xnetpsp,xnetpir,iqc_xnetpir,totalnet,
     3 iqc_totalnet,tc,iqc_tc,rh,iqc_rh,speed,iqc_speed,
     4 dir,iqc_dir,pres,iqc_pres
c
c change all values with bad qc flags to missing
c
         if(iqc_pres.eq.0) then
c
            n_pres = n_pres + 1
c
            p(n_pres) = pres
d           print *,'pres = ',p(n_pres)
c
c construct an i4time for each line
            write(ch_yyyy,50)iyear
            write(ch_jul,51)jday
            write(ch_hhmm,52)ihour,min
            fnam_p = ch_yyyy(3:4)//ch_jul//ch_hhmm
            call fnam_i4tim70(fnam_p,i4time,istatus)
c        adjust the i4time to local
            i4time_p_local = i4time - (ihr_to_gmt*3600)
            real_i4tim70_loc(n_pres) = float(i4time_p_local)
d           print *,'i4time_local = ',real_i4tim70_loc(n_pres)
         endif
c
      enddo
c
 400  print *,filename,' not found, not used'
c
 220  close (unit=lun_p)
c
      print *,n_pres,' records read'
c
      if (n_pres.le.1) then
         print *,' '
         print *,'no data in this file, return with bad read status'
         close (unit=lun_p)
         return
      endif
c
c  return with good status
      istatus_pres = 0
      return
c
      end
c------------------------------------------------------------------------------
      subroutine define_directories
c
c define the input and output directories
c
      parameter (nsta=10)
      parameter (nchannel=5)
c
      common/unit/lun_mfr,lun_mfrhead,lun_aod,lun_o3,lun_p,lun_cof,
     1lun_hom,lun_err,lun_snd
      common/directory/in_directory,in_directory_mfr,out_directory,
     1 in_directory_pres,in_directory_o3,home_directory,snd_directory,
     2 lstring_in,lstring_mfr,lstring_out,lstring_p,lstring_o3,len_hom,
     3 len_snd
c
      common/files/isite,extensions(nsta),names(nsta),lstring(nsta),
     1i4time_date,ista_to_gmt(nsta),pressure_nominal(nsta),xlat,xlon
c
      common/calibration/i_zero(nchannel),error(nchannel),
     1  constant_coef(nchannel),slope_coef(nchannel),
     2  sine_coef(nchannel),cos_coef(nchannel),
     3  err_slope_i_zero(nchannel),err_intercept_i_zero(nchannel)
c
      real*4 i_zero
      real*4 intercept_i_zero
c
      character*1 ans
      character*3 ch_ext
      character*24 names
      character*3 extensions
      character*4 ch_year
      character*80 in_directory,in_directory_mfr,out_directory
      character*80 in_directory_pres,in_directory_o3
      character*80 home_directory
      character*24 mfr_directory
      character*80 snd_directory
c
      integer*4 i4time_date
c
      home_directory = "/home/grad/surfrad/"
      len_hom = 19
c
c in_directory_pres now defined below, after ch_year is defined
c      in_directory_pres= '/q/surfrad/data/'
c      lstring_p = 16
c
      in_directory= '/home/grad/surfrad/aod/'
      lstring_in = 23
c
      in_directory_o3= '/home/grad/surfrad/aod/ozone/'
      lstring_o3 = 29
c
      snd_directory= '/q/surfrad/sounding/'
      len_snd = 20
c
c
c get the year for the root of the output file
      call i4tim70_int(i4time_date,nyear,nmm,ndd,nhh,nmin,nss,istatus)
c
c
c convert the 2-digit year to 4 digits
      if(nyear.gt.93)then
         nyyyy = 1900 + nyear
      else
         nyyyy = 2000 + nyear
      endif
c convert the 4-digit year to character
      write(ch_year,9)nyyyy
  9   format(i4.4)
c
      in_directory_pres= '/q/surfrad/data/'//ch_year//'/'
      lstring_p = 21
c
      print *,' '
      print *,'choose the input directory for MFRSR cccdata files'
c old      print *,' /u6/mfrsr/[sta]/mfrsr/ccc/               [<cr>]'
      print *,' /data/Inst/MFR/SURFRAD/[sta]/mfrsr/ccc/     [<cr>]'
      print *,' other (you will be requested to type it)    [o]'
      print *,' '
c
      accept 10, ans
  10  format(a1)
      if(ans.eq.' ')then
         in_directory_mfr='/data/Inst/MFR/SURFRAD/'//extensions(isite)//
     1'/mfrsr/ccc/'
         lstring_mfr= 37
      else
         print *,'enter the path to the mfrsr data directory'
         print *,' '
         accept 11, lstring_mfr,in_directory_mfr
  11     format(q,a80)
      endif
c
      print *,' '
      print *,'the input MFRSR directory is: ',
     1                   in_directory_mfr(1:lstring_mfr)
c
      print *,' '
      print *,'choose the output directory'
      print *,' /q/surfrad/aod/[sta]/[yyyy]/             [<cr>]'
      print *,' other (you will be requested to type it)    [o]'
      print *,' '
c
      ch_ext = extensions(isite)
      if(extensions(isite).eq.'bnd') then
         ch_ext = 'bon'
      else if(extensions(isite).eq.'fpe') then
         ch_ext = 'fpk'
      endif
c
      accept 10, ans
      if(ans.eq.' ')then
         out_directory = '/q/surfrad/aod/'//ch_ext//
     1'/'//ch_year//'/'
         lstring_out = 24
      else
         print *,'enter the complete path to the output directory'
         print *,'e.g. ./ for the working directory '
         accept 11, lstring_out,out_directory
      endif
c
      print *,' '
      print *,'the output directory is: ',out_directory(1:lstring_out)
      print *,' '
c
      return
      end
c-----------------------------------------------------------------------
      subroutine mfrsr_head(i4time_current,mfr_serial,lstring_serial)
c
c This routine determines the mfrsr head serial number for the current site 
c and time, and determines the fit coefficients for the associated optical
c depths of CO2, CH4 and H2O_5cm for 1625 nm if the current MFRSR has a
c 1625 channel.
c
      parameter (nsta=10)
      parameter (nchannel=5)
c
      common/channels/wave_l(6),rayleigh(nchannel),mfrsr_sn,od_co2,
     1od_ch4,od_h2o_5cm,pwv_snd(4),i4time_pwv(4)
      common/calibration/i_zero(nchannel),error(nchannel),
     1  constant_coef(nchannel),slope_coef(nchannel),
     2  sine_coef(nchannel),cos_coef(nchannel),
     3  err_slope_i_zero(nchannel),err_intercept_i_zero(nchannel)
      common/unit/lun_mfr,lun_mfrhead,lun_aod,lun_o3,lun_p,lun_cof,
     1lun_hom,lun_err,lun_snd
      common/directory/in_directory,in_directory_mfr,out_directory,
     1 in_directory_pres,in_directory_o3,home_directory,snd_directory,
     2 lstring_in,lstring_mfr,lstring_out,lstring_p,lstring_o3,len_hom,
     3 len_snd
      common/files/isite,extensions(nsta),names(nsta),lstring(nsta),
     1i4time_date,ista_to_gmt(nsta),pressure_nominal(nsta),xlat,xlon
      common/ch1625/iserial(20),number_batch(20),batch(20),channel_1625
      common/coefficients/ch4int,ch4x,ch4x2,co2int,co2x,co2x2,
     1h2oint,h2ox,h2ox2
c
      logical*1 channel_1625
c
      character*1 batch
      character*3 ch_ext
      character*24 names
      character*3 extensions
      character*11 mfr_head_hist_file
      character*9 ch_date_fmt_changed
      character*8 mfr_serial
      character*80 in_directory,in_directory_mfr,out_directory
      character*80 in_directory_pres,in_directory_o3
      character*80 home_directory
      character*80 snd_directory
c
      real*4 i_zero
      real*4 intercept_i_zero
c
      integer*4 date_fmt_changed
      integer*4 i4time_date
      integer*4 number_batch
c
      ch_ext = extensions(isite)
      if(extensions(isite).eq.'bnd') then
         ch_ext = 'bon'
      else if(extensions(isite).eq.'fpe') then
         ch_ext = 'fpk'
      endif
c
      mfr_head_hist_file = ch_ext//'_mfrhead'
      print *,' mfr_head_hist_file = ', mfr_head_hist_file
      open(unit=lun_mfrhead,file=in_directory(1:lstring_in)//
     1mfr_head_hist_file,status='old',err=100)
c
      print *,'MFRSR head file ',mfr_head_hist_file,' opened'
c
c Read loop for the MFRSR calibration file:
c
  10  read(lun_mfrhead,15,end=90)date_fmt_changed 
  15  format(i9.9)
      read(lun_mfrhead,20,end=90)ina,mfr_serial
  20  format(q,a)
      read(lun_mfrhead,*,end=90)(wave_l(i),i=1,6)
c
c  Read in slopes, intercepts, error slopes, and error intercepts
c  for all five channels of the MFRSR in the following order:
c  415nm,500nm,614nm or 1625nm, 670nm, and 868nm.
c
      do k = 1,nchannel
         read(lun_mfrhead,*,err=80)constant_coef(k),
     1             slope_coef(k),sine_coef(k),cos_coef(k)  
         read(lun_mfrhead,*,err=80)err_intercept_i_zero(k),
     1                             err_slope_i_zero(k)
      enddo
c
c  read past the 936 constants
      read(lun_mfrhead,*,err=80) dum,dum,dum,dum
      read(lun_mfrhead,*,err=80) dum dum
c
c Change character date to integer
      write(ch_date_fmt_changed,25)date_fmt_changed
  25  format(i9.9)
      call fnam_i4tim70(ch_date_fmt_changed,i4tim_change,istatus)
c
      if(i4time_current.gt.i4tim_change)then
c
         print 26,date_fmt_changed
  26     format('date format changed used for this date: ',i9.9)
         print *,'wavelengths read are ',(wave_l(i),i=1,6)
         lstring_serial=ina
c
c convert the MFRSR serial number character string to integer to
c load into common channel
c
c   Convert the character MFRSR serial number to integer
         read(mfr_serial,27)mfrsr_sn
  27     format(i<lstring_serial>)
c
         print 30,mfr_serial(1:lstring_serial)
  30     format(1x,'The mfrsr head is ',a,/) 
         print *,'I_zero constants used ',
     1         (constant_coef(k),k=1,nchannel)
         print *,'I_zero slope  coefficients used ',(slope_coef(k),
     1         k=1,nchannel)
         print *,'I_zero sine coefficients used ',(sine_coef(k),
     1         k=1,nchannel)
         print *,'I_zero cosine coefficients used ',(cos_coef(k),
     1         k=1,nchannel)
         print *,'error slopes used ',(err_slope_i_zero(k),k=1,nchannel)
         print *,'error intercept read ',(err_intercept_i_zero(k),k=1,
     1                                nchannel)
         print *,' '
         close(unit=lun_mfrhead)
         go to 50 ! to define 1625 channel coeffs if needed         
c
      else
c
         go to 10 ! try another MFRSR
c
      endif
c
  50  continue
c
c-----old 1625 od code-------------------------------------------
cc check whether this MFRSR has a 1625 channel; if so, define the 
cc absorption coefficients for co2, ch4, and h2o (valid for 5 cm)
cc
c      do i = 1,20
c         if(mfrsr_sn.eq.iserial(i))then
c            od_co2 = co2_392ppm(i)
c            od_ch4 = ch4_1_7ppm(i)
cc
c   Compute precipitible water vapor for four soundings, one 12h before
c   12 UTC of the day being analyzed, and the next three sounding times
c
c           print *,' '
c           print *,'Calling get_pwv for ',extensions(isite)
c           call get_pwv(isite,extensions(isite),i4time_date)
c           od_h2o_5cm = h2o_5cm(i)
c           return
c         endif
c      enddo
c------------------------------------------------------------------
c
c New code that gets 1625 optical depths for h20, ch4 and co2
c as a function of air mass. Here we only get the coefficients of
c the quadratic fits of optical depth as a function of air mass
c for each trace gas species. 
c Cannot compute optical depths here because air mass is needed
c for each measurement time.
c
c try to get the sequence number of the MFRSR in the iserial array
c this will only be successful if there is a 1625 channel
c
c     print *,' ' 
c     print *,'before testing if MFRSR has a 1625 channel'
c     print *,'mfrsr_sn = ', mfrsr_sn
      do i = 1,20
c        print *,'iserial(i) = ',iserial(i)
         if(iserial(i).eq.mfrsr_sn) then
c           print *,'after match found, iserial(i)',iserial(i)
c           print *,'mfrsr_sn = ', mfrsr_sn
            channel_1625 = .TRUE.
            index=i
            print *,' '
            print *,'This MFRSR has a 1625 channel'
            print *,'MFRSR ',mfr_serial,' found in iserial array'
            print *,'iserial (',index,') found as ',mfrsr_sn
c   MFRSR has a 1625 channel, get the ch4, co2 and h2o fit coefficients
c   and load them into common block "coefficients"
            if(batch(index).eq."a") call coeff_a(number_batch(index))
            if(batch(index).eq."b") call coeff_b(number_batch(index))
            if(batch(index).eq."c") call coeff_c(number_batch(index))
            if(batch(index).eq."d") call coeff_d(number_batch(index))
c
c   Compute integrated water vapor for four soundings, one 12h before
c   12 UTC of the day being analyzed, and the next three sounding times
            print *,' '
            print *,'Calling get_pwv for ',extensions(isite)
            call get_pwv(isite,extensions(isite),i4time_date)
            return
         endif
      enddo
c
      print *,' '
      print *,'This MFRSR must not have a 1625 channel'
      print *,' cannot find ' ,mfr_serial,' in iserial array'
      print *,' '
c
      print *,' '
      print *,'MFRSR ',mfrsr_sn,' does not have a 1625 ch, set optical',
     1' depths for H2O, CO2, and CH4 to 0, and sounding PWV info to miss
     2ing'
      print *,' '
c
      channel_1625 = .FALSE.
      od_co2 = 0.0
      od_ch4 = 0.0
      od_h2o_5cm = 0.0
      do i = 1,4
         pwv_snd(i) = -999.00
         i4time_pwv(4) = 0
      enddo
c
      return
c
  80  print *,' '
      print *,'Error: mfrhead file has wrong structure, no calibration'
      print *,'coefficient information, stop'
      stop
c
  90  print *,' '
      print *,'Error:  end of mfrsr history file reached before'
      print *,'expected, could not assign a head number ... stop'
      stop
c
  100 print *,' '
      print *,'Error: could not open the file mfr_head_hist_file'
      stop
c
      end
c-----------------------------------------------------------------------
      function trpt(p2,a,n,tp,xmiss)
c
c     trpt- interpolation within an array according to a second parameter
c
c        returns -999. if interpolation not possible
c
c     p2- second parameter array (real)
c     a- array of parameter to be interpolated 
c     n- number of points in p2 and a
c     tp- date/time (real i4time)  to interpolate to
c     xmiss - missing value in input data
c
      dimension p2(1),a(1)
      if(tp.gt.p2(n))go to 5
c
c
      nm1=n-1
      do 10 i=1,nm1
      if(tp.lt.p2(i)) go to 5
      if(tp.gt.p2(i+1)) go to 10
c   found the time interval to interpolate across
c      first check if one of the two is the exact time needed
	if(tp.eq.p2(i))then
	   trpt = a(i)
	   return
	else if(tp.eq.p2(i+1)) then
	   trpt = a(i+1)
	   return
	else
	   continue
	endif
c    check if the two levels each have bad data
	if((a(i) .eq. -999.) .or. (a(i+1) .eq. -999.)) then
	   trpt = -999.
	   return
	endif
	if((a(i) .eq. xmiss) .or. (a(i+1) .eq. xmiss)) then
	   trpt = -999.
	   return
	endif
c
c   interploation possible between p2(i) and p2(i+1)
      x1=(p2(i))
      x2=(p2(i+1))
      x=(tp)
      trpt=(a(i)*(x-x2)-a(i+1)*(x-x1))/(x1-x2)
      return
   10 continue
    5 trpt=-999.
      return
      end
c-----------------------------------------------------------------------
      subroutine get_pwv(isite,ch_sta,i4time_date)
c
c This subprogram computes integrated water vapor at station site for 
c 00Z and 12Z of the current date (i4time_date), 00Z of the next day, and
c for the following 12Z (designated 24 in variable names). The data 
c source is the interpolated soundings at the SURFRAD sites. 
c
c Construct the interpolated soundings filenames
c
      parameter (nsta=10)
      parameter (xmiss = -999.00)
      parameter (nchannel=5)
c
      common/channels/wave_l(6),rayleigh(nchannel),mfrsr_sn,od_co2,
     1od_ch4,od_h2o_5cm,pwv_snd(4),i4time_pwv(4)
      common/unit/lun_mfr,lun_mfrhead,lun_aod,lun_o3,lun_p,lun_cof,
     1lun_hom,lun_err,lun_snd
      common/directory/in_directory,in_directory_mfr,out_directory,
     1 in_directory_pres,in_directory_o3,home_directory,snd_directory,
     2 lstring_in,lstring_mfr,lstring_out,lstring_p,lstring_o3,len_hom,
     3 len_snd
c
      character*8 ch_yyyymmdd
      character*4 ch_yyyy
      character*15 snd_filename
      character*40 filename_in
      character*80 in_directory,in_directory_mfr,out_directory
      character*80 in_directory_pres,in_directory_o3
      character*80 home_directory
      character*80 snd_directory
      character*15 surfname(nsta)
      character*15 ch_site
      character*15 ch_site_read_in
      character*24 atime_in
      character*3 ch_sta
      character*2 ch_hh
c
      real*4 pres(38)
      real*4 temp(38)
      real*4 dpt(38)
      real*4 m_r(38)
c
      data (surfname(i),i=1,nsta)/'Bondville','Fort Peck',
     1 'Goodwin Creek','Table Mountain','Desert Rock','Penn State',
     1 'Sioux Falls','Canaan Valley','Bondville','Fort Peck'/
c
c define the site chosen
      if(isite.eq.1)then
        ch_site = surfname(1)
      else if (isite.eq.2)then
        ch_site = surfname(2)
      else if (isite.eq.3)then
        ch_site = surfname(3)
      else if (isite.eq.4)then
        ch_site = surfname(4)
      else if (isite.eq.5)then
        ch_site = surfname(5)
      else if (isite.eq.6)then
        ch_site = surfname(6)
      else if (isite.eq.7)then
        ch_site = surfname(7)
      else if (isite.eq.8)then
        ch_site = surfname(8)
      else if (isite.eq.9)then
        ch_site = surfname(9)
      else if (isite.eq.10)then
        ch_site = surfname(10)
      else
         print *,' '
         print *,'ERROR in get_pwv, site num ',isite,' not recognized'
         stop
      endif
c
      print *,'in get_pwv, ch_site = ',ch_site
c
c Start processing 
c
c Compute PWV from 4 interpolated soundings, one 12h
c before 12 UTC of the current date, and the three following soundings
c
c loop through the four sounding times. The begin time is i4time_date at 0000 UTC
c
d     print *,'in get_pwv, i4time_date = ',i4time_date
c
c    define the end (4th) sounding time
c
      i_end_time = i4time_date + 3*(12*3600)
      icount = 0 ! initialize the sounding counter
c
c loop through the 4 soundings, computing PWV for each
c
      do k = i4time_date,i_end_time,(12*3600)
c
         call i4tim70_int(k,iyy,mm,idd,ihh,min,isec,istatus)
         icount = icount + 1
c
c load the time of this sounding into the i4time_pwv array
         i4time_pwv(icount) = k
c
c   convert date parameters to character
         if(iyy.le.94)then
            iyyyy = 2000+iyy
         else
            iyyyy = 1900+iyy
         endif
c
         iyyyymmdd = (iyyyy*10000)+(mm*100)+idd
c
c convert integer to character
         write(ch_yyyymmdd,10)iyyyymmdd
  10     format(i8.8)
         write(ch_yyyy,15)iyyyy
  15     format(i4.4)
         write(ch_hh,16)ihh
  16     format(i2.2)
c
d        print *,'in get_pwv, ch_yyyy = ',ch_yyyy
d        print *,'        ch_yyyymmdd = ',ch_yyyymmdd
c
c Read and process the current interpolated sounding appropriate to the site chosen
c
         snd_filename = ch_yyyymmdd//'_'//ch_hh//'.int'
d        print *,'in get_pwv, snd_filename= ',snd_filename
         filename_in = snd_directory(1:len_snd)//ch_yyyy//'/'//
     1snd_filename
d        print *,'in get_pwv, filename_in= ',filename_in
c 
         open(lun_snd,name=filename_in,status='old',err=100,readonly)
         print *,' '
         print *,'sounding file ',filename_in,' opened in get_pwv'
c
         read(lun_snd,150) num_snd,num_levels,atime_in,npass,xl
 150     format(1x,i4,1x,i3,1x,a24,1x,i2,1x,f8.2)
c
c read the output, one interpolated sounding at a time
c
         do i = 1,num_snd
c
c read the header for this interpolated location
c
            read(lun_snd,152) ch_site_read_in,g_lat,g_lon,ielev
 152        format(1x,a15,1x,f6.3,1x,f8.3,1x,i4)
c
            if(ch_site_read_in.eq.ch_site) then
            print *,ch_site_read_in,' equals ',ch_site
c
               do j = 1,num_levels
                  read(lun_snd,*)p_g,hgt_g,temp_g,dwpt_g,u_g,v_g
c155              format(6(1x,f8.2)) ! do not use, use free format read instead
                  pres(j) = p_g
                  if(temp_g.ne.xmiss) then
                     temp(j) = temp_g + 273.16
                  else
                     temp(j) = xmiss
                  endif
                  if(dwpt_g.ne.xmiss) then
                     dpt(j) = dwpt_g + 273.16
                  else
                     dpt(j) = xmiss
                  endif
               enddo
               close (unit=lun_snd)
               go to 30 ! escape read-sounding loop,compute PWV 
c
            else ! not the requested site,
c               skip this sounding
c
               do j = 1,num_levels
                  read(lun_snd,*)p_g,hgt_g,temp_g,dwpt_g,u_g,v_g
               enddo
c
            endif
c
         enddo !end station loop
c
 100     print *,' '
         print *,'ERROR, ',ch_hh,' interpolated sounding for ',ch_site
         print *,'on ',ch_yyyymmdd,' not found in ',filename_in
         print *,'set ',ch_hh,' UTC PWV to missing'
c
c Assign missing precipitible water vapor
c
         pwv_snd(icount) = xmiss
         close (unit=lun_snd)
         go to 200 ! move to the next sounding
c
  30     continue
c
d        print *,' '
d        print *,'Compute the',ch_hh,' UTC PWV for ',ch_site
c
c Compute the integrated water vapor for the current sounding 
c
c compute the mixing ratio at each level
c
         do i = 1,num_levels
            if((dpt(i).ne.xmiss).and.(pres(i).ne.xmiss)) then
               m_r(i) = w(dpt(i),pres(i))
            else
               m_r(i) = xmiss
            endif
d           print *, pres(i),temp(i),dpt(i),m_r(i)
         enddo
c
c integrate the water vapor vertically
c
c   find a good place to start
         do i = 1,num_levels
            if((m_r(i).ne.xmiss).and.(pres(i).ne.xmiss)) then
               ibegin = i
               go to 40
            endif
         enddo
         print *,' '
         print *,'level to start vertical integration not found, skip'
         pwv_snd(icount) = xmiss
         go to 200
c
c Integrate vertically
c
  40     continue
c
d        print *,' '
d        print *,'ch_hh integration will start at level ',pres(ibegin)
c
         w_befor = m_r(ibegin)
         p_befor = pres(ibegin)
         sum = 0.0
         istart = ibegin+1
d              print *,'w_bottom w_top w_mean p_bottom p_top p_diff sum(
d    1p_diff*w_mean)'
         do i = istart,num_levels
            if((m_r(i).ne.xmiss).and.(pres(i).ne.xmiss)) then
               w_mean = (w_befor + m_r(i))/2.
               p_diff = p_befor - pres(i)
               sum = sum + (p_diff*w_mean)
d              print 95,w_befor,m_r(i),w_mean,p_befor,pres(i),p_diff,sum
d 95           format(3(1x,f8.5),2(1x,f5.1),1x,f4.1,1x,f7.2)
               w_befor = m_r(i)
               p_befor = pres(i)
            endif
         enddo
c
         pwv_snd(icount) = sum/980.
c
         print *,'The ',ch_hh,'UTC integrated water vapor is ',
     1           pwv_snd(icount),' cm'
         print *,' '
c
c----end PWV calculation-----
c
 200     continue
c
      enddo ! end the loop through the 4 sounding times
c
      return
      end
c------------------------------------------------------------------------------
c Thermodynamic functions
c-------------------------------------------------------------------------
      function  esat(t)
c  esat(millibars),t(kelvin)
      data abz/273.16/
      tc=t-abz
      esat=6.1078*exp((17.2693882*tc)/(tc+237.3))
      return
      end
c---------------------------------------------------------------------------
      function w(t,p)
c  w(grams water vapor/kilogram dry air ),temp(K), p(millibar )
      if((t.ge.999.).or.(t.lt.-990.))go to 10
         x  =  esat(t)
         w  =  621.97*x/(p-x)
      return
   10 w=0.0
      return
      end
c---------------------------------------------------------------------------
