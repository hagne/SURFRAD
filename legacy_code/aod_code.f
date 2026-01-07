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
