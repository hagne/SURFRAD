###############################################################################
# Subroutine sunae
#
# Calculates the Azimuth, Elevation, Solar Distance and Disk Diameter
#
# Fixed some errors in the refraction code on 20111108
#
###############################################################################
sub SUNAE () {
    use Math::Trig;
#    use strict;
    my $year = shift(@_);
    my $day  = shift(@_);
    my $hour = shift(@_); # Fractional hour
    my $lat  = shift(@_);
    my $long = shift(@_);

    my ($refrac);

    if ($year < 1950 || $year > 2050) {
        die "SUNAE--bad input variable YEAR ($year), stopped\n";
    }
    if ($day < 1 || $day > 366) {
        die "SUNAE--bad input variable DAY ($day), stopped\n";
    }
    if ($hour < -13.0 || $hour > 36.0) {
        die "SUNAE--bad input variable HOUR ($hour), stopped\n";
    }
    if ($lat < -90 || $lat > 90) {
        die "SUNAE--bad input variable LAT ($lat), stopped\n";
    }
    if ($long < -180 || $long > 180) {
        die "SUNAE--bad input variable LON ($long), stopped\n";
    }

    my $pi    = atan2(1,1) * 4;
    my $twopi = 2.0 * $pi;
    my $rpd   = $pi/180;

#                    ** current Julian date (actually add 2,400,000
#                    ** for true JD);  LEAP = leap days since 1949;
#                    ** 32916.5 is midnite 0 jan 1949 minus 2.4e6

    my $delta = $year - 1949;
    my $leap  = int($delta / 4);
    my $jd = 32916.5 + ($delta*365 + $leap + $day) + ($hour/24.);
#
#                    ** last yr of century not leap yr unless divisible
#                    ** by 400 (not executed for the allowed YEAR range,
#                    ** but left in so our successors can adapt this for
#                    ** the following 100 years)


    if (($year % 100 == 0) && ($year % 400 != 0)) {
        $jd = $jd - 1;
    }

#                     ** ecliptic coordinates
#                     ** 51545.0 + 2.4e6 = noon 1 jan 2000

    my $time = $jd - 51545.0;

#                    ** force mean longitude between 0 and 360 degs

    my $mnlong = 280.460 + (0.9856474 * $time);
#    $mnlong = ($mnlong % 360);
    while ($mnlong > 360) { $mnlong = $mnlong - 360; }

#                    ** mean anomaly in radians between 0 and 2*pi

    my $mnanom = 357.528 + (0.9856003 * $time);
#    $mnanom = ($mnanom % 360);
    if ($mnanom < 0) { $mnanom = $mnanom + 360; }
    $mnanom = $mnanom * $rpd;

#                    ** ecliptic longitude and obliquity
#                    ** of ecliptic in radians

    my $eclong = $mnlong + 1.915*sin($mnanom) + 0.020*sin(2*$mnanom);
    while ($eclong > 360) { $eclong = $eclong - 360; }

    my $oblqec = 23.439 - 0.0000004*$time;
    $eclong = $eclong*$rpd;
    $oblqec = $oblqec*$rpd;

#                    ** right ascension

    my $num  = cos($oblqec)*sin($eclong);
    my $den  = cos($eclong);
    my $temp = $num/$den;
    my $ra   = atan($num/$den);

#                    ** Force right ascension between 0 and 2*pi

    if ($den < 0) {
        $ra  = $ra + $pi;
    } elsif ($num < 0) {
        $ra  = $ra + $twopi;
    }

#                    ** declination

    my $dec = asin(sin($oblqec)*sin($eclong));

#                    ** Greenwich mean sidereal time in hours

    my $gmst = 6.697375 + 0.0657098242*$time + $hour;

#                    ** Hour not changed to sidereal time since
#                    ** 'time' includes the fractional day

    while ($gmst > 24) { $gmst = $gmst - 24; }

#                    ** local mean sidereal time in radians

    my $lmst  = $gmst + $long/15;
    while ($lmst > 24) { $lmst = $lmst - 24; }
    $lmst   = $lmst * 15 * $rpd;

#                    ** hour angle in radians between -pi and pi

    my $ha  = $lmst - $ra;
    if ($ha < -$pi) { $ha  = $ha + $twopi; }
    if ($ha > $pi) { $ha  = $ha - $twopi; }

#                    ** solar azimuth and elevation

    my $el = asin(sin($dec) * sin($lat * $rpd) +
               cos($dec) * cos($lat * $rpd) * cos($ha));

    my $az = asin(-cos($dec) * sin($ha) / cos($el));

#                    ** Put azimuth between 0 and 2*pi radians

    if ((sin($dec) - sin($el) * sin($lat * $rpd)) > 0) {
        if (sin($az) < 0) { $az = $az + $twopi; }
    } else {
        $az = $pi - $az;

    }

#                     ** Convert elevation and azimuth to degrees
    $az  =  $az / $rpd;
    $el  =  $el / $rpd;

#  ======== Refraction correction for U.S. Standard Atmos. ==========
#      (assumes elevation in degs) (3.51823=1013.25 mb/288 K)

    if ($el > 19.225 ) {
        $refrac = 0.00452*3.51823 / tan($el*$rpd);
    } elsif ($el > -0.766 && $el < 19.225 ) {
        $refrac = 3.51823 * (0.1594 + 0.0196*$el + 0.00002*$el**2) /
            (1 + 0.505*$el + 0.0845*$el**2);

    } elsif ($el <= -0.766 ) {
        $refrac = 0.0;
    }
    $el = $el + $refrac;

# ===================================================================

#                   ** distance to sun in A.U. & diameter in degs

    my $soldst = 1.00014 - 0.01671*cos($mnanom) - 0.00014*cos(2*$mnanom);
    my $soldia = 0.5332 / $soldst;

    if ($el < -90 || $el > 90) {
        die "SUNAE--Output argument EL ($el) out of range\n";
    }
    if ($az < 0 || $az > 360) {
        die "SUNAE--Output argument AZ out of range\n";
    }
    my @returnary = ($az, $el, $soldst, $soldia);
    return (@returnary);
}
1;
