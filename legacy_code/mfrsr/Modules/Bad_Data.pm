sub Bad_Data () {
    use strict;

    my $site = shift(@_);
    my $beg  = shift(@_); # Begin time
    my $end  = shift(@_); # End time

###############################################################################
# Begin and end times are in YYYYMMDD format.  Need to convert them to
# YYYYMMDDHHMNSS with the end set to the following day. 
#
    my $beg = sprintf "%8.8d%6.6d",$beg,"0";
    my @a;
    $a[0] = substr($end,0,4);
    $a[1] = substr($end,4,2);
    $a[2] = substr($end,6,2);
    @a = (gmtime(timegm(0,0,0,$a[2],$a[1]-1,$a[0])+86400))[5,4,3];
    $a[0]+=1900;
    $a[1]++;
    $end = sprintf"%4.4d%2.2d%2.2d%6.6d",@a,"0";

    my $file = "/home/grad/mfrsr/Calinfo/bad_data";
    open(BD, $file) or die "Can't open $file: $!\n";

    my (@allperiods, @subperiods);
    while (<BD>) {
        chomp;
        next if /^\s*$/;     # Skip blank lines                             
        next if /^\s*\#/;    # Skip comment lines                           
        s/\#.*$//;           # Strip tag comment  
	s/^\s+//;            # Strip preceding white space
	if (/^$site/) {
	    push @allperiods, $_;
	}
    }

    foreach my $line (@allperiods) {
	my $tB = "F"; # Test Begin
	my $tE = "F"; # Test End
	my $tS = "F"; # Test Span.  Does it span entire range?
	my @flds  = split(/\s+/,$line);
	my $b = $flds[2];
	my $e = $flds[3];

	if (($beg >= $b) && ($beg <= $e)) {$tB = "T";}
	if (($end >= $b) && ($end <= $e)) {$tE = "T";}
	if (($beg <= $b) && ($end >= $e)) {$tS = "T";}

	if (($tB eq "T") || ($tE eq "T") || ($tS eq "T")) {
	    push @subperiods, $line;
	} 
    }

    my %bd; # Hash to hold keys of bad data times
    foreach my $line (@subperiods) {
	my @flds  = split(/\s+/,$line);
	my $code = $flds[1];
	my $b = $flds[2]; # Beg time of bad data
	my $e = $flds[3]; # End time of bad data

###############################################################################
# At this point the begin and end times need to be determined for the while
# loop below.  The overall begin and end times are compared with the beg/end
# times from the bad_data file.  Those are compared to determine the most
# efficient beg/end times.
#
	my $wbeg = $beg;  # Beg/end vars for the while loop below
	my $wend = $end;
	if ($beg <= $b) {$wbeg = $b;}
	if ($end >= $e) {$wend = $e;}

# Now to set the sub-begin ($sb) and sub-end ($se) variables
	my $sb = $flds[4]; # Sub-Begin.
	my $se = $flds[5]; # Sub-End. To subset hours within day.
	if ($sb > $se) {die "Sub_Begin must be <= Sub_End";}

###############################################################################
# The while loop that tests the four beg/end times and assigns times (keys)
# as appropriate to the Bad Data hash (%bd).  Note the key is $wbeg
#
	while ($wbeg <= $wend) {
	    my $tsb = substr($wbeg,8,6);
	    if (($sb == 999999) || ($se == 999999)) {
		$bd{$wbeg} = $code;
	    } elsif (($tsb >= $sb) && ($tsb <= $se)) {
		$bd{$wbeg} = $code;
	    }

# Increment $wbeg by 20 sec each iteration of the while loop.
	    my $yyyy = substr($wbeg,0,4);
	    my $mm   = substr($wbeg,4,2);
	    my $dd   = substr($wbeg,6,2);
	    my $hh   = substr($wbeg,8,2);
	    my $mn   = substr($wbeg,10,2);
	    my $ss   = substr($wbeg,12,2);
	    my @t = (gmtime(timegm($ss,$mn,$hh,$dd,$mm-1,$yyyy)+20));
	    $t[5]+=1900;
	    $t[4]++;
	    $wbeg = sprintf "%4.4d%2.2d%2.2d%2.2d%2.2d%2.2d",reverse@t[0..5];
	}
    }
    return (%bd);
}
1;
