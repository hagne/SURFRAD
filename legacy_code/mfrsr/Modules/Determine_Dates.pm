sub Determine_Dates () {
    use strict;
    use Time::DaysInMonth;

    my $start = shift(@_);
    my $end   = shift(@_);

    $start =~ /(\d\d\d\d)(\d\d)(\d\d)/;
    my $beg_y = $1;
    my $beg_m = $2;
    my $beg_d = $3;
    
    $end   =~ /(\d\d\d\d)(\d\d)(\d\d)/;
    my $end_y = $1;
    my $end_m = $2;
    my $end_d = $3;

    my $year_calc  = (($end_y - $beg_y)*12);
    my $num_months = (($end_m + $year_calc) - $beg_m) + 1; # Num of months

    my $cnt  = 1;
    my $year = $beg_y;
    my $mon  = $beg_m;
    my @YYYYMM = ();
    while ($cnt < ($num_months + 1)) {
	if ($mon == 13) {
	    $year++;
	    $mon = 1;
	}
	push @YYYYMM, (($year*100)+$mon);
	$mon++;
	$cnt++;
    }

    my @dates;
    my $chk = "Y";
    foreach (@YYYYMM) {
	my ($cnt,$dim);
	
# Setting the first YYYYMMDD days to begin_days	
	if ($chk eq "Y") {
	    $cnt = $beg_d;
	    $chk = "N";
	} else {
	    $cnt = 1;
	}

	my $yyyy = substr($_,0,4);
	my $mm   = substr($_,4,2);

# Setting days_in_month (dim) based on YYYYMM or end_days
	if ( (($yyyy*100)+$mm) != (($end_y*100)+$end_m) ) {
	    $dim = days_in($yyyy,$mm);
	} else {
	    $dim = $end_d;
	}

	while ($cnt < $dim+1) {
	    push @dates, (($yyyy*10000)+($mm*100)+$cnt);
	    $cnt++;
	}
    }
    return @dates;
}
1;
