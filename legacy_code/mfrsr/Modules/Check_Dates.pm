
sub Check_Dates () {
    use strict;
    use Time::DaysInMonth;

    my $start = shift(@_);
    my $end   = shift(@_);

    my $days;

    if ($start !~ /^\d{8}$/) {
	die "\nStart date $start must be in YYYYMMDD format";
    }
    if ($end !~ /^\d{8}$/) {
	die "\nEnd date $end must be in YYYYMMDD format";
    }

    if ($end < $start) {die "\nStart date must be before or equal to end date";}

    $start =~ /(\d\d\d\d)(\d\d)(\d\d)/;
    my $beg_y = $1;
    my $beg_m = $2;
    my $beg_d = $3;
    
    $end   =~ /(\d\d\d\d)(\d\d)(\d\d)/;
    my $end_y = $1;
    my $end_m = $2;
    my $end_d = $3;

    if ($beg_y < 1990) {
	die "Begin year $beg_y out of range.  Check";
    }
    if ($end_y > 2100) { # Should change 2100 to current year
	die "End year $end_y out of range.  Check";
    }

    if ($beg_m < 1 || $beg_m > 12) {
	die "Begin month $beg_m out of range.  Check";
    }
    if ($end_m < 1 || $end_m > 12) {
	die "End month $end_m out of range.  Check";
    }
    
    $days = days_in($beg_y,$beg_m); # Calculating number days in begin YYYYMM
    if ($beg_d < 1 || $beg_d > $days) {
	die "\nBegin day $beg_d out of range.  There are $days days in $beg_y$beg_m";
    }

    $days = days_in($end_y,$end_m); # Calculating number days in end YYYYMM
    if ($end_d < 1 || $end_d > $days) {
	die "\nEnd day $end_d out of range.  There are $days days in $end_y$end_m";
    }
}
1;

