
sub Check_Date () {
    use strict;
    use Time::DaysInMonth;

    my $date = shift(@_);

    if ($date !~ /^\d{8}$/) {
	die "\n ERR: Date $date must be in YYYYMMDD format:\n Check";
    }

    chomp (my $today = `date +%Y%m%d`);

    if ($date > $today) {die "\n ERR: Cannot plot the future!\n Input $date must be <= ${today}:\n Check";}

    $date =~ /(\d\d\d\d)(\d\d)(\d\d)/;
    my $date_y = $1;
    my $date_m = $2;
    my $date_d = $3;
    
    if ($date_y < 1990) {
	die "\n ERR: Date year $date_y out of range, invalid prior to 1990.  \n\n Check";
    }

    if ($date_m < 1 || $date_m > 12) {
	die "\n ERR: Date month $date_m invalid.\n Check";
    }
    
    my $days = &days_in($date_y,$date_m); # Calc number days in date YYYYMM
    if ($date_d < 1 || $date_d > $days) {
	die "\n ERR: Day $date_d invalid. There are $days days in $date_y$date_m  \n\n Check";
    }
}
1;

