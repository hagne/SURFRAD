########################################################################
# Determines if time range entered is valid or exceeds specified number
# of days.  Would like to change the logic to check for number of whole
# months instead of days.
#
sub Valid_Date_Range () {
    use strict;
    use Time::Local;
    my $beg = shift(@_); # Begin date YYYYMMDD
    my $end = shift(@_); # End date YYYYMMDD
    my $val = shift(@_); # Valid number of days

    my $yyyy = substr($beg,0,4);
    my $mm   = substr($beg,4,2);
    my $dd   = substr($beg,6,2);
    my $begs = timegm(0,0,0,$dd,$mm-1,$yyyy);

    $yyyy = substr($end,0,4);
    $mm   = substr($end,4,2);
    $dd   = substr($end,6,2);
    my $ends = timegm(0,0,0,$dd,$mm-1,$yyyy);

    my $days = ($ends - $begs) / 86400;
    if ($days > $val) {
	print "\n                      *****  ERROR  *****\n\n";
	print "  Entered time range is: ${beg}-${end} which is $days days\n\n";
	print "  Date range is too large.  Files sometimes are not created\n";
	print "  correctly while processing much more than $val days at a\n";
	print "  time.  Please re-run with less than six months or $val days";
	print "\n\n";
	exit 1;
    }
}
1;
