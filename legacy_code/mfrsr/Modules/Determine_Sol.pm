# Subroutine Determine_Sol determines the correct cosine file to use
# for the site and date(s).  It does this by interrogating the .sol files
# located in /home/grad/mfrsr/Calinfo/SOL .  The information needed to 
# determine the correct file to use is imbedded in the filename.
#
sub Determine_Sol () {
    use strict;    

    my $site = shift(@_);
    my $type = shift(@_);
    my $date = shift(@_);

    my $path = "/home/grad/mfrsr/Calinfo/SOL";
    my @listing = `ls -1 ${path}/${site}_${type}_????????_????_????_???.sol`;
    chomp @listing;
    @listing = reverse@listing;

    my $solfile;
  OUTER: 
    foreach my $file (@listing) {
	my $soldate = (split '_',$file)[2]; # Extracting date of .sol file
	if ($soldate <= $date) {
	    $solfile = $file;
	    last OUTER;
	}
    }

    if (-f $solfile) {
	my @solfile = split '\/',$solfile;
	$solfile = pop @solfile; # To return only filename, not path too
	return "$solfile";
    } else {
        print "\n              *** ERROR! ***\n";
        print " There is no .sol file for this instrument/date\n";
        print " Check ${path}/ for a .sol file\n";
        print " Make sure it's in the proper format:\n\n";
        print "  <site>_<type>_<yyyymmdd>_<hexb>_<hexh>_<sn>.sol\n\n";
	print " Also confirm you\'re not trying to process data prior";
	print " to first .sol file\n\n";
        die;
    }

}
1;


