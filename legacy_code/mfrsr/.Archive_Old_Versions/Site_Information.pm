# Subroutine site_information gathers information about each site such as
# latitude, longitude, time zone, and whether it is a MFRSR or Tower inst.
#
# "new_date" is the date the new style head went into operation
#
sub Site_Information () {
    my $site = shift(@_);
    my $type = shift(@_);
    my $file = "/home/grad/mfrsr/Calinfo/site_attributes";
    open(SI, $file) or die "Can't open $file: $!\n";

    my @info;
    while (<SI>) {
	chomp;
	next if /^\s*$/;     # Skip blank lines                             
	next if /^\s*\#/;    # Skip comment lines                           
	s/\#.*$//;           # Strip tag comment                            
	if (/^${site}_${type}/) {
	    @info = split,$_;
	}
    }
    close (SI);
    if (scalar @info == 0) {
	print "\n              *** ERROR! ***\n";
	print " There is no site attribute information for this instrument\n";
	print " Check the $file file \n\n";
	die;
    }
    shift @info; # Strip off site_type, leaving active/lat/lon/tz/new_date
    return @info;
    
}
1;
