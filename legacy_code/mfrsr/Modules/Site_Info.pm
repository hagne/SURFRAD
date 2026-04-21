############################################################
# Subroutine site_information gathers information about each site such as
# latitude, longitude, time zone, and whether it is a MFRSR or Tower inst.
#
# "new_date" is the date the new style head went into operation
#
use strict;
use warnings;

sub Site_Info () {
    my $tli_type = shift(@_);
    my ($tli,$type) = split /_/,$tli_type;
    my $file = "/home/grad/mfrsr/Calinfo/site_attributes";
    open(SI, $file) or die "Can't open $file: $!\n";

    my %info = ();
    my $dir;


# regex variables for site_attributes file
    my $a = qr/[a-z0-9]{3}_[mfrsr|mfr10]{5}/;  # tli_type
    my $b = qr/A|I|P/;                         # Active, Inactive, or Process
    my $c = qr/-?[\d|\.]+/;                    # Lat/Lon
    my $d = qr/-?\d+/;                         # Time Zone
    my $e = qr/\d{8}/;                         # YYYYMMDD
    my $f = qr/.+/;                            # Long name/description

    while (<SI>) {
	chomp;
	next if /^\s*$/;     # Skip blank lines                             
	next if /^\s*\#/;    # Skip comment lines                           
	s/\#.*$//;           # Strip tag comment                            
	if (/^SURF/) {$dir = "SURFRAD";  next;}
	if (/^SOLR/) {$dir = "SOLRAD";   next;}
	if (/^Camp/) {$dir = "Campaign"; next;}
	if (/^Proj/) {$dir = "Project";  next;}
	if (/^${tli_type}/) {
	     $_ =~ m/^($a)\s+($b)\s+($c)\s+($c)\s+($d)\s+($e)\s+\"($f)\"/;
	     %info = (
		 DIR => $dir,
		 T_T => $1,
		 ACT => $2,
		 LAT => $3,
		 LON => $4,
		 TZ  => $5,
		 NEW => $6,
		 LOC => $7
		 );
	} 
    }
    close (SI);
    return(%info);
}
1;
__END__






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
