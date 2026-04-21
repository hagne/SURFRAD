#!/usr/bin/perl
#
#========================================================
# Examples on how to run...
#
#     perl raw2tu.pl <site> <sensor> <start_YYYYMMDD> <end_YYYYMMDD>
#     
#     perl raw2tu.pl gwn mfrsr 20100101 20100228
#
# Some example inputs.  Check site_attributes file for complete list.
#
# SURFRAD SITES (surfrad)
# Table Mountain -> tbl
# Desert Rock    -> dra
# Fort Peck      -> fpe
#
# Campaigns
# Rutland VT -> rut
# San Luis Valley CO -> slv 
# Wasco OR -> was
# Red Lake Playa AZ -> red
#
# Sensor types:  mfrsr or mfr10
#
# NEED TO (Some may have resolved after March 2020 update):
#
# 1. Sanity checks for: headIDs, loggerIDs, valid from dates in .sol files
# 2. Figure out a 1FF vs 11FF check for .sol files
# 3. Bug with "Processing..." dra_mfrsr 20151113.  Missing line
# 4. Fix sanity checks and traps in cal file reading
#========================================================
use strict;
use lib "/home/grad/mfrsr/bin";
use Modules::Check_Site;
use Modules::Check_Dates;
use Modules::Determine_Dates;
use Modules::Site_Info;
use Modules::Determine_Sol;
use Modules::Valid_Date_Range;

# Make sure it is user mfrsr running the program
my $user = getlogin();
if ($user !~ /^mfrsr/) {
    print "\n Derp! You are running as user: $user\n";
    print " This program needs to be run by user \'mfrsr\' in order";
    print " to create\n the directories and files with proper ownership.\n";
    die "\n";
}

# Ensure correct number of command line arguments
if (scalar@ARGV != 4) {usage();}

my $tli   = $ARGV[0];
my $type  = $ARGV[1];
my $start = $ARGV[2];
my $end   = $ARGV[3];

# Some initial checks to make sure inputs are valid
my $tli_type = join('_',$tli,$type); # Because Check_Site needs TLI_TYPE
&Check_Site($tli_type);              # Does TLI_TYPE exist?
&Check_Dates($start,$end);           # Make sure dates are valid
&Valid_Date_Range($start,$end,185);  # Ensure range does not exceed 185 days

# Input data are good.  Here we go!
my @dates = &Determine_Dates ($start,$end); # Returns array of dates (yyyymmdd)

my $lastdate = pop @dates;
push @dates, $lastdate;

# Need to get a few dates past the last date
my $cnt  = 1;
my $year = substr($lastdate,0,4);
my $mm   = substr($lastdate,4,2);
my $dd   = substr($lastdate,6,2);
if ($mm+1 == 13) {
    $year++;
    $mm = 1;
} else {
    $mm++;
}
my $finaldate = ($year*10000) + ($mm*100) + $dd;

# Now to get an array of a few dates past the lastdate.
# Remove first date so as not to duplicate a value.
# Keep the next fourteen dates.
#
my @extra = &Determine_Dates ($lastdate,$finaldate);
shift @extra;    
my $numextra = scalar @extra;
foreach (1 .. ($numextra - 14)) {pop @extra;}

# Copy raw files for date range requested into temp directory
#
my $tdir = "/tmp/temp.mfrsr_processing";
`rm -rf $tdir`;
`mkdir $tdir`;

print "\nCopying \"${tli_type}\" raw files ${start}-${end} to $tdir\n\n";

# Gather site attributes from site_attributes file
#
my %info = &Site_Info($tli_type);   # Returns lat, lon, tz, etc.
my $dir  = $info{DIR}; # directory or "category"
my $tz   = $info{TZ};  # time zone, used with rsrsplit below...
my $path = "/data/Inst/MFR/${dir}/${tli}/${type}";  # Common path

foreach (@dates,@extra) {
    my $yyyy = substr($_,0,4); # Need to add year to $path
    `/bin/cp ${path}/raw/${yyyy}/${tli}_${type}_${_}\*\.xmd ${tdir}/`;
}

# Run rsrsplit on those files
#
chdir "$tdir";
`rsrsplit -z $tz -g 1439 -i ${tli}_${type}*.xmd`;

# Run tu w/correct options on each mtm file
# 
$tz = $tz * -1; # Time zone sign change for tu below (+ west of GMT)
foreach (@dates) {
    my $solfile = &Determine_Sol ($tli,$type,$_);
    print " ---> Using cosine file: $solfile\n";
    my ($hexb,$hexh,$sn) = (split /_/,$solfile)[3,4,5];
    $sn = (split '\.',$sn)[0];

    my $ofile = "${tli}_${type}_${_}_${hexb}_${hexh}_${sn}.tu";
    my $YYMMDD = substr($_,2,6);
    my $solfile = "/home/grad/mfrsr/Calinfo/SOL/$solfile";

    chomp (my $f = `ls ${tdir}/*.${YYMMDD}.mtm`);
    if (-s $f) {
	print "Processing $f  $_\n";

	chomp (my $ph  = `ph $f | grep Unit`);
	my $bID = (split /\$/,$ph)[1];
	if ($hexb ne $bID) {
	    print "   ==> BOARD ID \"$hexb\" from:\n";
	    print "       $solfile\n";
	    print "       Does not match ID \"$bID\" from raw file $f\n";
	    die "\n";
	}

	chomp (my $ph  = `ph $f | grep Head`);
	my $hID = (split /\$/,$ph)[1];
	if ($hexh ne $hID) {
	    print "   ==> HEAD ID \"$hexh\" from:\n";
	    print "       $solfile\n";
	    print "       Does not match ID \"$hID\" from raw file $f\n";
	}

# Process MFRSR or MFR10 data with tu.  Checks the "Flags" field of 
# the .mtm files using ph.  20 is shadowband enabled so MFRSR, 00 is 
# shadowband disabled so MFR10.  Note that time_zone's ($tz) sign
# was changed above as required by tu (west of GMT tu wants poistive
# values.
#
	if ($type eq "mfrsr") {  
	    chomp (my $ph  = `ph $f | grep Flags`);
	    my $flags = (split /\$/,$ph)[1];
#	    if ($flags != 20) {
	    if ($flags !~ /20|21|34/) {  # Added 21 on 20230613 to deal with
                                         # halted instruments.
		print "\n==> This is supposed to be a MFRSR .xmd file, but it seems\n";
		print "    to be something else based on \"$ph\"\n\n";
		print "    The file in question is: $f\n\n";
		die;
	    }
	    `tu -d joe -r start -z $tz -x maez -c $solfile -H -o $ofile $f > /dev/null 2>&1`;
	} else {   
	    chomp (my $ph  = `ph $f | grep Flags`);
	    my $flags = (split /\$/,$ph)[1];
	    if ($flags != 0) {
		print "\n==> This is supposed to be a MFR10 .xmd file, but it seems\n";
		print "    to be something else based on \"$ph\"\n\n";
		print "    The file in question is: $f\n\n";
		die;
	    }
	    `tu -d joe -r start -z $tz -x maez -H -o $ofile $f > /dev/null 2>&1`;
#	    `tu -d joe -r start -x maez -H -o $ofile $f > /dev/null 2>&1`;
	}
    }
}

print "\nMove *.tu ${path}/tu/<YYYY>\n\n";

my @tulist = `ls *.tu`; # Need listing of files to extract year
foreach my $tuf (@tulist) {  # Foreach tu file...
    chomp($tuf);
    my $yyyy = substr($tuf,10,4);

# Move .tu files to regular directory structure
    if (! -e "${path}/tu/${yyyy}/") {   # Dir exist?
	`mkdir -p ${path}/tu/${yyyy}/`; # Create if necessary
    }
    `mv $tuf ${path}/tu/${yyyy}`;
}

# Clean up
`rm -rf $tdir`;


# Subroutines below
#
sub usage {
    use strict;
    print "\n  To run, type \"perl raw2tu.pl site sensor begin_date end_date\"\n\n";
    print "  Site should be of TLI format, e.g. bnd dra fpe... \n";
    print "  Reference site_attributes file for valid list of sites\n";
    print "  Dates should be in YYYYMMDD format\n\n";
    die;
}
1;
