#project Check_Site;
#
#####################################
# Module Check_Site.pm 
#
# Purpose is to do away with the need of including a site list, and the
# associated maintenance, in every program that deals with MFRSRs or
# MFR10s.  It was getting cumbersome to have to continually edit programs
# when a station/instrument was added or decomissioned.  Check_Site.pm
# uses a single source, site_attributes, to determine stations.  Now
# only the site_attributes file needs to be edited when an instrument
# is added or shut down.
#
# Inputs can be:
#    1: No arguments.  In this situation the module interrogates
#       the site_attributes file and returns the TLI_TYPE and 
#       associated directory foir each active station.  No arguments 
#       are most often used for the regular polling of all active
#       MFRSRs and MFR10s
#
#    2: With arguments in the form of TLI_TYPE.  Typically used if 
#       interested in collecting from one or more stations, maybe for 
#       a quick test.  Also to initialize a MFRSR or MFR10.
#
#    3: With arguments in the form of TLI.  Used to check the validity
#       of a station.  It provides a sanity check in case a station TLI
#       was entered incorrectly
#
# Example usage:
#
#    my @return_list = &Check_Site(@arg_list); # Used with poll.pl
#
#    &Check_Site(@arg_list)  # Check if TLIs are valid
#
# Known Shortcomings:
#
#  If input arguments are mixed types, e.g. TLIs and TLI_TYPEs, madness
#  will follow.
######################################
#
sub Check_Site () {
    use strict;
    use experimental 'smartmatch'; # Suppresses "~~" experimental warnings

# Open file with site attribute information
    my $file = "/mnt/mfrsr/Calinfo/site_attributes";
#    my $file = "/home/grad/Calinfo/site_attributes";
    open(SI, $file) or die "Can't open $file: $!\n";

# Sanity tests on inputs ensure validly formed tli or tli_type
    my $tli_flag = "N" ;
    if (scalar@_ > 0) {
	foreach my $site (@_) {
	    my $n = scalar(split /_/,$site); # Test if TLI or TLI_TYPE
	    if ($n == 2) {                   # TLI only, no underscore
		if ($site !~ /^[a-z0-9]{3}_mfrsr$|mfr10$/) {
		    &usage0($site,$file);
		    die "\n";
		}
	    } elsif ($n == 1) {
		if ($site !~ /^[a-z0-9]{3}$/) {
		    &usage0($site,$file);
		    die "\n";
		}
		$tli_flag = "Y";
	    } else {
		print "\n Sumpin\' weird be entered. Try again Homie!\n";
	    }
	}
    }

#######################################
# Here comes the meat of the matter
#
# First we read in the contents of site_attribute file with while loop
#
# @vlist contains only TLI_TYPEs listed as active in attributes file.
#        along with the appropriate directrory.  Used for MFRSR polling
# @tlist contains only TLI_TYPEs listed as active in attributes file.
#        without directory information.  Used for MFRSR polling
# @alist contains every site in attributes file.  Used to determine
#        if an entered tli is valid
#
    my @vlist; # Array of valid sites with associated directory
    my @tlist; # Array of valid TLI_TYPE, no associated directory
    my @alist; # Array of all sites
    my $dir;   # Directory determined from site_attributes
    while (<SI>) { # SI is file handle for site_attributes
	chomp;
	next if /^\s*$/;     # Skip blank lines
	next if /^\s*\#/;    # Skip comment lines
	s/\#.*$//;           # Strip tag comment
	if (/^SURF/) {$dir = "SURFRAD";  next;} #\
	if (/^SOLR/) {$dir = "SOLRAD";   next;} # \ These lines set appropriate
	if (/^Camp/) {$dir = "Campaign"; next;} # / directory names
	if (/^Proj/) {$dir = "Project";  next;} #/  
	if (/^[a-zA-Z]|[\d\d\d]/) {
	    my @info = split,$_;
	    if ($info[1] eq "A") { # If site is active...
		push (@tlist, $info[0]);
		my $line = sprintf "%s %s", $info[0],$dir;
		push (@vlist, $line);
	    }
	    my $line = sprintf "%s %s", $info[0],$dir;
	    push (@alist, $line);
	}
    }
    close (SI);

# Check if entered TLI_TYPE is an actual site. @tlist holds only valid TLI_TYPEs
    if ($tli_flag eq "N") {
	my @badarg = (); # Initialize array for invalid TLI_TYPEs
	foreach my $arg (@_) {
	    if ($arg ~~ @tlist) { # Is $arg "in" @tlist?
		next;
	    } else {
		push (@badarg, $arg);
	    }
	}
	if (scalar@badarg > 0) {
	    my $cnt = 0;
	    print "\n";
	    foreach my $ba (@badarg) {
		print " This TLI_TYPE not found: $ba\n";
		$cnt++;
	    }
	    &usage1($file,$cnt);  # If bad arguments print nasty gram and die
	    die "\n";
	}
    }

# Next section extracts all TLIs from @alist and returns unique TLIs.
# Does similar for TLI_TYPEs
    my @tli;      # Three letter identifier
    my @tli_type; # Three letter identifier w/instrument type
    foreach my $info (@alist) {
	my @site = split /\s+/, $info;
	push (@tli, (split /_/, $site[0])[0]); # Keep only TLI site name
	push (@tli_type, $site[0]);            # Keep TLI_TYPE
    }	
    @tli = do {my %seen; grep {!$seen{$_}++} @tli}; # Return Unique TLIs

# Here we pause regularly scheduled processing to ensure no duplicate
# entries exist in site_attributes file.
    my @verify = do {my %seen; grep {!$seen{$_}++} @tli_type}; # Unique tli_type
    if (scalar@verify - scalar@tli_type != 0) { # Test array sizes for dupe test
	print "\n Duplicate entry in $file\n";
	my %h = ();
	my $r;
	my $flag = 0;
	foreach $r (@tli_type) {    # Not entirely sure this foreach section
	    if (!exists($h{$r}))  { # is bullet proof.  Seems to work though.
		# First time we've seen this one
		$h{$r} = 0;
	    } elsif (exists $h{$r}) {
		# We've seen this one before and reported
		$h{$r}++;
		print " This site has a duplicate entry: $r\n";
		$flag = 1;
	    } else {
		continue;
	    }
	}
	%h = (); # Destroy %h sine we're done with it
	if ($flag == 1) {
	    die "\n Fix duplicate entries in $file\n\n";
	}
    }	

# Check if entered TLI is an actual site.
    if ($tli_flag eq "Y") {
	my @badarg = (); # Initialize array for invalid TLIs
	foreach my $arg (@_) {
	    if ($arg ~~ @tli) { # Is $arg "in" @tli?
		next;
	    } else {
		push (@badarg, $arg);
	    }
	}
	if (scalar@badarg > 0) {
	    my $cnt = 0;
	    print "\n";
	    foreach my $ba (@badarg) {
		print " This TLI not found: $ba\n";
		$cnt++;
	    }
	    &usage2($file,$cnt);  # If bad arguments print nasty gram and die
	    die "\n";
	}
    }

# Return all valid TLI_TYPE and DIRECTORY if @_ (passed as ARGV) is ZERO
    if (scalar@_ < 1) {return (@vlist);}

# Return command line TLI_TYPE with appropriate DIRECTORY if @_ > ZERO
    if ((scalar@_ > 0) && ($tli_flag eq "N")) { # Return nothing if only TLI
	my @returnlist;
	foreach my $clsite (@_) {         # Command line entered site(s)
	    foreach my $vline (@vlist) {  # Valid line from valid list
		my $vsite = (split /\s+/, $vline)[0]; # Extract only valid site
		if ($clsite eq $vsite) {
		    push (@returnlist, $vline);
		}
	    }
	}
	return (@returnlist);
    }
}


##################################
# Subroutines below...
#
sub usage0 () {
    my $site = shift(@_);
    my $file = shift(@_);
    print "\n Died due to the funky command line entry: $site\n";
    print " \n This module, Check_Site.pm, expects one of the following:\n";
    print "  1: No passed in arguments\n";
    print "  2: Arguments in the form of tli_type, e.g. tbl_mfrsr or bnd_mfr10\n";
    print "  3: Arguments in the form of tli\n";
    print "     tli  = Three Letter Identifier (lowercase).\n";
    print "            Also accepts numeric or alphanumeric and\n";
    print "            must be three characters in length\n";
    print "     type = mfrsr or mfr10 \(Note the underscore in between tli and type\)\n";
    print "\n Tips:\n";
    print "\n  \"tli_type\" must be found in: $file\n";
    print "  or another nasty gram will be issued.\n\n";
    print "  If passing in a site \"tli\" to check validity it must match a\n  tli found in: $file\n";
}


sub usage1 () {
    my $file = shift(@_);
    my $cnt  = shift(@_);
    
    if ($cnt == 1) {
	$p1 = "TLI_TYPE was";
	$p2 = "is not \n an active station.";
	$p3 = "station was";
    } else {
	$p1 = "TLI_TYPEs were";
	$p2 = "are not \n active stations.";
	$p3 = "stations were";
    }

    print "\n Entered $p1 either entered incorrectly, or $p2";
    print " Check: $file \n";
    print " If sure $p3 entered correctly change status to Active (A).\n";
}


sub usage2 () {
    my $file = shift(@_);
    my $cnt  = shift(@_);
    
    if ($cnt == 1) {
	$p1 = "TLI is";
	$p2 = "Three Letter Indentifyer";
	$p3 = "station was";
    } else {
	$p1 = "TLIs are";
	$p2 = "Three Letter Identifyers";
	$p3 = "stations were";
    }

    print "\n Entered $p1 not a valid $p2 for any known\n";
    print " station. ";
    print " Check: $file \n";
    print " If sure $p3 entered correctly edit site_attributes file\n";
    print " and add station.\n";
}


1;







__END__
    
    if (scalar@_ > 0) {
	my $site = $_[0];
	
	my @info;
	
	if (scalar @info == 0) {
	    print "\n                  *** ERROR! ***\n";
	    print " There is no site attribute information for ${site}.\n";
	    print " Check the $file file.\n\n";
	    die;
	}
    } else {
	
	
    }
    
}
1;

#    my @all_sites = qw/bnd dra fpe gwn psu sxf tbl lab slv rut was red/;
    my $site = shift(@_);
    
    if ($site !~ /bnd|dra|fpe|gwn|psu|sxf|tbl|lab|slv|rut|was|red|dsr|633|iss/) {
	die "\n  \"${site}\" does not match any site.  Check";
    }
    
    
#    my %site_codes = abbrev @all_sites;
    
#    my $tli = $site_codes{$site};
# || (input=~/^\d+$/ && $input<=$#all_sites && $all_sites[$input]);

#    if (! $tli) {
#	die "\n$site does not match any site.  Check";
#    }
}
1;
    
