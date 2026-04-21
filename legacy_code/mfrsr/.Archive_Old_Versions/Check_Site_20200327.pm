#project Check_Site;

sub Check_Site () {
    use strict;

    my $file = "/home/grad/mfrsr/Calinfo/site_attributes";
    open(SI, $file) or die "Can't open $file: $!\n";

    if (scalar@_ == 1) {
	my $site = $_[0];
	
	my @info;
	while (<SI>) {
	    chomp;
	    next if /^\s*$/;     # Skip blank lines
	    next if /^\s*\#/;    # Skip comment lines
	    s/\#.*$//;           # Strip tag comment
	    if (/^${site}/) {
		@info = split,$_;
	    }
	}
	close (SI);
	
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
__END__
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
    
