#!/usr/bin/perl -w
#
# gets version from $Makefile and then creates a new header file $header
# with the version information and a date stamp.
# if the makefile contains VERSION=1.0.0d (i.e., ends with a 'd' then the
#    version is assumed to be a development release and the time stamp in
#    the $header file is the current time otherwise it gets the last this is a development
#    version and the date becomes 
#
my $makefile = '../Makefile';  #want the parent Makefile
my $datefile = 'version.h';
my $header   = 'version.c';
my $version  = '0-0-0';

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings);

open IN, $makefile or die "could not open <$makefile>";
while (<IN>){
	last if m/VERSION\s*=/;
}
if (/=\s*([-0-9]+d?)/) {$version=$1};
close IN;

open OUT, ">$header" or die "could not open <$header> for output";

my $time = '';
if ($version =~ /d/) {
	($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	$time = sprintf "%02d:%02d on ", $hour, $minute;
} else {
	($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime((stat $datefile)[9]);
}
$year = 1900 + $yearOffset;
$time .= "$dayOfMonth $months[$month] $year";

print OUT "/*";
print OUT "   AUTOMATICALLY generated version header by  <$0>\n";
print OUT "         The version number obtained from       <$makefile>\n";
if ($version =~ /d/) {
    print OUT "         The date is the last compile time.\n";
} else {
	print OUT "         The date is the last change to         <$datefile>\n";
}
print OUT "*/\n";
print OUT "char *Version = \"$version ($time)\";\n";
close OUT;
