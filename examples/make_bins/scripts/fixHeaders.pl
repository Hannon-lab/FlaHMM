#!/usr/bin/perl -w
# Author: Susanne Bornelöv
#
use strict;

open IN,$ARGV[1] or die "Cannot open header file\n";

while (my $row = <IN>) {
	chomp $row;
	my @line = split "\t",$row;

	next if $line[0] eq $line[1];

	open IN2,$ARGV[0];
	open OUT,">$ARGV[0].tmp";
	while (my $in2 = <IN2>) {
		my $search = $line[1];
		my $replace = $line[0];
		
		$in2 =~ s/\Q$search\E/$replace/g;
		print OUT $in2;
	}
	close OUT;
	close IN2;
	system("mv $ARGV[0].tmp $ARGV[0]");
}
