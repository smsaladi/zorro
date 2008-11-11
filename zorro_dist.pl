#!/usr/bin/perl

use warnings;
use strict;

if(@ARGV != 2){
    die "Format:calc_dist.pl file Nseq\n"
}


my $File = shift(@ARGV);
my $Nseq = shift(@ARGV);

my $pars = "$File.trim.pars";
my $phyfile = "$File.trim.phy";
my $fafile = "$File.trim";


my $cmd = "cp $File $fafile";
system($cmd);

if(!-e "outfile"){
    $cmd = "echo >outfile";
    system($cmd);
}

$cmd = "rm -rf infile $phyfile";
system($cmd);


$cmd = "clustalw2 -convert -outfile=$phyfile -infile=$fafile -output=phylip >& /dev/null";
system($cmd);

if(!-e $phyfile){
    $cmd = "mv $phyfile* $phyfile";
    system($cmd);
}

$cmd = "rm -rf $File.protdist $pars $fafile";
system($cmd);

open(FOUT,">$pars");

print FOUT "$phyfile\n";
print FOUT "F\n";
print FOUT "$File.protdist\n";
print FOUT "2\n";
print FOUT "y\n";

close(FOUT);


$cmd = "protdist < $pars >& /dev/null";
system($cmd);




my $matrix = "$File.trim.mat";
open(FIN,"$File.protdist");
open(FOUT,">$matrix");
<FIN>;


my $tseq = 0;
while(<FIN>){
    chop;
    my @list = split(/\s+/,$_);
    if($tseq == 0){
	shift(@list);
    }
    if($list[0] eq ""){
	shift(@list);
    }
    $tseq += @list;
    if($tseq >= $Nseq){
	print FOUT "@list\n";
	if($tseq > $Nseq){
	    die "Incorrect PROTDIST output\n";
	}
	$tseq = 0;
    }
    else{
	print FOUT "@list ";
    }
}

close(FIN);
close(FOUT);
$cmd = "rm -rf $pars $phyfile $File.protdist";
system($cmd);




