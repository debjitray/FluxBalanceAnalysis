#!/usr/local/bin/perl -w

# Debjit Ray
# 06/28/2012

use strict;
use Statistics::Descriptive;

my $inputfile      = "2.Calculated_Z_All(GLPK).txt";
my $inputfile2      = "/Users/debjit/Desktop/MyModel/SubNetwork/Model2/3.Validation/3.DeficientReaction13.txt";
my $x = 0;
my %hash;
my %list;
my %new;
my $essen1=0;
my $non_essen1=0;
my $uncap=0;
# File containing the main data
open(FDR,"<$inputfile") or die "Can't open $inputfile: $!\n";
open(FDR2,"<$inputfile2") or die "Can't open $inputfile2: $!\n";
my $head;
my @header=();
foreach(<FDR>){
  (my $line = $_) =~ s/\r*\n//;
  if ($line =~ /^\t/) {
    $line=~s/WILD//;
    ($head,@header) = split(/\t/,$line);
  }
  next if ($line =~ /^\t/);
  my ($names,@tokens) = split(/\t/,$line); 
  for (my $i=0; $i<scalar(@tokens);$i++) {
    if ($tokens[$i]>=2) {
      $hash{$header[$i]}++;
    }
    elsif ($tokens[$i]<=-2) {
      $list{$header[$i]}++;
    }
  }   
  $x++;
}
close(FDR);

foreach(<FDR2>){
  (my $line = $_) =~ s/\r*\n//;
  $new{$line}++;
}
  
foreach my $key (keys %list)   
{
  if (exists $hash{$key}) {
    $essen1++;
    if (exists $new{$key})
      {
	$uncap++;
      }
  }
  else {
    $non_essen1++;
  }
}
print "OverLap of Reactions between >=2 and <=-2 (A&B): ".$essen1."\n";
print "OverLap of Reactions between >=2 and <=-2 and essential genes (A&B&C): ".$uncap."\n";

