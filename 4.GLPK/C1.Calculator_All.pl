#!/usr/local/bin/perl -w

# Debjit Ray
# 06/08/2012
# This code calculates the Z score for the 62 (Recation) Mutants with the wild type
use strict;
use Statistics::Descriptive;

my $inputfile      = "1.GLPK Based_All.txt";
my $outputfile     = "2.Calculated_Z_All(GLPK).txt";
my $x = 0;
# File containing the main data
open(FDR,"<$inputfile") or die "Can't open $inputfile: $!\n";
open(FDW,">$outputfile") or die "Can't open $outputfile: $!\n";

foreach(<FDR>){
  (my $line = $_) =~ s/\r*\n//;
  if ($line =~ /^V1/) {
    $line=~s/WILD//;
    print FDW "\t".$line."\n";
  }
  next if ($line =~ /^V1/);
  my @tokens = split(/\t/,$line); 
  my $Obj = pop(@tokens);
  my @b=();
  foreach (@tokens) {
    next if $_=~/NaN/;
    push(@b,$_);
  }
  print FDW $Obj."\t";
  for (my $i=0; $i<scalar(@tokens)-1;$i++) {
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@b);
    my $std  = $stat->standard_deviation();
    print FDW sprintf("%.5f",($tokens[$i]-$tokens[(scalar(@tokens)-1)])/$std)."\t";
  }  
  print FDW "\n"; 
  $x++;
}
close(FDR);
close (FDW);
print "\nThe number of Time Frames:".$x."\n\n";
