#!/usr/local/bin/perl -w

# Debjit Ray
# 06/23/2012
# This code calculates separately Z >=2 and <=-2 for the 62 (Reactions) Mutants with the wild type
use strict;
use Statistics::Descriptive;

my $inputfile      = "2.Calculated_Z_All(GLPK).txt";
my $inputfile2      = "/Users/debjit/Desktop/MyModel/SubNetwork/Model2/3.Validation/3.DeficientReaction13.txt";
my $outputfile1     = "3A.Analyzed_Z_All(>=2).txt";
my $outputfile2     = "3B.Unique_Z_All(>=2).txt";
my $outputfile3     = "3C.Analyzed_Z_All(<=-2).txt";
my $outputfile4     = "3D.Unique_Z_All(<=-2).txt";
my $x = 0;
my %hash;
my %other;
my %list;
my $essen1 =0;
my $essen2 =0;
my $non_essen1 =0;
my $non_essen2 =0;
my $uncap1 =0;
my $uncap2 =0;

# File containing the main data
open(FDR,"<$inputfile") or die "Can't open $inputfile: $!\n";
open(FDW,">$outputfile1") or die "Can't open $outputfile1: $!\n";
open(FDW2,">$outputfile2") or die "Can't open $outputfile2: $!\n";
open(FDW3,">$outputfile3") or die "Can't open $outputfile3: $!\n";
open(FDW4,">$outputfile4") or die "Can't open $outputfile4: $!\n";
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
  print FDW $names."\t";
  print FDW3 $names."\t";
  for (my $i=0; $i<scalar(@tokens);$i++) {
    if ($tokens[$i]>=2) {
      #print FDW sprintf("%.5f",$tokens[$i])."\t";
      print FDW $header[$i]."\t";
      $hash{$header[$i]}++;
    }
    elsif ($tokens[$i]<=-2) {
      print FDW3 $header[$i]."\t";
      $other{$header[$i]}++;
    }
      
  }  
  print FDW "\n"; 
  print FDW3 "\n";
  $x++;
}
close(FDR);
close (FDW);
print "\nThe number of Time Frames:".$x."\n";

open(FDR2,"<$inputfile2") or die "Can't open $inputfile2: $!\n";
foreach(<FDR2>){
  (my $line = $_) =~ s/\r*\n//;
  $list{$line}++;
}

foreach my $key (keys %hash)   # Finding the Validated Reactions with the Deficeint Ones for >=2 Z Values
{
  if (exists $list{$key}) {
    print FDW2 $key." [ESSENTIAL REACTION]"."\n";
    $essen1++;
  }
  else {
    print FDW2 $key."\n";
    $non_essen1++;
  }
}

foreach my $key (keys %other) # Finding the Validated Reactions with the Deficeint Ones for <=-2 Z Values
{
  if (exists $list{$key}) {
    print FDW4 $key." [ESSENTIAL REACTION]"."\n";
    $essen2++;
  }
  else {
    print FDW4 $key."\n";
    $non_essen2++;
  }
}


print FDW2 "\n<< ESSENTIAL REACTION UNCAPTURED >>\n";
print FDW4 "\n<< ESSENTIAL REACTION UNCAPTURED >>\n";

foreach my $key (keys %list) # Finding the Validated Reactions with the Deficeint Ones for <=-2 Z Values
{
  if (!exists $hash{$key}) {
    print FDW2 $key."\n";
    $uncap1++;
  }
}

foreach my $key (keys %list) # Finding the Validated Reactions with the Deficeint Ones for <=-2 Z Values
{
  if (!exists $other{$key}) {
    print FDW4 $key."\n";
    $uncap2++;
  }
}
print "The number of Essential Reactions Captured for >=2 Z-Value  :".$essen1."\n";
print "The number of Essential Reactions UnCaptured for >=2 Z-Value  :".$uncap1."\n";
print "The number of Other Reactions Captured for >=2 Z-Value:  ".$non_essen1."\n";
print "The number of Essential Reactions Captured for <=-2 Z-Value:  ".$essen2."\n";
print "The number of Essential Reactions UnCaptured for <=-2 Z-Value:  ".$uncap2."\n";
print "The number of Other Reactions Captured for <=-2 Z-Value:  ".$non_essen2."\n\n";
