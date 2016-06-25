#!/usr/local/bin/perl -w

# Debjit Ray
# 06/23/2012
# This code calculates the Z >=2 and <=-2 for the 62 (Reactions) Mutants with the wild type
use strict;
use Statistics::Descriptive;

my $inputfile      = "2.Calculated_Z_All(GLPK).txt";
my $inputfile2      = "/Users/debjit/Desktop/MyModel/SubNetwork/Model2/3.Validation/3.DeficientReaction13.txt";
my $outputfile     = "2A.Analyzed_Z_All.txt";
my $outputfile2     = "2B.Unique_Z_All.txt";
my $x = 0;
my %hash;
my %list;
my $essen1=0;
my $non_essen1=0;
my $uncap=0;
# File containing the main data
open(FDR,"<$inputfile") or die "Can't open $inputfile: $!\n";
open(FDW,">$outputfile") or die "Can't open $outputfile: $!\n";
open(FDW2,">$outputfile2") or die "Can't open $outputfile2: $!\n";
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
  for (my $i=0; $i<scalar(@tokens);$i++) {
    if ($tokens[$i]>=2 or $tokens[$i]<=-2) {
      #print FDW sprintf("%.5f",$tokens[$i])."\t";
      print FDW $header[$i]."\t";
      $hash{$header[$i]}++;
    }
  }  
  print FDW "\n"; 
  $x++;
}
close(FDR);
close (FDW);
print "\nThe number of Time Frames:".$x."\n\n";

open(FDR2,"<$inputfile2") or die "Can't open $inputfile2: $!\n";
foreach(<FDR2>){
  (my $line = $_) =~ s/\r*\n//;
  $list{$line}++;
}

foreach my $key (keys %hash)   # Finding the Validated Reactions with the Deficeint Ones 
{
  if (exists $list{$key}) {
    print FDW2 $key." [ESSENTIAL REACTIONS]"."\n";
    $essen1++;
  }
  else {
    print FDW2 $key."\n";
    $non_essen1++;
  }
}


print FDW2 "\n<< UnCaptured Reactions>>\n";
foreach my $key (keys %list) # Finding the Validated Reactions with the Deficeint Ones 
{
  if (!exists $hash{$key}) {
    print FDW2 $key."\n";
    $uncap++;
  }
}

print "The number of Essential Reactions in small Model Captured  :".$essen1."\n";
print "The number of Essential Reactions in small Model UnCaptured  :".$uncap."\n";
print "The number of Other Reactions Captured:  ".$non_essen1."\n\n";
