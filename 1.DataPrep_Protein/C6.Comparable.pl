#!/usr/local/bin/perl -w

# Debjit Ray
# 05/21/2012
# This program is for Normalizing the expression values by dividing them with the maximum value in the dataset

use strict;
use List::Util 'max';
use List::Util qw(sum);

my $input  = "7.62Reaction_Expression.txt";
my $output = "8.62Finale_Expression.txt";
my $greatest = 0;
my $Overall_Sum = 0;
my $count = 0;
my $length = 0;
my $Avg = 0;
my $factor = 0;
my @all = ();
my @summation = ();

open(IN, "< $input") || die "could not open $input\n";
while (<IN>) {
  /^GeneName/ and next;
  (my $line = $_) =~ s/\r*\n//;
  my($gene,@arr)=split(/\t/,$line);
  push(@all, max(@arr));
  push(@summation, sum(@arr));
  $length = scalar(@arr);
  $count++;
}
$greatest = max(@all);
$Overall_Sum = sum(@summation);
close(IN);

open(IN, "< $input") || die "could not open $input\n";
open (OUT , "> $output");
while (<IN>) {
  if ($_=~/^GeneName/) {
    print OUT "$_";
    next;
  }
  (my $line = $_) =~ s/\r*\n//;
  my(@arr) = split(/\t/,$line);

  $Avg = sprintf "%.5f",$Overall_Sum/($count*$length);
  $factor = sprintf "%.5f",1000/$Avg;

  print OUT $arr[0]."\t".sprintf("%.5f",$arr[1]*$factor)."\t".sprintf("%.5f",$arr[2]*$factor)."\t".sprintf("%.5f",$arr[3]*$factor)."\t".sprintf("%.5f",$arr[4]*$factor)."\t".sprintf("%.5f",$arr[5]*$factor)."\t".sprintf("%.5f",$arr[6]*$factor)."\t".sprintf("%.5f",$arr[7]*$factor)."\t".sprintf("%.5f",$arr[8]*$factor)."\t".sprintf("%.5f",$arr[9]*$factor)."\t".sprintf("%.5f",$arr[10]*$factor)."\t".sprintf("%.5f",$arr[11]*$factor)."\t".sprintf("%.5f",$arr[12]*$factor)."\t".sprintf("%.5f",$arr[13]*$factor)."\t".sprintf("%.5f",$arr[14]*$factor)."\t".sprintf("%.5f",$arr[15]*$factor)."\t".sprintf("%.5f",$arr[16]*$factor)."\t".sprintf("%.5f",$arr[17]*$factor)."\t".sprintf("%.5f",$arr[18]*$factor)."\n";
}
close(IN);
close(OUT);
print "\nTotal number of genes in Expression File\t".$count."\n";
print "The Greatest Value in the selected dataset is\t".$greatest."\n";
print "The Sum in the selected dataset is\t".$Overall_Sum."\n";
print "The Averaged in the selected dataset is\t".$Avg."\n";
print "The Multiplication factor to be used is\t".$factor."\n\n";

exit;


