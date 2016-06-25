#!/usr/local/bin/perl -w

# Debjit Ray
# 01/18/2013
# This program is for Normalizing the expression values by dividing them with the maximum value in the dataset

use strict;
use List::Util 'max';
use List::Util qw(sum);

my $input = "5.Selected69Genes.txt";
my $output = "6.Selected69Genes_normalized.txt";
my $greatest=0;
my $Overall_Sum=0;
my $count = 0;
my @all=();
my @summation=();

open(IN, "< $input") || die "could not open $input\n";
while (<IN>) {
  /^GeneName/ and next;
  (my $line = $_) =~ s/\r*\n//;
  my($gene,@arr)=split(/\t/,$line);
  push(@all, max(@arr));
  push(@summation, sum(@arr));
}
$greatest = max(@all);
$Overall_Sum = sum(@summation);
close(IN);

open(IN, "< $input") || die "could not open $input\n";
open (OUT , "> $output");
while (<IN>) {
  if ($_=~/^GeneName/) {
    print OUT "GeneName\tT1\tT2\tT3\tT4\tT6\tT8\tT10\n";
    next;
  }
  (my $line = $_) =~ s/\r*\n//;
  my(@arr)=split(/\t/,$line);
  print OUT $arr[0]."\t".sprintf("%.5f",$arr[2]/$greatest)."\t".sprintf("%.5f",$arr[3]/$greatest)."\t".sprintf("%.5f",$arr[4]/$greatest)."\t".sprintf("%.5f",$arr[5]/$greatest)."\t".sprintf("%.5f",$arr[6]/$greatest)."\t".sprintf("%.5f",$arr[7]/$greatest)."\t".sprintf("%.5f",$arr[8]/$greatest)."\n";
  $count++;
}
close(IN);
close(OUT);
print "\nTotal number of genes in Expression File\t".$count."\n";
print "The Greatest Value in the selected dataset is\t".$greatest."\n";
my $var = sprintf "%.5f", $Overall_Sum/$greatest;
print "The Sum of all expression values before normalization is\t".$Overall_Sum."\n";
print "The Sum of all expression values after normalization is\t".$var."\n\n";

exit;


