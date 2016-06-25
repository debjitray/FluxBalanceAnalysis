#!/usr/local/bin/perl -w

# Debjit Ray
# 01/18/2013
# This program is for selecting the expression values for the 83 genes that are considered in our model

use strict;

my $inputfile      =  "2.SK1_ratio.txt";
my $inputfile2     =  "4.GeneList.txt";
my $outputfile     =  "5.Selected69Genes.txt";
my %hash        = ();
my $count1       = 0;
my $count2       = 0;
my $count3       = 0;

# Input file containing the main data
open(FDR,"<$inputfile") or die "Can't open $inputfile: $!\n";
foreach(<FDR>) {
    (my $line = $_) =~ s/\r*\n//;
    if ($line =~ /^\s*$/ or /^ORF/) {next;}
    my ($gene,@tokens) = split(/\t/,$line); 
    my $tokens=join("\t",@tokens);   
    $hash{uc($gene)}{$tokens}++;
    $count1++;
  }
close(FDR);

#Input file containing the names of genes of our interest
open(FDR2,"<$inputfile2") or die "Can't open $inputfile: $!\n";
#Output file
open(FDW,">$outputfile") or die "Can't open $outputfile: $!\n";
print FDW "GeneName\tT0\tT1\tT2\tT3\tT4\tT6\tT8\tT10\n";
foreach(<FDR2>) {
  chomp;
  /^\s*$/ and next;
  my ($genes,$num) = split(/\t/,$_);
   if (exists $hash{$genes}) {
     foreach my $i( keys %{$hash{$genes}}) {
	  print FDW $genes."\t".$i."\n";
	  $count2++;
	}
    }
  $count3++;
}
close(FDR2);
close (FDW);

print "\nNumber of genes in the main file\t".$count1."\n";
print "Number of unique genes in the model\t".$count3."\n";
print "Number of gene expression selected in the model\t".$count2."\n\n";
