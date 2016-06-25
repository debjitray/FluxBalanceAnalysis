#!/usr/local/bin/perl -w

# Debjit Ray
# 05/16/2012
# This program is for selecting the expression data from B-J and 4-18 for the rib fp, it also considers only the gene name for the first column

use strict;

my $inputfile      = "1.GeneRPKM_allmrna,fps_notnormalized.txt";
my $outputfile     = "2.selected_18TF.txt";
my $x=0;
my %hash;
my @notunique;
# File containing the main data
open(FDR,"<$inputfile") or die "Can't open $inputfile: $!\n";
open(FDW,">$outputfile") or die "Can't open $outputfile: $!\n";

foreach(<FDR>){
  #/^description/ and next;
  (my $line = $_) =~ s/\r*\n//;
  next if ($line =~ /^\s/);
  $line =~ s/"//g;
  $line =~ s/fp RPKM//g;
  my ($first,@tokens) = split(/\t/,$line); 
  my ($genename) = split(/\s/,$first); 
  print FDW uc($genename)."\t".$tokens[29]."\t".$tokens[30]."\t".$tokens[31]."\t".$tokens[32]."\t".$tokens[33]."\t".$tokens[34]."\t".$tokens[40]."\t".$tokens[41]."\t".$tokens[42]."\t".$tokens[43]."\t".$tokens[44]."\t".$tokens[45]."\t".$tokens[46]."\t".$tokens[47]."\t".$tokens[48]."\t".$tokens[49]."\t".$tokens[50]."\t".$tokens[51]."\n";
  if ($hash{$genename}++ ) {
    push @notunique, $genename;
  }
  $x++;
}
close(FDR);
close (FDW);

foreach my $value(@notunique) {
  print $value."\n";
}
print "\nTotal number of non-unique genes \t".scalar(@notunique)."\n";
print "Total genes in the file \t".$x."\n\n";






