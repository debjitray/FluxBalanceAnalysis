#!/usr/local/bin/perl -w

# Debjit Ray
# 05/21/2012
# This program is for creating a list of the genes that are in the model (mainly derived from the grRules.txt)

# use strict;

my %new_hash;

my $inputfile1      = "3.grRules.txt";
my $outputfile      = "4.GeneList.txt";

open(FDR1,"<$inputfile1") or die "Can't open $inputfile1: $!\n";
open(FDW,">$outputfile") or die "Can't open $outputfile: $!\n";

my $total = 0;
my $genes = 0;
my $OR = 0;
my $AND = 0;
my $AVG = 0;
my $OTHERS = 0;

foreach(<FDR1>){
  (my $line = $_) =~ s/\r*\n//;
  $line =~ s/'//g;
  $line =~ s/\s+//g;
  # FOR TAKING THE MAXIMUM VALUE FOR REACTIONS WITH GENES IN "OR" LOGIC
  if ($line =~ m/or/){
    my (@genename) = split(/or/,$line);
    foreach (@genename) {
      $new_hash{uc($_)}++;
      $genes++;
    }
    $OR++;
  }

  elsif ($line =~ m/and/){
    my (@genename) = split(/and/,$line);
    foreach (@genename) {
      $new_hash{uc($_)}++;
      $genes++;
    }
    $AND++;
  }

  elsif ($line =~ m/avg/){
    my (@genename) = split(/avg/,$line);
    foreach (@genename) {
      $new_hash{uc($_)}++;
      $genes++;
    }
    $AVG++;
  }
  
  else {
    $new_hash{uc($line)}++;
    $genes++;
    $OTHERS++;
  }
  $total++;
}
close(FDR1);
foreach (keys %new_hash) {
  print FDW $_."\t".$new_hash{$_}."\n";
}
close(FDW);

print "\nThe total gene rules in the grRules file\t".$total."\n";
print " The total reactions with OR logic \t".$OR."\n";
print " The total reactions with AND logic \t".$AND."\n";
print " The total reactions with AVERAGE logic \t".$AVG."\n";
print " The total reactions with NO logic \t".$OTHERS."\n";
print " The total genes in the Model \t".$genes."\n";
print " The total unique genes in the Model \t".scalar (keys %new_hash)."\n\n";
