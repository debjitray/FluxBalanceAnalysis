#!/usr/local/bin/perl -w

#use strict;

# Debjit Ray
# 06/07/2010
# This program is for....

my $inputfile      = "1.SK1aalpha_meiosis_modified.txt";

my $outputfile     = "2.SK1_ratio.txt";

my %hash        = ();
my %hash_avg        = ();

my $d=0;
my $x=0;




# Input file containing the main data
open(FDR,"<$inputfile") or die "Can't open $inputfile: $!\n";

foreach(<FDR>){
  /^ORF/ and next;
  (my $line = $_) =~ s/\r*\n//;
  next if ($line =~ /^\s*$/ && /^ORF/);
  my ($gene,@tokens) = split(/\s+/,$line); 

  if (!exists $hash{uc($gene)}){
    $hash{uc($gene)}{DATA}=[@tokens];
    $hash{uc($gene)}{COUNT}++;
    $d++;
  }
  elsif (exists $hash{uc($gene)}){
    foreach my $i(0...$#tokens){
      ${$hash{uc($gene)}{DATA}}[$i]+=$tokens[$i];
      
    }
    $hash{uc($gene)}{COUNT}++;
    
  }
  $x++;

}
close(FDR);

open(FDW,">$outputfile") or die "Can't open $outputfile: $!\n";
foreach my $i(sort keys %hash){

  my @arr=map {$_/$hash{$i}{COUNT}} @{$hash{$i}{DATA}};

  my $str=join("\t",@arr);

  print FDW $i."\t".$str."\n";
}
close (FDW);

print " Total genes genes \t".$x."\n";
print " No of genes after removal of the duplicate genes \t".$d."\n";





