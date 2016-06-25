#!/usr/local/bin/perl -w

# Debjit Ray
# 05/22/2012
# This program is for taking care of reactions with OR, AND and AVERAGE logic

use strict;

use List::Util 'min';
use List::Util 'max';
use List::Util qw(sum);

sub mean { return @_ ? sum(@_) / @_ : 0 };

my %new_hash;

my $inputfile1      = "6.Selected69Genes_normalized.txt";
my $inputfile2      = "3.grRules.txt";
my $outputfile      = "7.62Reaction_Expression.txt";

open(FDR1,"<$inputfile1") or die "Can't open $inputfile1: $!\n";
foreach(<FDR1>){
  (my $line = $_) =~ s/\r*\n//;
  my ($gene,@expression) = split(/\t/,$line);
  $new_hash{$gene} = [@expression];
}
close(FDR1);

open(FDR2,"<$inputfile2") or die "Can't open $inputfile2: $!\n";
open(FDW,">$outputfile") or die "Can't open $outputfile: $!\n";
my $m=1;
my $Rname='V_';
my $counter=0;

foreach(<FDR2>) {
  (my $line = $_) =~ s/\r*\n//;
  $line =~ s/'//g;
  $line =~ s/\s+//g;
# FOR TAKING THE MAXIMUM VALUE FOR REACTIONS WITH GENES IN "OR" LOGIC
  if ($line =~ m/or/) {
    my @AE=(),my @AF=(),my @AG=(),my @AH=(),my @AI=(),my @AJ=(),my @AP=(),my @AQ=(),my @AR=(),my @AS=(),my @AT=(),my @AU=(),my @AV=(),my @AW=(),my @AX=(),my @AY=(),my @AZ=(),my @BA=();
    my (@genename) = split(/or/,$line);
    foreach my $i (@genename){
      push(@AE, @{$new_hash{$i}}[0]);
      push(@AF, @{$new_hash{$i}}[1]);
      push(@AG, @{$new_hash{$i}}[2]);
      push(@AH, @{$new_hash{$i}}[3]);
      push(@AI, @{$new_hash{$i}}[4]);
      push(@AJ, @{$new_hash{$i}}[5]);

      push(@AP, @{$new_hash{$i}}[6]);
      push(@AQ, @{$new_hash{$i}}[7]);
      push(@AR, @{$new_hash{$i}}[8]);
      push(@AS, @{$new_hash{$i}}[9]);
      push(@AT, @{$new_hash{$i}}[10]);
      push(@AU, @{$new_hash{$i}}[11]);
      push(@AV, @{$new_hash{$i}}[12]);
      push(@AW, @{$new_hash{$i}}[13]);
      push(@AX, @{$new_hash{$i}}[14]);
      push(@AY, @{$new_hash{$i}}[15]);
      push(@AZ, @{$new_hash{$i}}[16]);
      push(@BA, @{$new_hash{$i}}[17]);
    }
    print "There is a 'OR' logic for reaction\t".$Rname.$m."\n";
    print FDW $Rname.$m."\t".max(@AE)."\t".max(@AF)."\t".max(@AG)."\t".max(@AH)."\t".max(@AI)."\t".max(@AJ)."\t".max(@AP)."\t".max(@AQ)."\t".max(@AR)."\t".max(@AS)."\t".max(@AT)."\t".max(@AU)."\t".max(@AV)."\t".max(@AW)."\t".max(@AX)."\t".max(@AY)."\t".max(@AZ)."\t".max(@BA)."\n";
  }

# FOR TAKING THE MINIMUM VALUE FOR REACTIONS WITH GENES IN "AND" LOGIC
  elsif ($line =~ m/and/) {
    my @AE=(),my @AF=(),my @AG=(),my @AH=(),my @AI=(),my @AJ=(),my @AP=(),my @AQ=(),my @AR=(),my @AS=(),my @AT=(),my @AU=(),my @AV=(),my @AW=(),my @AX=(),my @AY=(),my @AZ=(),my @BA=();
    my (@genename) = split(/and/,$line);
    foreach my $i (@genename){
      push(@AE, @{$new_hash{$i}}[0]);
      push(@AF, @{$new_hash{$i}}[1]);
      push(@AG, @{$new_hash{$i}}[2]);
      push(@AH, @{$new_hash{$i}}[3]);
      push(@AI, @{$new_hash{$i}}[4]);
      push(@AJ, @{$new_hash{$i}}[5]);

      push(@AP, @{$new_hash{$i}}[6]);
      push(@AQ, @{$new_hash{$i}}[7]);
      push(@AR, @{$new_hash{$i}}[8]);
      push(@AS, @{$new_hash{$i}}[9]);
      push(@AT, @{$new_hash{$i}}[10]);
      push(@AU, @{$new_hash{$i}}[11]);
      push(@AV, @{$new_hash{$i}}[12]);
      push(@AW, @{$new_hash{$i}}[13]);
      push(@AX, @{$new_hash{$i}}[14]);
      push(@AY, @{$new_hash{$i}}[15]);
      push(@AZ, @{$new_hash{$i}}[16]);
      push(@BA, @{$new_hash{$i}}[17]);
    }
    print "There is a 'AND' logic for reaction\t".$Rname.$m."\n";
    print FDW $Rname.$m."\t".min(@AE)."\t".min(@AF)."\t".min(@AG)."\t".min(@AH)."\t".min(@AI)."\t".min(@AJ)."\t".min(@AP)."\t".min(@AQ)."\t".min(@AR)."\t".min(@AS)."\t".min(@AT)."\t".min(@AU)."\t".min(@AV)."\t".min(@AW)."\t".min(@AX)."\t".min(@AY)."\t".min(@AZ)."\t".min(@BA)."\n";
  }
# FOR TAKING THE AVERAGE VALUE FOR REACTIONS WITH GENES IN "AVG" LOGIC
  elsif ($line =~ m/avg/) {
    my @AE=(),my @AF=(),my @AG=(),my @AH=(),my @AI=(),my @AJ=(),my @AP=(),my @AQ=(),my @AR=(),my @AS=(),my @AT=(),my @AU=(),my @AV=(),my @AW=(),my @AX=(),my @AY=(),my @AZ=(),my @BA=();
    my (@genename) = split(/avg/,$line);
    foreach my $i (@genename){
      push(@AE, @{$new_hash{$i}}[0]);
      push(@AF, @{$new_hash{$i}}[1]);
      push(@AG, @{$new_hash{$i}}[2]);
      push(@AH, @{$new_hash{$i}}[3]);
      push(@AI, @{$new_hash{$i}}[4]);
      push(@AJ, @{$new_hash{$i}}[5]);

      push(@AP, @{$new_hash{$i}}[6]);
      push(@AQ, @{$new_hash{$i}}[7]);
      push(@AR, @{$new_hash{$i}}[8]);
      push(@AS, @{$new_hash{$i}}[9]);
      push(@AT, @{$new_hash{$i}}[10]);
      push(@AU, @{$new_hash{$i}}[11]);
      push(@AV, @{$new_hash{$i}}[12]);
      push(@AW, @{$new_hash{$i}}[13]);
      push(@AX, @{$new_hash{$i}}[14]);
      push(@AY, @{$new_hash{$i}}[15]);
      push(@AZ, @{$new_hash{$i}}[16]);
      push(@BA, @{$new_hash{$i}}[17]);
    }
    print "There is a 'AVG' logic for reaction\t".$Rname.$m."\n";
    print FDW $Rname.$m."\t".mean(@AE)."\t".mean(@AF)."\t".mean(@AG)."\t".mean(@AH)."\t".mean(@AI)."\t".mean(@AJ)."\t".mean(@AP)."\t".mean(@AQ)."\t".mean(@AR)."\t".mean(@AS)."\t".mean(@AT)."\t".mean(@AU)."\t".mean(@AV)."\t".mean(@AW)."\t".mean(@AX)."\t".mean(@AY)."\t".mean(@AZ)."\t".mean(@BA)."\n";
  }
# FOR TAKING THE ACTUAL VALUES FOR REACTIONS WITH NO LOGIC 
  else{
    my ($gene) = split(/\n/,$line);
    print FDW $Rname.$m."\t".@{$new_hash{$gene}}[0]."\t".@{$new_hash{$gene}}[1]."\t".@{$new_hash{$gene}}[2]."\t".@{$new_hash{$gene}}[3]."\t".@{$new_hash{$gene}}[4]."\t".@{$new_hash{$gene}}[5]."\t".@{$new_hash{$gene}}[6]."\t".@{$new_hash{$gene}}[7]."\t".@{$new_hash{$gene}}[8]."\t".@{$new_hash{$gene}}[9]."\t".@{$new_hash{$gene}}[10]."\t".@{$new_hash{$gene}}[11]."\t".@{$new_hash{$gene}}[12]."\t".@{$new_hash{$gene}}[13]."\t".@{$new_hash{$gene}}[14]."\t".@{$new_hash{$gene}}[15]."\t".@{$new_hash{$gene}}[16]."\t".@{$new_hash{$gene}}[17]."\n";
  }
$m++;
$counter++;
}
close(FDR2);
close(FDW);

print "The total reaction rules in the grRules file\t".$counter."\n";


