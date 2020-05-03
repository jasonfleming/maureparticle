#!/usr/bin/env perl
use strict;
use warnings;

open IN, "<fort.14";


my $x0=101;  # center of vortes
my $y0=101;
my @K=(1..10);  # different proportianaly constants


my $line=<IN>;
$line=<IN>;
chomp ($line);

$line =~ s/^\s+//;
my ($ne,$np)=split(/\s+/,$line);

my @X;
my @Y;
foreach my $n (1..$np){
   $line=<IN>;
   chomp ($line);
   $line =~ s/^\s+//;
   my ($i,$x,$y)=split(/\s+/,$line);
   print "ixy $n $i $x $y\n";
   $X[$i]=$x;
   $Y[$i]=$y;
}
   
close IN;

my $nrecs = $#K +1;
open OUT, ">fort.64";
print OUT "An irrotatoinal Vortex\n";
print OUT "$nrecs $np 1000 1000 2\n";

my $t=1000;
my $dt=1000;

foreach my $k (@K){
   print OUT "$t $t\n";
   foreach my $n (1..$np){
      my $dx=$X[$n]-$x0;
      my $dy=$Y[$n]-$y0;
      my $r=( $dx**2 + $dy**2 )**0.5;
      my $vx=0;
      my $vy=0;
      if ($r > 0.0){
        my $speed = $k/$r;
        $vx=-$speed*$dy/$r;
        $vy=$speed*$dx/$r;
      }
      print OUT "$n $vx $vy\n";
     
   }
   $t+=$dt;
}
close (OUT);     
   




