#!/usr/bin/env perl

use strict;
use warnings;

use lib 'c:/myPerl';
use AdcircUtils::AdcGrid;


# get the name of the grid file (simple single command line argument)
my $gridFile=$ARGV[0];

unless (defined $gridFile){
   print "enter the name of the grid file (fort.14)\n"; 
   $gridFile=<>;
   chomp $gridFile;
}

print "Building neighbor Tables for $gridFile\n";

# load the grid
my $adcGrid=AdcGrid->new($gridFile);

# number of elements and nodes
my $ne=$adcGrid->{NE};
my $np=$adcGrid->{NP};


# get the node connectivity table
my @NOC=();

foreach my $i (1..$ne){
    my @noc=$adcGrid->getElement($i);
    $NOC[$i]=\@noc;  # note index zero will be  undef
}

######################################################
# determine the node2el table
my @node2el;
foreach my $i (1..$ne){
    foreach my $node (@{$NOC[$i]}){   
       push @{$node2el[$node]}, $i;   
    }  
}


# write the data
open OUT, ">node2el.tbl";
foreach my $i (1..$np){
   my $cnt=0;
   foreach my $nid (@{$node2el[$i]}){
      print OUT "$nid ";
      $cnt++;
   }
   while ($cnt < 12){
      print OUT "0 ";
      $cnt++;
   }
   print OUT "\n";
}
close OUT;





####################################################
# determie the el2el table
my @el2el;
print "generating el2el table\n";

foreach my $nid (1..$np){
    my @els=@{$node2el[$nid]};
    my @N2;
    my @N3;
 # print "nid, els is $nid @els\n";
    foreach my $el (@els){    #make a list of the elements with their nodes shifted so that $nid is the first node
         my @el_tmp=@{$NOC[$el]};
         while ($el_tmp[0] != $nid){
            my $n=shift(@el_tmp);
            push(@el_tmp, $n);
         }
         push @N2, $el_tmp[1];
         push @N3, $el_tmp[2];
  #   print "n2 @N2\n";
 #    print "n3 @N3\n";
     }
     foreach my $i (0..$#els){
         foreach my $j (0..$#els){
              if ($N2[$i]==$N3[$j]){
                 my $el1=$els[$i];
                 my $el2=$els[$j];
#         print "el1 el2 $el1 $el2\n";                
                 push(@{$el2el[$el1]},$el2);
                 push(@{$el2el[$el2]},$el1);
              }             
         }
     }
}         
        

# replace duplicates with zeros
foreach my $eid (1..$ne){
   my @els = @{$el2el[$eid]};
   
   foreach my $i (0..$#els){
        foreach my $k ($i+1..$#els){
             $el2el[$eid][$k]=0 if ($els[$i]==$els[$k]);   
        }
   }
}

# remove the zeros, except pad with zero to three elements
foreach my $eid (1..$ne){
   my @els = @{$el2el[$eid]};
   my @els2;
   foreach my $el (@els){
      push @els2, $el unless ($el ==0);
   }
   while ($#els2 < 2){
     push @els2, 0
   }
   $el2el[$eid]=\@els2;
}




# write the data with weirs
open OUT, ">el2el.tbl";
foreach my $i (1..$ne){
   foreach my $eid (@{$el2el[$i]}){
      print OUT "$eid ";
   }
   print OUT "\n";
}
close OUT;







# get the weir boundaries and add on weir neighbors
# have to add the weir matching elements after generating
# el2el otherwise you end up with an infinite loop 
# when generating the el3el table
my $nbou=$adcGrid->getVar('NBOU');
print "there are $nbou flux bnds\n";

foreach my $bnd (1..$nbou){
     my @list = $adcGrid->getFluxBnd($bnd);  # get open bnd number 1
     my $nvell=shift(@list);
     my $ibtype=shift(@list);
     my @nbvv=@{shift(@list)};
     if ($ibtype ==4 or $ibtype == 24){
        my @ibconn=@{shift(@list)} ;
        foreach my $n1 (@nbvv){
             my $n2=shift(@ibconn);
             my @n2els=@{$node2el[$n2]};
             my @n1els=@{$node2el[$n1]};
             foreach my $el (@n2els){
                 push @{$node2el[$n1]}, $el;
             }
             foreach my $el (@n1els){
                 push @{$node2el[$n2]}, $el;
             }
        }
     }
}



# write the data with weirs
open OUT, ">node2el_withWeirs.tbl";
foreach my $i (1..$np){
   my $cnt=0;
   foreach my $nid (@{$node2el[$i]}){
      print OUT "$nid ";
      $cnt++;
   }
   while ($cnt < 12){
      print OUT "0 ";
      $cnt++;
   }
   print OUT "\n";
}
close OUT;

        




