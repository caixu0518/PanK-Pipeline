#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- AA_CFv5.candidate.pos.list
my $in1 = $ARGV[1]; ##- Bra_CFv5.fasta.sizes

my ($winsize, $stepsize) = (200000, 5000);

my $winfile = $in1.".w200k.s5k.txt";
   #system("bedtools makewindows -g   $in1   -w 200000  -s 5000 |  awk '{print \$1\"\t\"\$2+1\"\t\"\$3}' > $winfile");
   &size2bed($in1, $winsize, $stepsize, $winfile);


my %kmerindex = ();
   &readkmerPos($in0, \%kmerindex);

my $bedoutfile = $in0.".w200k.s5k.list";   
   &countkmersIntervals($winfile, \%kmerindex, $bedoutfile);
   system("rm -rf $winfile");
   #system("Rscript distributionPlot.r   $bedoutfile");


sub countkmersIntervals {
    
    my ($bedIn, $kmerIndex, $bedout) = @_;


    open (my $FR1, $bedIn);
    open (my $FO, ">$bedout");
    print $FO join("\t", "Chr", "Start", "End", "Count"), "\n";
    while(<$FR1>){
      chomp;
      my @temp     = split(/\t/, $_);
      my $countnum = 0;
   
      next, if($temp[0] !~ /^C/);  ##- it depends_

      for my $keys(keys %{$kmerIndex ->{$temp[0]}}){
          for my $keye(keys %{$kmerIndex ->{$temp[0]} ->{$keys}}){
              if($keys >= $temp[1] && $keye <= $temp[2]){
	         $countnum += 1;
	      }
	  }
      }
      #print $FO join("\t", @temp, $countnum), "\n"; 

      if($countnum > 0){
         print $FO join("\t", @temp, $countnum), "\n";
      }
      else{
        if($temp[1] == 1){      ##- first interval
      	   print $FO join("\t", @temp, $countnum), "\n";
      	}
      	if($temp[2]%100 != 0){  ##- last interval
      	   print $FO join("\t", @temp, $countnum), "\n";
      	}
      
      }
    }
    close ($FO);
    close ($FR1);

    ##---- calculate threshold (top 5%; top 10%)
    my @data  = ();
    open (my $FR000, $bedout);
    <$FR000>; 
    while(<$FR000>){
      chomp;
      my @temp = split(/\t/, $_);
      push(@data, $temp[3]), if($temp[3] >0);
    }
    close ($FR000);
    
    my @sorteddata  = sort {$b<=>$a} @data;
    my $top5index   = int((scalar(@sorteddata)/100)*5) -1;
    my $top10index  = int((scalar(@sorteddata)/100)*10) -1;

    my $top5  = $sorteddata[$top5index];
    my $top10 = $sorteddata[$top10index];

    print "top 5%  threshold: ", $top5, "\n";
    print "top 10% threshold: ", $top10, "\n";

}

sub size2bed {

    my ($fileIn, $winsize, $stepsize, $fileout) = @_;

    open (my $FR00, $fileIn);
    open (my $FO00, ">$fileout");
    while(<$FR00>){
      chomp;
      my @temp   = split(/\t/, $_);
      my $chrid  = $temp[0];
      my $chrlen = $temp[1];

      my $r=0;
      for($r=0; ($r*$stepsize+$winsize) <= $chrlen; $r++){
          print $FO00 join("\t", $chrid, $r*$stepsize+1, $r*$stepsize+$winsize), "\n"; 
      }
      if($r*$stepsize+$winsize > $chrlen && $r*$stepsize < $chrlen){
         print $FO00 join("\t", $chrid, $r*$stepsize+1, $chrlen), "\n";
      }
    }
    close ($FR00);
    close ($FO00);

}

sub readkmerPos {

    my ($FileIn, $kmerPosindex) = @_;

    open (my $FR0, $FileIn);
    while(<$FR0>){
      chomp;
      my @temp = split(/\t/, $_);
         $kmerPosindex ->{$temp[0]} ->{$temp[1]} ->{$temp[2]} = "Y";
    }
    close ($FR0);

}
