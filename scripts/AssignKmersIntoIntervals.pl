#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0];  ##- AA_CFv5.coordsorted.list.gz
my $in1 = $ARGV[1];  ##- AA_CFv5.fasta.sizes

my $out = $in0.".out";

my $winsize     = 10000;
my $winsizeFile = $in1.".windata.bed";
   &winfastasize($in1, $winsize, $winsizeFile);

##- calculate overlaps ------
my $tmpPosfile = $in0.".kmerpos.bed"; 
   `pigz -dc $in0 | awk  '{print \$1"\\t"\$2"\\t"\$2}'  >  $tmpPosfile`;
my $tmpPosfileintersect = $tmpPosfile.".intersect"; 
   `bedtools  intersect  -a  $tmpPosfile -b  $winsizeFile -wao  > $tmpPosfileintersect`;

my $KmerPosinIntervalFile = $winsizeFile.".add.Pos.list";
   &KmerPosinInterval($tmpPosfileintersect, $KmerPosinIntervalFile);

   system("rm -rf $winsizeFile  $tmpPosfile  $tmpPosfileintersect");

##----- all subs --------------------------------------------------------------------------------   
sub KmerPosinInterval {

    my ($fileIn, $fileout) = @_;

    my %indexpostmp = ();

    open (my $FR01, $fileIn);
    while(<$FR01>){
      chomp;
      my @temp = split(/\t/, $_);
      my $info = $temp[0].":".$temp[1];

      if(not exists $indexpostmp{$temp[3]}{$temp[4]}{$temp[5]}){
         $indexpostmp{$temp[3]}{$temp[4]}{$temp[5]} = $info;
      }
      else{
         $indexpostmp{$temp[3]}{$temp[4]}{$temp[5]} .= ";".$info;
      }
    }
    close ($FR01);

    open (my $FO01, ">$fileout");
    for my $chr(sort keys %indexpostmp){
	next, if($chr =~ /^\./);    
        for my $start(sort {$a<=>$b} keys %{$indexpostmp{$chr}}){
	    for my $end(sort {$a<=>$b}  keys %{$indexpostmp{$chr}{$start}}){
	        print $FO01 join("\t", $chr, $start, $end, $indexpostmp{$chr}{$start}{$end}), "\n";
	    }
	}
    }
    close ($FO01);

    %indexpostmp = ();

}

sub winfastasize {

    my ($fastaSize, $winsize, $fileout) = @_;

    if(-e $fastaSize){
   
       open (my $FO00, ">$fileout");	    
       open (my $FR00, $fastaSize);
       while(<$FR00>){
         chomp;
	 my @temp = split(/\t/, $_);

	 ##- generate intervals
	 my $m=1;
         my ($start, $end) = (0, 0);
	 for($m=1; $m*$winsize <=$temp[1]; $m++){
	     ($start, $end) = (($m-1)*$winsize+1, $m*$winsize);
	      print $FO00 join("\t", $temp[0], $start, $end), "\n";
	 }
         if(($m-1)*$winsize < $temp[1] && $m*$winsize > $temp[1]){
	     ($start, $end) = (($m-1)*$winsize+1, $temp[1]);
	     print $FO00 join("\t", $temp[0], $start, $end), "\n";
	 }
       } 	       
       close ($FR00); 
       close ($FO00);
    }
    else{
       die "cannot find: $fastaSize. \n";
    }

}
