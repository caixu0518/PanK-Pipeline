#!/usr/bin/perl -w
use strict;
use List::Util qw/shuffle/;
use Cwd;
use Getopt::Long;

#--Usage-----------------------------------------

my $usage=<<USAGE;

   The current script $0 is used to detected candidate domestication related k-mers between derived and control groups

   -deriviedGroup   [required]  The file contains the IDs of all members of the derived group.
   -vcfIn           [required]  The VCF file containing the genotypes of all representative k-mers across all individuals in the resequencing population 
   -output          [required]  The output file name

   Bug reports: Xu Cai
                caixu\@caas.cn
		18,8,2025

USAGE

my $currentPath = getcwd();
if (@ARGV == 0){die $usage}

my $in0; ##- Zicaitai.ids.txt      ##- target ids
my $in1; ##- AA.maf.0.05.vcf.gz    ##- vcf in
my $in2; ##- outfile

GetOptions(
 
   "deriviedGroup:s"              =>\$in0,
   "vcfIn:s"                      =>\$in1,
   "output:s"                     =>\$in2,
);

if(not defined $in0 || not defined $in1 || not defined $in2){
   die $usage;
}


my $randemRuns = 20;
my $cutoffup   = 0.8;
my $cutoffdown = 0.2;

my $posnum = 0;
   $posnum = &calcuateTotalPos($in1);

my %targetsamid  = ();
my %controlsamid = ();
   &readtargetSamid($in0, $in1, \%targetsamid, \%controlsamid, $randemRuns, $in2);


##---  all subs -------------------------------------------------------------------------------------------------
sub readtargetSamid {

    my ($Filein, $vcfFileIn, $targetsamid, $controlid, $randemRuns, $out) = @_;

    open (my $FO, ">$out");
    my $timestrings = &Times();
    my %givenSamid  = ();
    open (my $FR0, $Filein);
    while(<$FR0>){
      chomp;
      my @temp = split(/\t/, $_);
      $givenSamid{$temp[0]} = "Y";
    }
    close ($FR0);

    my %index2samid = ();
    my $postmpcpunt = 0;

    my $FR1; 
    if($vcfFileIn =~ /gz$/){
       open ($FR1, "pigz -dc $vcfFileIn | ");
    }
    else{
       open ($FR1, "$vcfFileIn");
    }

    while(<$FR1>){
      if(/^#CHROM/){
         chomp;
         my @samLine = split(/\t/, $_);
         for(my $m=9; $m<=$#samLine; $m++){
	     $index2samid{$m} = $samLine[$m];
	     if(exists $givenSamid{$samLine[$m]}){
	        $targetsamid ->{$samLine[$m]} = "Y";
	     }
	     else{
	        $controlid ->{$samLine[$m]} = "Y";
	     }
	 } 
         my @targetgroups = keys %{$targetsamid};
	 my @controlid    = keys %{$controlid};

	 $timestrings = &Times();
	 print STDERR "Total positions: ", $posnum, "\n";
	 print STDERR "[Current system time: $timestrings] .... start  .....\n";
	 print STDERR "The target group individuals: ",  scalar(@targetgroups), "\t", join(";", @targetgroups), "\n";
         print STDERR "The control group individuals: ", scalar(@controlid), "\t", join(";", @controlid), "\n";

	 $timestrings = &Times();
	 print STDERR "[Current system time: $timestrings] .... start index gt about the given vcf .....\n";
      } 
      elsif(/^##/){
         next;
      }
      else{
         chomp;

	 ##---------process --------------------------------------------------------------------------------------
         $postmpcpunt += 1;
         if($postmpcpunt%500 == 0){
	    $timestrings = &Times();
	    my $runprocessratio = sprintf "%.2f", ($postmpcpunt/$posnum)*100;
	    print STDERR "[Current system time: $timestrings] .... Process: $runprocessratio % .....\n";
	 } 
	 elsif($postmpcpunt == $posnum){
            print STDERR "[Current system time: $timestrings] .... Process: 100 % .....\n";		 
	 }  
         else{
	    sleep 0.000001;
	 }
         ###----------------------------------------------------------------------------------------------------

	 my @temp = split(/\t/, $_);

         my $targetinfo  = "NaN";
	 my @controlInfo = ();

	 my @targetGTs   = ();
	 my @controlGTs  = ();
         for(my $n=9; $n<=$#temp; $n++){
	    my $tmpsamname = $index2samid{$n};
	    if(exists $targetsamid ->{$tmpsamname}){
	       push(@targetGTs, $temp[$n]); 
	    }
	    else{
	       push(@controlGTs, $temp[$n]);
	    }
	 }        

	 ##- record target gt info 
	 my $targetgts   = 0;
	 my $targettotal = scalar(@targetGTs);
	 for my $each(@targetGTs){
	     if($each eq '0/0'){
	        $targetgts += 1;
	     }
	 }
         $targetinfo     = $targetgts.":".$targettotal;
         my $targetRatio = sprintf "%.2f", $targetgts/$targettotal;

	 ##- record control gt info    <= $randemRuns
         my $totalcontrol = scalar(@controlGTs);
         my $n            = $targettotal;
         if($targettotal > $totalcontrol){
	    $n = int(0.8*$totalcontrol);
	 }	 
      
	 my $maxgt = 0; 
	 for(my $i=1; $i<=$randemRuns; $i++){
	     my @shuffled_indexes = shuffle(0..$#controlGTs);
	     my @pick_indexes     = @shuffled_indexes[ 0 .. $n - 1 ];  
	     my @picks            = @controlGTs[ @pick_indexes ];

	     my $gtcounttmp = 0;
	     my $totaltmp = scalar(@picks);
	     for my $tmpeach(@picks){
	         if($tmpeach eq '0/0'){
		    $gtcounttmp += 1;
		 }
	     }  
	     if($maxgt < $gtcounttmp){
	        $maxgt = $gtcounttmp;
	     }   
             push(@controlInfo, $gtcounttmp.":".$totaltmp); 
	 }
         my $controlgtratiomax = sprintf "%.2f", $maxgt/$n;  

	 ##- print 
	 if($targetRatio >= $cutoffup && $controlgtratiomax <= $cutoffdown){
            my @inputs    = ($targetgts, $targettotal, $maxgt, $n);
            my $prefixtmp = $temp[0].".".$temp[1]; 
            my $pvalue    = &fisherTest(\@inputs, $prefixtmp);
	    print  $FO  $temp[0], "\t", $temp[1], "\t",  $pvalue, "\t",  $targetRatio, "\t", $controlgtratiomax, "\t", $targetinfo, "\t", join(";", @controlInfo), "\n";
	 }

      }
    }
    close ($FR1);
    close ($FO);

    $timestrings = &Times();
    print STDERR "[Current system time: $timestrings] .... finish index gt about the given vcf .....\n";
    %givenSamid  = ();

}




sub Times {
 
    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();

    $year += 1900;
    $mon  += 1;
    my $current_time = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year, $mon, $mday, $hour, $min, $sec);
    return($current_time);

}


sub fisherTest {

    my ($data, $outPre) = @_;

    my $outRscript = $outPre.".r";
    my $Rout       = $outRscript.".r.out";

    my $tmpA1 = $data ->[0];
    my $tmpA2 = $data ->[1];
    my $tmpB1 = $data ->[2];
    my $tmpB2 = $data ->[3];

    open  OUTRscript,  ">$outRscript";
    print OUTRscript   "sink(\"$Rout\")", "\n";
    print OUTRscript   "x <- matrix(c($tmpA1, $tmpA2, $tmpB1, $tmpB2), ncol=2, nrow=2)", "\n";
    print OUTRscript   "fisher.test(x)\$p.value", "\n";
    print OUTRscript   "sink()", "\n";
    close OUTRscript;
    system("Rscript  $outRscript");

    my $pvalue = "NA";
    open INX, $Rout;
    while(<INX>){
       chomp;
       my @temp   = split(/\s+/, $_);
          $pvalue = $temp[1];
    }
    close INX;

    system("rm  $outRscript  $Rout");
    return($pvalue);

}

sub calcuateTotalPos {
   
    my ($fileIn) = @_;

    my $FR00;
    my $count = 0;
    if($fileIn =~ /gz$/){
       open ($FR00, "pigz -dc $fileIn | ");
    }
    else{
       open ($FR00, "$fileIn");
    }
    while(<$FR00>){
      next, if(/^#/);
      $count += 1;
    }
    close ($FR00);
    return($count);

}
