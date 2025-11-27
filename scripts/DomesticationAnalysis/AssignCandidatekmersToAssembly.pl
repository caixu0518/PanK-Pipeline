#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;

#--Usage-----------------------------------------

my $usage=<<USAGE;

   The current script $0 is used to  assign the candidate k-mers to different assemblies

   -candidateKmers        [required] Candidate k-mer index, output from  Generate_enrichmentKmersInTargetGroup.pl
   -fastaKmers            [required] The file includes the filenames of the k-mer files corresponding to each assembly. These k-mer files are generated in Step 1 and are located in the 'FastaKmerList' directory.
   -representativeKmers   [required] Representative k-mer list file
   
   Bug reports: Xu Cai
                caixu\@caas.cn
		18,8,2025

USAGE

my $currentPath = getcwd();
if (@ARGV == 0){die $usage}

my $in0; ##- Zicaitai.candidates.txt
my $in1; ##- clean.list.filenames.txt
my $in2; ##- representative.kmers.list

GetOptions(

	"candidateKmers:s"              =>\$in0,
	"fastaKmers:s"                  =>\$in1,
	"representativeKmers:s"         =>\$in2,

);

if(not defined $in0 || not defined $in1 || not defined $in2){
   die $usage;
}


my $candidateKseqfile = $in0.".kmers.list";
   if(not -e $candidateKseqfile){
      &extractKmers($in0, $in2, $candidateKseqfile);
   }

my %targetKmers = ();
   &readtargetKmers($candidateKseqfile, \%targetKmers);

   &generateOut($in1, \%targetKmers);
    

##-- all subs -------------------------------------------------------   
sub extractKmers {

    my ($candidatefile, $rawkmers, $kmerout) = @_;

    my %posindex = ();
    open (my $FR0, $candidatefile);
    while(<$FR0>){
      chomp;
      my @temp = split(/\t/, $_);
      $posindex{$temp[1]} = "Y";
    }
    close ($FR0);

    my $count = 0;
    open (my $FR1, $rawkmers);
    open (my $FO, ">$kmerout");
    while(<$FR1>){
      $count += 1;
      if(exists $posindex{$count}){
         print $FO $_;
      }
    } 
    close ($FO);
    close ($FR1); 

}

sub generateOut {

    my ($fileList, $targetKmers) = @_;

    if(not -e "results"){
       system("mkdir results");
    } 
    else{
       system("rm -rf results");
       system("mkdir results");
    }

    my $totalKmers = scalar(keys %{$targetKmers});

    open (my $LOG, ">Kmers.sample.log");
    open (my $FR1, $fileList);
    while(<$FR1>){
      chomp;
      my $filename =  $_;
      my $samid    =  $filename;
         $samid    =~ s/\.k\S+\.clean\.list\.gz//;
  
      my $hitkmercount = 0;	  
      my $timestrings = &Times();
      print STDERR "[Current system time: $timestrings] .... start $samid .....\n";

      my $regionFile = $samid.".candidate.pos.list";	 
      open (my $FR00, "pigz -dc $filename | ");	
      open (my $FO00, ">$regionFile");
      while(<$FR00>){
        chomp;
	my @temp = split(/\t/, $_);
	if(exists $targetKmers ->{$temp[0]}){
	   $hitkmercount += 1;
	   my @infos = split(/;/, $temp[2]);
	   for my $eachtmp(@infos){
	       my @detials = split(/:/, $eachtmp);
	       print $FO00 join("\t", $detials[0], $detials[1]+1, $detials[1]+17), "\n";
	   }
	}
      }
      close ($FO00);
      close ($FR00);     

      system("mv  $regionFile  results");
      my $mappedratio = sprintf("%.2f", $hitkmercount/$totalKmers);  
      print $LOG join("\t", $filename, $samid, $totalKmers, $hitkmercount, $mappedratio), "\n";

    
    }
    close ($FR1);
    close ($LOG);


}   



sub readtargetKmers {

    my ($Filein, $targetKmers) = @_;

    open (my $FR0, $Filein);
    while(<$FR0>){
      chomp;
      my @temp = split(/\t/, $_);
      $targetKmers ->{$temp[0]} = "Y";
    }
    close ($FR0);

}

sub Times {

    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();

    $year += 1900;
    $mon  += 1;
    my $current_time = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year, $mon, $mday, $hour, $min, $sec);
    return($current_time);

}
