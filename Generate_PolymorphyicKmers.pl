#!/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;

#--Usage-----------------------------------------

my $usage=<<USAGE;
 
   The current script $0 is used for constructing polymorphic k-mers using Pan-genome assemblies

   -species        [required]    Prefix for output files    i.e. Brapa  
   -ksize          [required]    k-mer size  i.e 17
   -pangenome      [required]    A file containing the abbreviation of each pan-genome member together with the absolute path to its genome sequence.
   -PipelinePath   [required]    The absolute path to the PanK-Pipeline directory   i.e. /mydata/caix/PanK-Pipeline 

   Bug reports: Xu Cai
                caixu\@caas.cn

USAGE

my $currentPath = getcwd(); 

if (@ARGV == 0){die $usage}
my $in0; ##- species
my $in1; ##- pangenome
my $in2; ##- PipelinePath
my $in3; ##- kmer size

my $threads = 20;      ##- Number of threads (1)
my $memory  = "15G";   ##- Initial hash size

GetOptions(
    "species:s"              =>\$in0,
    "ksize:i"                =>\$in3,
    "pangenome:s"            =>\$in1,
    "PipelinePath:s"         =>\$in2,
);
if(not defined $in0 || not defined $in1 || not defined $in2 || not defined $in3){
   die $usage;
}

my $mergedKmers = $in0.".merged.kmer.k$in3.list";
   &main(); 

##---- all subs --------------------------------------------------------
sub main {

    ##- check dependencies
    &checkdependence();
  
    ##- generate merged k-mer list 
    my @indexarray = ();
    my @fastaarray = (); 
    &generateMergedKmers($in1, $in0, $in3, $mergedKmers, \@indexarray, \@fastaarray);    

    ##- query k-mer in each fasta
    my $fastaCount     = "KmersInEachGenome";
    my $fastaKmercount = "kmercountFilesIneachGenome.txt";
    if(not -e $fastaCount){
       system("mkdir  $fastaCount");
    }

    open (my $FO, ">$fastaKmercount");
    for(my $m=0; $m<=$#indexarray; $m++){
        my $indexname = $indexarray[$m];
	my $fastaname = $fastaarray[$m];
        my $timestrings = &Times();
        if(not -e "$currentPath/KmersInEachGenome/$indexname.kmerCount.list.gz"){
	   ##- query k-mers in each assembled genome fasta
           my $cmdstring = "perl  $in2/scripts/CountKmersInFasta.pl  $mergedKmers  $fastaname  $in3  $in2/scripts  $indexname";
	   print STDERR "[Current system time: $timestrings] $cmdstring \n";
	   system("$cmdstring");
           system("mv  $indexname.histo.gz  $indexname.kmerCount.list.gz  $fastaCount");    
        }
	else{
	  print STDERR "[Current system time: $timestrings] skipt, as $indexname.kmerCount.list.gz already exists.\n";   
	}

        ##- output kmer count file name
	if(not -e "$currentPath/$indexname.kmerCount.list.gz"  &&  -e "$currentPath/KmersInEachGenome/$indexname.kmerCount.list.gz"){
	   system("ln -s $currentPath/KmersInEachGenome/$indexname.kmerCount.list.gz  .");
	}
        if(-e "$currentPath/$indexname.kmerCount.list.gz"){
	   print $FO "$indexname.kmerCount.list.gz", "\n";
	}
    }
    close ($FO);
  
    ##- select polymorphic k-mers
    my $cmd         = "perl $in2/scripts/Select_Polymorphic_kmers.pl  $fastaKmercount   $mergedKmers";
    my $timestrings = &Times();
    print STDERR "[Current system time: $timestrings] $cmd \n";
    system("$cmd");
    for(my $m=0; $m<=$#indexarray; $m++){
        my $indexname = $indexarray[$m];
	if(-e "$currentPath/$indexname.kmerCount.list.gz"){
	   system("rm $indexname.kmerCount.list.gz");
	}
    }
    system("rm -rf $fastaKmercount");

}

sub generateMergedKmers {
   
    my ($pangenome, $outpre, $ksize, $kmerlist, $indexarray, $fastaarray) = @_;

    my $tmpfaOut = $outpre.".tmp.fa";
    my @indexs   = ();
    my @fasta    = ();
    open (my $FR0, $pangenome);
    while(<$FR0>){
      chomp;
      my @temp = split(/\t/, $_);
      push(@indexs, $temp[0]);
      push(@fasta, $temp[1]);
    }
    close ($FR0);
    @{$indexarray} = @indexs;
    @{$fastaarray} = @fasta;
    
    my $timestrings = &Times();
    if(not -d  $mergedKmers){
       for my $eachfa(@fasta){
           system("cat $eachfa >> $tmpfaOut");
       }
       my $timestrings = &Times();
       print STDERR "[Current system time: $timestrings] jellyfish  count  -m  $ksize  -t $threads  -s  $memory  -C   $tmpfaOut -o  $outpre.tmp.jf ... \n";
       system("jellyfish  count  -m  $ksize  -t $threads  -s  $memory  -C   $tmpfaOut -o  $outpre.tmp.jf");
       $timestrings = &Times();
       print STDERR "[Current system time: $timestrings] jellyfish  dump -c -t $outpre.tmp.jf ... \n";
       system("jellyfish  dump -c -t $outpre.tmp.jf > $kmerlist");
       system("rm -rf $outpre.tmp.jf  $tmpfaOut"); 
    }
    else{
      print STDERR "[Current system time: $timestrings] skip, as $mergedKmers already exists.\n"; 
    }
}

sub checkdependence {

    ##- check jellyfish
    my $timestrings = &Times();
    print STDERR "[Current system time: $timestrings] Check if the system has installed the jellyfish software ... \n";
    my $jellyfishinfo = `jellyfish  --help`; 
    if($jellyfishinfo =~ /Usage:\s+jellyfish/){
       $timestrings = &Times();
       print STDERR "[Current system time: $timestrings] The jellyfish software has been detected.\n";
    }
    else{
       $timestrings = &Times();	     
       die "[Current system time: $timestrings] Please install jellyfish first.\n";
    }

}

sub Times {

    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();

    $year += 1900;
    $mon  += 1;
    my $current_time = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year, $mon, $mday, $hour, $min, $sec);
    return($current_time);
}
