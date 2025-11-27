#!/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;

#--Usage-----------------------------------------

my $usage=<<USAGE;
 
   The current script $0 is used to extract representative k-mers

   -species            [required]    Prefix for output files    i.e. Brapa  
   -ksize              [required]    k-mer size  i.e 17
   -pangenome          [required]    A file containing the abbreviation of each pan-genome member together with the absolute path to its genome sequence.
   -PipelinePath       [required]    The absolute path to the PanK-Pipeline directory   i.e. /mydata/caix/PanK-Pipeline 
   -PolymorphicKmer    [required]    rapa.merged.kmer.k17.list.Polymorphic_kmers.List

   Bug reports: Xu Cai
                caixu\@caas.cn

USAGE

my $currentPath = getcwd(); 

if (@ARGV == 0){die $usage}
my $in0; ##- species
my $in1; ##- pangenome
my $in2; ##- PipelinePath
my $in3; ##- kmer size
my $in4; ##- Polymorphic_kmers.List

GetOptions(
    "species:s"              =>\$in0,
    "ksize:i"                =>\$in3,
    "pangenome:s"            =>\$in1,
    "PipelinePath:s"         =>\$in2,
    "PolymorphicKmer:s"      =>\$in4,
);
if(not defined $in0 || not defined $in1 || not defined $in2 || not defined $in3 || not defined $in4){
   die $usage;
}

my $ksize   = $in3;
my $threads = 20;      ##- Number of threads (1)
my $memory  = "15G";   ##- Initial hash size


    &main(); 

##---- all subs --------------------------------------------------------
sub main {

    ##- check dependencies
    &checkdependence();
 
    ##- store index/fasta
    my @indexarray = ();
    my @fastaarray = ();
       &readPangenome($in1, \@indexarray, \@fastaarray);

    ##- query k-mer in each fasta
    my $fastaKmerList       = "FastaKmerList";
    my $fastaKmerListLOG    = "CleanKmerList.txt";
    if(not -d $fastaKmerList){
       system("mkdir  $fastaKmerList");
    }

    open (my $FO, ">$fastaKmerListLOG");
    for(my $m=0; $m<=$#indexarray; $m++){
        my $indexname   = $indexarray[$m];
	my $fastaname   = $fastaarray[$m];
        my $timestrings = &Times();

        if(not -e "$currentPath/$fastaKmerList/$indexname.k$ksize.clean.list.gz"){
	   ##- query k-mers in each assembled genome fasta
           my $cmdstring = "perl  $in2/scripts/FastaIntoCanonicalKmerWithCoords.pl  $fastaname  $ksize  $threads  $indexname";
	   print STDERR "[Current system time: $timestrings] $cmdstring \n";
	   system("$cmdstring");
           system("mv  $indexname.k$ksize.clean.list.stat.xls  $indexname.k$ksize.clean.list.gz   $fastaKmerList");    
        }
	else{
	   print STDERR "[Current system time: $timestrings] skip, as $indexname.k$ksize.clean.list.gz already exists.\n";   
	}

        ##- output kmer count file name
	if(not -e "$currentPath/$indexname.k$ksize.clean.list.gz"  &&  -e "$currentPath/$fastaKmerList/$indexname.k$ksize.clean.list.gz"){
	   system("ln -s $currentPath/$fastaKmerList/$indexname.k$ksize.clean.list.gz   .");
	}
        if(-e "$currentPath/$indexname.k$ksize.clean.list.gz"){
	   print $FO "$indexname.k$ksize.clean.list.gz", "\n";
	}
    }
    close ($FO);
 

    ##- select representative k-mers
    my $cmd         = "perl $in2/scripts/Select_Representativekmers.pl  $in4  $fastaKmerListLOG  $in2/scripts";
    my $timestrings = &Times();
    print STDERR "[Current system time: $timestrings] $cmd \n";
    system("$cmd");
    for(my $m=0; $m<=$#indexarray; $m++){
        my $indexname = $indexarray[$m];
    	if(-e "$currentPath/$indexname.k$ksize.clean.list.gz"){
    	   system("rm $indexname.k$ksize.clean.list.gz");
    	}
    }
    system("rm -rf $fastaKmerListLOG");

}

sub readPangenome {
  
    my ($filein, $indexarray, $fastaarray) = @_;

    ##- generate fastasize dir
    if(not -d "Fastasizes"){
       system("mkdir Fastasizes");
    }

    open (my $FR00, $filein);
    while(<$FR00>){
      chomp;
      my @temp = split(/\t/, $_);
      push(@{$indexarray}, $temp[0]);
      push(@{$fastaarray}, $temp[1]);
     
      my %id2seqtmp = ();
         &readFasta($temp[1], \%id2seqtmp);
 
      open (my $FOtmp, ">$temp[0].sizes");
      for my $key(sort keys %id2seqtmp){
          my $tmpseq = $id2seqtmp{$key};
	  print $FOtmp  join("\t", $key, length($tmpseq)), "\n";
      }
      close ($FOtmp);
      system("mv  $temp[0].sizes  Fastasizes");
        
    }
    close ($FR00);

}

sub readFasta {

    my ($in,$id2seq) = @_;

    open(my $SFR,$in);
    my $id;
    while($_=<$SFR>) {
      if(/^>([^\s^\n]+)\s*\n*/){
         $id = $1;
	 $id2seq->{$id} = "";
      }
      else{
         chomp;
	 $id2seq->{$id} .= $_;
      }
    }
    close($SFR);

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
