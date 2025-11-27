#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- merged.kmers.list.Polymorphic_kmers.List
my $in1 = $ARGV[1]; ##- T01.k17.clean.list.gz
my $in2 = $ARGV[2]; ##- outpre

my %kmer2freqindex = ();
   &kmer2freq($in0, \%kmer2freqindex);

my $out =  $in2.".coordsorted.list";

   &output($in1, \%kmer2freqindex, $out);

##---  all subs ---------------------------------------------
sub output {

    my ($filein, $kmerfreqindex, $out) = @_;

    my $outtmp = $out.".tmp";
    open (my $FO, ">$outtmp");
    open (my $FR1, "pigz -dc $filein | ");
    while(<$FR1>){
      chomp;
      my @temp = split(/\t/, $_);
      if($temp[1] == 1){
         my @info = split(/:/, $temp[2]);
	 if(exists $kmerfreqindex ->{$temp[0]}){
	    print $FO join("\t", $info[0], $info[1], $temp[0], $kmerfreqindex ->{$temp[0]}), "\n";
	 }
      } 
    }
    close ($FR1);
    close ($FO);
   
    system("sort -k1,1 -k2,2n $outtmp > $out");
    system("rm -rf $outtmp");
    system("pigz $out");

}   

sub kmer2freq {

    my ($fileIn, $kmer2freqindex) = @_;

    open (my $FR0, $fileIn);
    while(<$FR0>){
      chomp;
      my @temp = split(/\t/, $_);
         $kmer2freqindex ->{$temp[0]} = $temp[2];
    }
    close ($FR0);

}
