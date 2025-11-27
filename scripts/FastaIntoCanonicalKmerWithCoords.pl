#!/usr/bin/perl -w
use strict;
use warnings;
use threads;

##- update 17,8,2025

my $in0    = $ARGV[0]; ##- Bra_PCF.fasta
my $in1    = $ARGV[1]; ##- ksize
my $in2    = $ARGV[2]; ##- threads
my $outpre = $ARGV[3]; ##- Bra_PCF

my $k       = $in1;
my $threads = $in2;

my $outraw      = $outpre.".k$k.list";
my $outclean    = $outpre.".k$k.clean.list";

##--------------------------------------------------------------------------
my %id2seq = ();
   &readFasta($in0, \%id2seq);

##---- generate kmers ------------------------------------------------------
  &generateSeqKmers(\%id2seq, $k, $threads, $outraw); 

##---- remove redundance ---------------------------------------------------
  &cleanKmers($outraw, $outclean);
  system("rm -rf $outraw");
  system("pigz  $outclean");

##---- all subs ------------------------------------------------------------
sub cleanKmers {

    my ($fileIn, $fileOut) = @_;

    my %kmerindex = ();
    open In0Clean, $fileIn;
    while(<In0Clean>){
      chomp;
      my @temp = split(/\t/, $_);
      if($temp[0] =~ /N/){
         next;
      }
      else{
         my $coords = $temp[1].":".$temp[2];
         if(not exists $kmerindex{$temp[0]}){
            $kmerindex{$temp[0]} = $coords;
         }
         else{
            $kmerindex{$temp[0]} .= ";".$coords;
         }
      } 
    }
    close In0Clean;
 
    my %statrepeatedkmers = ();
    open OUTClean, ">$fileOut";
    for my $key(keys %kmerindex){
        my @array      = split(/;/, $kmerindex{$key});
        my $kmercounts = scalar(@array);
        $statrepeatedkmers{$kmercounts} += 1;
        print  OUTClean  join("\t", $key, scalar(@array), $kmerindex{$key}), "\n";
    }
    close OUTClean;
    %kmerindex = ();

    open OUTCleanStat, ">$fileOut.stat.xls";
    print OUTCleanStat join("\t", "KmerFreq", "KmerCount"), "\n";
    for my $key1(sort {$a<=>$b} keys %statrepeatedkmers){
        print OUTCleanStat join("\t", $key1, $statrepeatedkmers{$key1}), "\n";
    }
    close OUTCleanStat;
    %statrepeatedkmers = ();

}
  
sub generateSeqKmers {

    my ($id2seq, $k, $threads, $output) = @_;
  
    system("rm -rf $output"),  if(-e $output);

    my @chrs = sort keys %{$id2seq};
    my %finished = ();
    my $m = 1;
    for($m=1; ($m*$threads)<=$#chrs; $m++){

        my ($start, $end) = (($m-1)*$threads, $m*$threads);
        my @thr        = ();
        my @fileoutTmp = ();
        for(my $i=$start; $i<$end; $i+=1){
            $finished{$chrs[$i]} = "Y";
            my $tmpOut = $output.".a.".$chrs[$i].".txt";
            push(@fileoutTmp, $tmpOut);
            my $chrid  = $chrs[$i];
            my $tmpSeq = $id2seq ->{$chrid};
               $thr[$i] = threads->create(\&generateSeqKmerFile, $tmpSeq, $chrid, $k, $tmpOut);    
        }   
        for(my $i = $start; $i < $end; $i+=1){
            $thr[$i]->join;
        }
        for my $each(@fileoutTmp){
            my $timestrings = &Times();
	    print STDERR "[Current system time: $timestrings] .... start $each .....\n";
            system("cat $each >> $output");
	    system("rm  $each");
        }        
    }
    ##- process left
    my @left = ();
    for(my $m=0; $m<=$#chrs; $m++){
        if(not exists $finished{$chrs[$m]}){
           push(@left, $chrs[$m]);    
        }
    }  

    my @thr = ();
    my @lastBinFiles = ();    
    for(my $i=0; $i<=$#left; $i+=1){
        my $tmpOut = $output.".b.".$left[$i].".txt";
        push(@lastBinFiles, $tmpOut);      
        my $chrid  = $left[$i];
        my $tmpSeq = $id2seq ->{$chrid};
           $thr[$i] = threads->create(\&generateSeqKmerFile, $tmpSeq, $chrid, $k, $tmpOut);
    }
    

    for(my $i=0; $i<=$#left; $i+=1){
        $thr[$i]->join;
    }
     
    for my $each(@lastBinFiles){
	   my $timestrings = &Times();
	   print STDERR "[Current system time: $timestrings] .... start $each .....\n";  
           system("cat $each >> $output");
	   system("rm  $each");
    }  
 
}

sub generateSeqKmerFile {

    my ($seq, $chrid, $k, $fileOut) = @_;
     
    my @canonical_kmers = ();
    my @canonical_kmersInfo = ();
       &generate_canonical_kmers($seq, $k, \@canonical_kmers, \@canonical_kmersInfo);  

       open OUTX, ">$fileOut";
       for(my $m=0; $m<=$#canonical_kmers; $m++){
           print OUTX join("\t", $canonical_kmers[$m], $chrid, $canonical_kmersInfo[$m]), "\n";
       }
       close OUTX;

}

sub generate_canonical_kmers {

    my ($sequenceIn, $k, $canonical_kmers, $canonical_kmersInfo) = @_;

    my $sequence = uc($sequenceIn);

    for (my $i = 0; $i <= length($sequence) - $k; $i++) {
        my $kmer = substr($sequence, $i, $k);
        my $canonical = canonical_kmer($kmer);
        push(@{$canonical_kmers}, $canonical);
        push(@{$canonical_kmersInfo}, $i);
    }
}

sub canonical_kmer {

    my ($kmer) = @_;
    my $rev_comp      =  reverse $kmer;
       $rev_comp      =~ tr/ACGT/TGCA/;
    my @sortString    =  sort($kmer, $rev_comp);
       return($sortString[0]);

}

sub readFasta {

  my ($in,$id2seq) = @_;
  open(my $SFR,$in);

  my $id;
  while($_=<$SFR>) {
    if(/^>([^\s^\n]+)\s*\n*/) {
      $id = $1;
      $id2seq->{$id} = "";
    }
    else {
      chomp;
      $id2seq->{$id} .= $_;
    }
  }
  close($SFR);

}

sub Times {

    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();

    $year += 1900;
    $mon  += 1;
    my $current_time = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year, $mon, $mday, $hour, $min, $sec);
    return($current_time);

}
