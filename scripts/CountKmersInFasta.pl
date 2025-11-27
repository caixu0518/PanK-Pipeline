#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- Bra_PCF.kmers.txt
my $in1 = $ARGV[1]; ##- Bra_PCF.fasta
my $in2 = $ARGV[2]; ##- ksize
my $in3 = $ARGV[3]; ##- scripts path
my $in4 = $ARGV[4]; ##- outpre

my $k             = $in2;
my $scripts       = $in3;
my $filePartsNum  = 30;


##-step 1: run jellyfish 
my $jfFile    = $in4.".jf";
my $histoFile = $in4.".histo";
   if(! -e $jfFile || ! -e $histoFile){
      &RunJellfish($in1, $k, $jfFile, $histoFile);   
   }

##-step 2: devided query kmer file into small parts
my @fileParts    = ();
   &devidedIntoParts($in0, $filePartsNum, $in4, \@fileParts);
   
##-step 3: count target kmer
my $out = $in4.".kmerCount.list";
   &RunBatchCMD(\@fileParts, $jfFile, $out); 

##-clean
   system("rm  $jfFile");
   system("pigz $histoFile  $out");

##------------ all subs -----------------------------------------------------------
sub RunJellfish {

    my ($fasta, $ksize, $jfFile, $histoFile) = @_;
  
    if(! -e $jfFile){ 
       system("jellyfish  count  -m  $ksize  -t 20 -s 5G  -C  $fasta  -o  $jfFile");        ##- threads: 20 Mem: 5G
    }
    if(! -e $histoFile){
       system("jellyfish  histo  -t  20  $jfFile  >  $histoFile");                        ##- threads: 20
    }

}

sub RunBatchCMD {

    my ($fileParts, $jfFile, $out) = @_;

    my $count  = 0;
    my $outCMD = $jfFile.".batch.cmds";
    system("rm -rf $outCMD.completed"), if(-e  "$outCMD.completed");  

    open OUT0, ">$outCMD";
    for my $key(@{$fileParts}){
        $count += 1;
        print OUT0  "perl  $scripts/CalTargetKmerCounts.pl  $key  $jfFile  $key.out", "\n";
    }
    close OUT0;
    system("ParaFly -c $outCMD -CPU $count");

    system("rm -rf $out"), if(-e $out);
    for my $each(@{$fileParts}){
        my $outFile = $each.".out";
        system("cat  $outFile  >> $out");
        system("rm -rf $each  $outFile");
    }
    system("rm -rf  $outCMD  $outCMD.completed");

}

sub devidedIntoParts {

    my ($fileIn, $parts, $prefix, $files) = @_;

    my $partPara =  "l/".$parts;
       $prefix   =~ s/\.fasta//;
    my $outPre   =  $fileIn."_".$prefix;

    system("split  -n  $partPara  -d   $fileIn   $outPre");
    system("ls  $outPre*  >  $prefix.batch.inputs");

    my %indexFiles = ();
    open In0, "$prefix.batch.inputs";
    while(<In0>){
      chomp;
      my @fileInfo = split(/_/, $_);
      my $tmpindex = $_;
         $tmpindex =~ s/$outPre//;       
         $indexFiles{$tmpindex} = $_;
    }
    close In0;
    system("rm -rf $prefix.batch.inputs");
 
    for my $key(sort {$a<=>$b} keys %indexFiles){
        push(@{$files}, $indexFiles{$key});
    }

}
