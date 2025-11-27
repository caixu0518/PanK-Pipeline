#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- Bra_PCF.kmers.txt
my $in1 = $ARGV[1]; ##- sample id
my $in2 = $ARGV[2]; ##- k-mer size
my $out = $ARGV[3]; ##- output 
my $in3 = $ARGV[4]; ##- /mydata/caix/PanK-Pipeline/scripts

my $k            = $in2;

my $filePartsNum = 130;
my $scripts      = $in3;


##-step 1: run jellyfish 
my $jfFile    = $in1.".jf";
my $histoFile = $in1.".histo.gz";
   if(! -e $jfFile){
      &RunJellfish($in1, $k, $jfFile, $histoFile);   
   }

##-step 2: devided query kmer file into small parts
my @fileParts    = ();
   &devidedIntoParts($in0, $filePartsNum, $in1, \@fileParts);
   
##-step 3: count target kmer
  &RunBatchCMD(\@fileParts, $jfFile, $out); 

##-clean
   #system("rm  $jfFile");
   system("pigz $out");

##------------ all subs -----------------------------------------------------------
sub RunJellfish {

    my ($samid, $ksize, $jfFile, $histoFile) = @_;
 
    my $leftRead  = $samid."_1.fq.ft.gz";
    my $rightRead = $samid."_2.fq.ft.gz"; 

    my $tmpFastq = $samid.".tmp.fastq";
    if(! -e $jfFile){
       system("pigz -dc $leftRead  >  $tmpFastq");
       system("pigz -dc $rightRead >> $tmpFastq");
    }
    if(! -e $jfFile){ 
       system("jellyfish  count  -m  $ksize  -t 20 -s 5G  -C  $tmpFastq  -o  $jfFile");        ##- threads: 20 Mem: 5G
    }
    if(! -e $histoFile){
	    #system("jellyfish  histo  -t  20  $jfFile    >  $histoFile");                           ##- threads: 20
    }
    if( -e $tmpFastq){
       system("rm $tmpFastq");
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
