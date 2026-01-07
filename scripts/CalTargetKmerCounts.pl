#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- Bra_PCF.kmers.txt
my $in1 = $ARGV[1]; ##- mer_counts.jf
my $out = $ARGV[2]; ##- output file

my $batchSize = 50000;

#my $out = $in1.".targetKmerCount.txt";
   &calculateKcounts($in0, $in1, $batchSize, $out);


sub calculateKcounts {
  
    my ($kmerFile, $jfFile, $batchSize, $output) = @_;

    ##- file line num
    my $lineCount = `wc -l  $kmerFile`;
       chomp($lineCount);
       $lineCount =~ s/\s+\S+$//;

    open OUTKmerC, ">$output";
    my $r = 0;
    my $count = 0;
    my @kmers = ();
    open IN001, $kmerFile;
    while(<IN001>){
      $count += 1;
      my @temp = split(/\t/, $_);
      push(@kmers, $temp[0]);    

      if(int($count/$batchSize) == $count/$batchSize || $count == $lineCount){
         my @kmerCounts = ();
         my $kmerInfo = join("  ", @kmers);
            &kmerQuery($jfFile, $kmerInfo, \@kmerCounts); 
            print OUTKmerC join("\n", @kmerCounts), "\n";
            @kmerCounts = ();
            @kmers = ();
             
            if($r % 5 == 0){
                my $finishedRatio = sprintf "%.2f", ($count/$lineCount)*100;
                my $currentTime = &Times();
                print "[$currentTime]--Process $kmerFile --:  $finishedRatio %\n";
            }
            $r += 1;
      } 
    }       
    close IN001;
    close OUTKmerC;

}

sub kmerQuery {

    my ($jfFile, $kmerArry, $kcountArray) = @_;

    my $a = `jellyfish query  $jfFile  $kmerArry`;
    chomp($a);
    my @ainfo = split(/\n/, $a);   
    for(@ainfo){
        my @each = split(/\s+/, $_);
        push(@{$kcountArray}, $each[1]);
    }  
    
}

sub Times {

    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();

    $year += 1900;
    $mon += 1;
    my $current_time = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year, $mon, $mday, $hour, $min, $sec);
    return($current_time);

}
