#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- kmerCount files
my $in1 = $ARGV[1]; ##- merged.kmers.list

my @kmercountArrays = ();
   &readKmercounts($in0, \@kmercountArrays);

my %lineIndex = ();
   &getPolymorphicKmerLinenum(\%lineIndex, \@kmercountArrays);

my $out = $in1.".Polymorphic_kmers.List";
   &output(\@kmercountArrays, $in1, \%lineIndex, $out);


##----- all subs -------------------------------------------------------------------------------
sub output {

    my ($filearray, $fileIn, $lineIndex, $fileout) = @_;

    my $samnum = scalar(@{$filearray});

    my $lineNum = 0;
    open In00, $fileIn;
    open OUT00, ">$fileout";
    while(<In00>) {
      chomp;
      $lineNum += 1;
      if(exists $lineIndex ->{$lineNum}){
         my $ratio = ($lineIndex ->{$lineNum})/$samnum;
   
         if($ratio > 0 && $ratio < 1){
            print OUT00 $_, "\t", $ratio, "\n";
         }

      }

    }
    close OUT00;
    close In00;

}

sub getPolymorphicKmerLinenum {

    my ($lineIndex, $fileArrays) = @_;

    for my $eachFile(@{$fileArrays}){
 
        ##- start each file
        my $timestrings = &Times();
        print STDERR "[Current system time: $timestrings] .... Start the file: $eachFile .....\n";  

        my $newfiletmp = $eachFile.".tmp";
        system("pigz -dc $eachFile  >  $newfiletmp");
        
        my $count = 0;
        open IN00, $newfiletmp;
        while(<IN00>){
          $count += 1;
          chomp;
          if($_ > 0){
             $lineIndex ->{$count} += 1; 
          }
        }
        close IN00; 
        system("rm -rf $newfiletmp");
    }   

}

sub readKmercounts {

    my ($fileIn, $arrays) = @_;
  
    open IN0, $fileIn;
    while(<IN0>){
      chomp;
      push(@{$arrays}, $_);
    }
    close IN0;

}

sub Times {

    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();

    $year += 1900;
    $mon += 1;
    my $current_time = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year, $mon, $mday, $hour, $min, $sec);
    return($current_time);

}

