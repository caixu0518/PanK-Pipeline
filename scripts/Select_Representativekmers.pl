#!/usr/bin/perl -w
use strict;
use Cwd;

my $in0 = $ARGV[0]; ##- merged.kmers.list.Polymorphic_kmers.List
my $in1 = $ARGV[1]; ##- all.cleanfiles.txt
my $in2 = $ARGV[2]; ##- /mydata/caix/PanK-Pipeline/scripts 

my $out = $in0.".representative.list";

my $currentPath        = getcwd();

my $cleanListDir       = "$currentPath/FastaKmerList";
my $fastsiesDir        = "$currentPath/Fastasizes";
my $CoordSortedKmerDir = "$currentPath/CoordSortedKmers";  

#my $scriptsPath        = "/mydata/caix/Brapa_k17_analysis/cleanList/scripts";
my $scriptsPath        = $in2;

my $winsize = 10000;
    
###-----------------------------------------------------------------
my $timestrings = &Times();
   print STDERR "[Current system time: $timestrings] .... Start reading kmers: $in0 .....\n";
my %kmerindex   = ();
   &readKmerList($in0, \%kmerindex);
   $timestrings = &Times();
   print STDERR "[Current system time: $timestrings] .... Finished index kmers: $in0 .....\n\n";

   $timestrings = &Times();
   print STDERR "[Current system time: $timestrings] .... Start reading clean kmer list files: $in1 .....\n";
my @sams = ();
   &createSingleCopyList($in1, \@sams, \%kmerindex);    
   $timestrings = &Times();
   print STDERR "[Current system time: $timestrings] .... Finished reading clean kmer list files: $in1 .....\n\n";

##-read kmers
   %kmerindex   = ();
   $timestrings = &Times();
   print STDERR "[Current system time: $timestrings] .... Start generate representative kmers .....\n";
my %representativeKmers = ();
   &selectrepresentativeKmers(\@sams, \%representativeKmers);
   $timestrings = &Times();
   print STDERR "[Current system time: $timestrings] .... Finished generate representative kmers .....\n\n";

##-generate output
   $timestrings = &Times();
   print STDERR "[Current system time: $timestrings] .... Start generate output file $out .....\n";
   &output($in0, \%representativeKmers, $out);    
   $timestrings = &Times();
   print STDERR "[Current system time: $timestrings] .... Finished generate output file $out .....\n";  

####----- all subs --------------------------------------------------
sub output {

    my ($fileIn, $representativekmers, $fileOut) = @_;

    open (my $FR0000, $fileIn);
    open (my $FOoooo, ">$fileOut");
    while(<$FR0000>){
      chomp;
      my @temp = split(/\t/, $_);
      if(exists $representativekmers ->{$temp[0]}){
         print $FOoooo $_, "\n";
      }
    }
    close ($FOoooo);
    close ($FR0000);

}
   
sub selectrepresentativeKmers {

    my ($sams, $representativeKmers) = @_;

    ##-- assgin sample-kmers into intervals --------------------------------------------
    if(not -d "AssignKmersIntoIntervals"){
       system("mkdir AssignKmersIntoIntervals");
    }
   
    my $batchCMDfile = "AssignKmersIntoIntervals.batch.cmds";  
    if(-e "$batchCMDfile"){
       system("rm -rf  $batchCMDfile");
    } 
    if(-e "$batchCMDfile.completed"){
       system("rm -rf $batchCMDfile.completed");
    }


    open (my $FOCMDs, ">$batchCMDfile");
    my $totalruns = 0;
    my @outfiles  = ();
    for my $eachsam(@{$sams}){
        my $fastasizeFile        = "$fastsiesDir/$eachsam.sizes";
        my $CoordSortedKmerFile = "$CoordSortedKmerDir/$eachsam.coordsorted.list.gz"; 

        if(not -e "./AssignKmersIntoIntervals/$eachsam.sizes.windata.bed.add.Pos.list") {   
           if(-e $fastasizeFile && -e $CoordSortedKmerFile){
              $totalruns += 1;		   
	      system("ln -s $fastasizeFile        .");
	      system("ln -s $CoordSortedKmerFile  .");
	      push(@outfiles, "$eachsam.sizes.windata.bed.add.Pos.list");
              print $FOCMDs "perl  $scriptsPath/AssignKmersIntoIntervals.pl   $eachsam.coordsorted.list.gz  $eachsam.sizes\n"; 
	   }
	   else{
	      die "cannot find ... $fastasizeFile or $CoordSortedKmerFile ... \n";
	   }
	}
    }
    close ($FOCMDs);

    ##- batch run ------
    system("ParaFly -c $batchCMDfile -CPU 30"); 

    ##-- save
    for my $eachfile(@outfiles){
        if(-e $eachfile){
	   system("mv $eachfile AssignKmersIntoIntervals");
	}
    } 

    ##- clean
    for my $eachsam(@{$sams}){
        my $fastasizeFilea       = "$eachsam.sizes";
	my $CoordSortedKmerFilea = "$eachsam.coordsorted.list.gz";
	my $intervalPos          = "$eachsam.sizes.windata.bed.add.Pos.list";
        if(-e $fastasizeFilea){
	   system("rm -rf $fastasizeFilea");
	}
	if(-e $CoordSortedKmerFilea){
	   system("rm -rf $CoordSortedKmerFilea");
	}

        if(-e $intervalPos){
	   system("mv  $intervalPos  AssignKmersIntoIntervals");
	}
    }

    for my $eachsam(@{$sams}){

        my $CoordSortedKmerFile = "$CoordSortedKmerDir/$eachsam.coordsorted.list.gz";

	if(-e $CoordSortedKmerFile){
           
           $timestrings = &Times();
           print STDERR "[Current system time: $timestrings] .... Start process $CoordSortedKmerFile .....\n";

           ##- AA_LongYou.coordsorted.list.gz
	   my %pos2Kmer  = ();
	   my %pos2Value = (); 
	   open (my $FR11, "pigz -dc $CoordSortedKmerFile | ");
           while(<$FR11>){
	     chomp;
	     my @temp                       = split(/\t/, $_);
	     my $tmpfreqValue               = sprintf("%.2f", $temp[3]);
             $pos2Kmer{$temp[0]}{$temp[1]}  = $temp[2];
	     $pos2Value{$temp[0]}{$temp[1]} = $tmpfreqValue; 
	   }
	   close ($FR11);

           ##- mark representative kmers
           my $intervalkmerfile = "./AssignKmersIntoIntervals/$eachsam.sizes.windata.bed.add.Pos.list";  	
           if(-e $intervalkmerfile){
	      $timestrings = &Times();
	      print STDERR "[Current system time: $timestrings] .... Start process $eachsam windata $intervalkmerfile.....\n";   
             
	      open (my $FR002, $intervalkmerfile); 
              while(<$FR002>){
	        chomp;
		my @temp = split(/\t/, $_);
		if(defined $temp[3] && $temp[3] =~ /:/){
		   my @locus = split(/;/, $temp[3]);
		  
		   my %freq2kmer = ();
		   for my $eachtmp(@locus){
		       my @ainfo = split(/:/, $eachtmp);
		       if(exists $pos2Value{$ainfo[0]}{$ainfo[1]}){
		          my $freqvaluetmp = $pos2Value{$ainfo[0]}{$ainfo[1]};
			  my $kmerseqtmp   = $pos2Kmer{$ainfo[0]}{$ainfo[1]};

			  if(not exists $freq2kmer{$freqvaluetmp}){
			     $freq2kmer{$freqvaluetmp} = $kmerseqtmp;
			  }
			  else{
			     $freq2kmer{$freqvaluetmp} .= ";".$kmerseqtmp;
			  }
		       }
		   }  
	           for my $tmpfreq(keys %freq2kmer){
		       my @kmersTmp = split(/;/, $freq2kmer{$tmpfreq});
		       if(scalar(@kmersTmp) == 1){
		          $representativeKmers ->{$kmersTmp[0]} = "Y";
		       }
		       else{
		         for(my $j=1; $j<=$#kmersTmp; $j++){
                             if(exists $representativeKmers ->{$kmersTmp[$j]}){
			        next;
			     }
			     else{
			       $representativeKmers ->{$kmersTmp[$j]} = "Y";
			       last;
			     }
		         }		 
		       }
		   }
		}
	      }
	      close ($FR002);   
	   }
	   else{
            die "cannot find: $intervalkmerfile.\n";
           }		   
	}
	else{
	   die "cannot find: $CoordSortedKmerFile";
	}
    }

}

sub createSingleCopyList {

    my ($cleankmerList, $samArray, $kmerindex) = @_;

    ##- create dir for storing coords list file
    if(not -d "CoordSortedKmers"){
       system("mkdir CoordSortedKmers");
    }
    

    open (my $FR1, $cleankmerList);
    while(<$FR1>){
      chomp;
      my $cleanListfilename = $_;

      my $tmpsamid =  $cleanListfilename;
         $tmpsamid =~ s/\.k\S+\.clean\.list\.gz//;    
    
      push(@{$samArray}, $tmpsamid);  

      my $coordsortedfile = $tmpsamid.".coordsorted.list";   
      if(not -e "./CoordSortedKmers/$coordsortedfile.gz"){
     
	 my $coordsortedfiletmp = $coordsortedfile.".tmp";      
	 open (my $FO0, ">$coordsortedfiletmp");
         open (my $FR11, "pigz -dc  $cleanListDir/$cleanListfilename | ");
         while(<$FR11>){
	   chomp;
	   my @temp = split(/\t/, $_);
	   if($temp[1] == 1){    ##- single copy
	      my @info = split(/:/, $temp[2]);
              if(exists $kmerindex ->{$temp[0]}){
		 my $freqValue = $kmerindex ->{$temp[0]};      
                 print $FO0 join("\t", $info[0], $info[1], $temp[0], $freqValue), "\n";
              }


	   } 
	 }
         close ($FR11); 	 
         close ($FO0);
         
         system("sort -k1,1 -k2,2n $coordsortedfiletmp > $coordsortedfile"); 
	 system("rm -rf $coordsortedfiletmp");
	 system("pigz $coordsortedfile");
	 system("mv $coordsortedfile.gz  CoordSortedKmers");
      
      }
      else{
         $timestrings = &Times();
	 print STDERR "[Current system time: $timestrings] .... Skip, the file $coordsortedfile exists ....\n";
      }    
    }
    close ($FR1);

}

sub readKmerList {

    my ($fileIn, $kemerIndex) = @_;

    open (my $FR0, $fileIn);
    while(<$FR0>){
      chomp;
      my @temp  = split(/\t/, $_);
      my $value = sprintf("%.2f", $temp[2]);
         $kemerIndex ->{$temp[0]} = $value;
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


