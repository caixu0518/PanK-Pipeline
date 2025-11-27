#!/usr/bin/perl -w
use strict;
use threads;
use Cwd;
use Getopt::Long;
#use File::Which;

#--Usage-----------------------------------------

my $usage=<<USAGE;
  
   The current script $0 is used to aggregate the genotyping results from each resequenced accession into a single VCF file.
   
   -Sam         [Required]  The file contains all Resequenced sample ID. Resequenced sample ID (must start with a letter, recommended length <8 characters).
   -KmerGTDir   [Required]  Directory for storing the genotyping results of each resequenced accession
   -KmerList    [Required]  The file was used as the input for genotyping k-mers in each resequencing accession. 
                            Be sure to use the same file, as all genotyping result files have the same number of lines, which (starting from 1) serve as positions in the VCF. 
   -threads     [Optional]  The number of threads used in the current script. default: 60

   Bug reports: Xu Cai
                caixu\@caas.cn
		18,8,2025

USAGE


my $currentPath = getcwd();
if (@ARGV == 0){die $usage}

my $in0     = $ARGV[0]; ##- all samples list 
my $in1     = $ARGV[1]; ##- Directory for storing the genotyping results of each resequenced accession
my $in2     = $ARGV[2]; ##- Brassica_rapa.representatvie.OneDepthKmers.list
my $threads = 60;

GetOptions(
    "Sam:s"                      =>\$in0,	
    "KmerGTDir:s"              =>\$in1,
    "KmerList:s"               =>\$in2,
    "threads:i"                =>\$threads,
);

if(not defined $in0 || not defined $in1 || not defined $in2){
   die $usage;
}

my $dataDir = $in1;

my $bgzip   = "bgzip";
my $tabix   = "tabix";
   &checkdependence();

my @samArray  = ();
my @fileArray = ();
   &checkexsist($dataDir, $in0, \@samArray, \@fileArray);

##-- generate header and position file
my %posindex = ();
my $header   = "gt.header"; 
my $pos      = "gt.pos";
    &generateHeader($in2, $header, $pos, \%posindex);

my %samPart  = ();
my %filePart = ();
   &devidedIntoParts(\@samArray, \@fileArray, $threads, \%samPart, \%filePart);


##-- multi-threads 
my $revisedthreads = scalar(keys %samPart);   

my @outgts = ();
   my @thr = ();
   for(my $i=1; $i<=$revisedthreads; $i++){
       my $gtfile = "gt.part".$i.".List";
       system("rm -rf $gtfile"), if(-e $gtfile); 
       #system("touch $gtfile"); 

       push(@outgts, $gtfile);
       my $samcount    = scalar(split(";", $samPart{$i}));
       my $timestrings = &Times();
       print STDERR  "[Current system time: $timestrings] Thread $i: $gtfile, Samples: $samcount --- ", $samPart{$i}, "--- \n";
       $thr[$i] = threads->create(\&generatePartGT, $samPart{$i}, $filePart{$i}, $gtfile);
   }
   for(my $i=1; $i<=$revisedthreads; $i++){
       $thr[$i]->join;
   }

##- merge gt file
   my $infos = join("\t", @outgts);
   my $cmds  = "paste  $pos  $infos > All.gt.pos.List";
   system("$cmds");

##- clean   
   for my $each(@outgts){
       system("rm $each");
   } 
   system("cat  $header All.gt.pos.List  > kmer.gt.vcf"); 
   system("rm $pos $header All.gt.pos.List"); 

my $vcflog = "kmer.gt.vcf.log";
   &calculatemaf("kmer.gt.vcf", $vcflog);     

   system("$bgzip kmer.gt.vcf");
   system("$tabix -p vcf kmer.gt.vcf.gz");
   #system("VCF2Dis_multi -InPut  kmer.gt.vcf   -Threads 80   -OutPut  kmer.gt.vcf.mat");

###---- all subs -------------------------------------------------
sub calculatemaf {
   
    my ($vcfIn, $vcfmaflog) = @_;

    open (my $FR000, $vcfIn);
    open (my $FO, ">$vcfmaflog");
    print $FO join("\t", "#CHROM", "POS", "Refcount/Total", "MAF"), "\n";
    while(<$FR000>){
      chomp;
      next, if(/^#/);
      my @temp     = split(/\t/, $_);
      my $count    = 0;
      my $refCount = 0;
      for(my $m=9; $m<=$#temp; $m++){
          $count += 1;
	  if($temp[$m] eq '0/0'){
	     $refCount += 1;
	  }
      }
      my $maf = $refCount/$count;
      if($refCount > $count/2){
         $maf = ($count - $refCount)/$count;
      }
      print $FO join("\t", $temp[0], $temp[1], "$refCount\/$count", $maf), "\n";
    }
    close ($FO);
    close ($FR000);

}

sub devidedIntoParts {

    my ($samArray, $fileArray, $threads, $samPart, $filePart) = @_;

    my @tmpsam   = @{$samArray};
    my $totalsam = scalar(@tmpsam);
    my $maxindex = $totalsam -1;
    my $eachPart = int($totalsam/$threads) + 1;

    my $m = 1; 
    my %finishedindex = (); 
    for($m=1; $m*$eachPart <= $maxindex; $m++){
        
	my ($start, $end) = (($m-1)*$eachPart, $m*$eachPart);   ##- [ )   
        
	#if(($m-1)*$eachPart < $maxindex && $m*$eachPart > $maxindex){
	#   ($start, $end) = (($m-1)*$eachPart, $#tmpsam + 1);
	#}

        for(my $i=$start; $i<$end; $i++){
            $finishedindex{$i} = "Y";		
	    my $samid  = $samArray ->[$i];
	    my $fileid = $fileArray ->[$i];

	    if(not exists $samPart ->{$m}){
	       $samPart ->{$m} = $samid;
	    }
	    else{
	       $samPart ->{$m} .= ";".$samid;
	    }
 
            if(not exists $filePart ->{$m}){
	       $filePart ->{$m} = $fileid;
	    }
	    else{
	       $filePart ->{$m} .= ";".$fileid;
	    }
	}
    }
    
    ##- my last bin
    my @lastbin = ();
    for(my $p=0; $p<=$maxindex; $p++){
        if(not exists $finishedindex{$p}){
	   push(@lastbin, $p);
	}
    }
   
    if(scalar(@lastbin) > 0){

       for(my $j=0; $j<=$#lastbin; $j++){	
           my $samtmp = $samArray ->[$lastbin[$j]];   
           if(not exists $samPart ->{$m}){
              $samPart ->{$m} = $samtmp;
           }
	   else{
	      $samPart ->{$m} .= ";".$samtmp;
	   }

           my $filetmp = $fileArray ->[$lastbin[$j]];
	   if(not exists $filePart ->{$m}){
	      $filePart ->{$m} = $filetmp;
	   }
	   else{
	      $filePart ->{$m} .= ";".$filetmp;
	   }
       }
    }   

}

sub generatePartGT {
 
    my ($samList, $fileList, $gtFile) = @_;

    my @samList  = split(";", $samList);
    my @fileList = split(";", $fileList);

    my @tmpouts = ();
    for(my $m=0; $m<=$#samList; $m++){
        my $samid  = $samList[$m];
	my $fileid = $fileList[$m];

	my $tmpout = $samid."tmpxxx";
	push(@tmpouts, $tmpout);
        open (my $FR2, "pigz -dc $fileid | ");
        open (my $FO2, ">$tmpout");
	print $FO2 $samid, "\n";
        while(<$FR2>){
	  chomp;
	  ###----  value > 0 will genotyped as C 0/0
	  if($_ > 0){
	     print $FO2 '0/0', "\n";    ##----   it depends_
	  }
	  else{
	     print $FO2 '1/1', "\n";
	  }
	}
	close ($FO2);
	close ($FR2);

    }

    ##- output  
    my $infotmp = join("\t", @tmpouts);
    my $tmpcmd  = "paste $infotmp > $gtFile";
    system("$tmpcmd");
    for my $each(@tmpouts){
	   system("rm -rf $each");
    }

}

sub generateHeader {

    my ($fileIn, $headerout, $posout, $posindex) = @_;

    my $count = 0;
    open (my $FR1, $fileIn);
    open (my $FO1, ">$posout");
    print $FO1 join("\t", "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"), "\n";
    while(<$FR1>){
      $count += 1; 
      $posindex ->{$count} = 0;
      print $FO1 join("\t", "A1", $count, ".", "C", "G", ".", ".", ".", "GT"), "\n";
    }
    close ($FO1);
    close ($FR1);

    ##- output vcf header
    open (my $FO0, ">$headerout");
    print $FO0 "##fileformat=VCFv4.2", "\n";
    print $FO0 "##contig=<ID=A1,length=$count>", "\n";
    close ($FO0);

}
   
sub checkexsist {

    my ($dataDir, $samlist, $samArray, $fileArray) = @_;

    open (my $FR0, $samlist);
    while(<$FR0>){
      chomp;
      my $tmpsamid = $_;
      my $filename = "$dataDir/$tmpsamid.gz";

      if(-e $filename){
         push(@{$fileArray}, $filename);
	 push(@{$samArray}, $tmpsamid);
      }
      else{
         die "cannot find: --- $filename --- \n";
      }
    }
    close ($FR0);

}

sub checkdependence {
 
    ##- tabix 
    my $timestrings = &Times();
    print STDERR "[Current system time: $timestrings] Check if the system has installed the tabix software ... \n";
    my $tabixinfo = `tabix  2>&1`;
    $timestrings  = &Times();
    if($tabixinfo =~ /Version/){
       print STDERR "[Current system time: $timestrings] The tabix software has been detected.\n";
    }
    else{
       die "[Current system time: $timestrings] Please install tabix first.";
    }

    ##- bgzip 
    my $bgzipinfo = `bgzip -h  2>&1`;
    $timestrings  =  &Times();
    if($bgzipinfo =~ /Options/){
       print STDERR "[Current system time: $timestrings] The bgzip software has been detected.\n";
    }
    else{
       die "[Current system time: $timestrings] Please install bgzip first.";
    }
     
}

sub Times {

    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
  
    $year += 1900;
    $mon  += 1;
    my $current_time = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year, $mon, $mday, $hour, $min, $sec);
    return($current_time);

}
