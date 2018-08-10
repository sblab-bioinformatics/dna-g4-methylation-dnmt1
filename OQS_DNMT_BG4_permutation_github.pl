#### Monte Carlo simulation: 
## 1- The peaks provided by set1 argument will be shuffled, 
## 2- Shuffled peaks are overlapped with peaks provided in set2 argument, 
## After all overlap files are generated, The number of peaks overlapping set2 is then counted and written to an output file.
## This output file is then used to calculate Monte Carlo's P-value.
## The number of random peak to choose from the total set, set1, is provided through number option
## The number of peak in the total set, set1, is provided through max option
## shuffD option is the path to the directory to save the shuffled datasets
## overlD option is the path to the directory to save the Ooverlapping datasets

#### Example:
## perl OQS_DNMT_BG4_permutation.pl -set1 open_Na_K_Na_PDS_merged.bed 
##      -set2 optimal_ENCFF549TVW_sorted.bed -number 7491 -max 43506 
##      -shuffD randomised_sets -overlD overlapping_bed_files
 

#!/usr/bin/perl
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my ($set1, $set2, $numN, $maxNum, $shuffDir, $overlDir);

GetOptions ("set1=s" 	=> \$set1,
	    "set2=s"	=> \$set2,
	    "number=i"	=> \$numN,
	    "max=i"	=> \$maxNum,
	    "shuffD=s"	=> \$shuffDir,	   
	    "overlD=s"	=> \$overlDir)
	    or die("Error in command line arguments!\n");

### Sub to create the script for sbatch submission and then submit the job to computer cluster with slurm scheduler
### It returns the job ID assigned to it after the submission to the cluster
### In case of dependencies: depJobID is the job id of the job this job depend on	    
sub submit_job {
        my ($dir, $command, $jobName, $depJobID, $cdFlag) = @_;

	## By default this job does not depend on anything else
	$depJobID //=0;
	## By default it changes the directory to the one provided
	$cdFlag //=1;

	my $script = "\#!/bin/bash \n";
	$script .= "\# set the number of nodes \n";
	$script .= "\#SBATCH --nodes=1 \n";
	$script .= "\#SBATCH --ntasks-per-node=2 \n";
	$script .= "\#SBATCH -o $jobName\.out # File to which STDOUT will be written \n";
	$script .= "\#SBATCH -e $jobName\.err # File to which STDERR will be written \n";


	$script .= "echo \"Job name \" $jobName \n";
	$script .= "echo hostname=\`hostname\` \n"; 
	$script .= "echo \"Starting at \`date\`\" \n";
	$script .= "echo \"Running on hosts\: \$SLURM_NODELIST\" \n";
	$script .= "echo \"Running on \$SLURM_NNODES nodes\.\" \n";
	$script .= "echo \"Running on \$SLURM_NPROCS processors\.\" \n";
	$script .= "echo wd=\$\(pwd\) \n";
	$script .= "echo \$working \n";
	$script .= $command . "\n";

	open WRT , ">", $dir."\/" . $jobName . "_sbatch\.sh" or die "Coulod not open the file to write the batch script!";
	print WRT $script;

	## Change the dir unless the flag is set
	chdir($dir) if($cdFlag);
 
	## Now submit the job to the scheduler on the cluster
	my $cmd = "sbatch ";
	## Does it depend on another job?
	$cmd .= "--dependency=afterok:$depJobID " if($depJobID);
	$cmd .= $dir."\/" . $jobName . "_sbatch\.sh"; 
	my $line = `$cmd`;

	my ($jobID); 
	if($line =~ /(\d+)$/){
		return($1);
	}else{
		print("Error in job submission: ".$line."\n");
	}
	
}


## Constants:
my $iteration = 8000; ## Monte Carlo #iteration, N

my (@index);
## Create the index array to be shuffled
for (my $i=0; $i < $numN; $i++) {
	$index[$i]=1;
}
for (my $i=$numN; $i <=$maxNum; $i++) {
	$index[$i]=0;
}


## Iteration loop: Shuffle and Overlap 
for (my $i=0; $i <= $iteration; $i++) {

	## Shuffle index array. If the index array is 1 the corresponding region will be included in the shuffled set
	my @shuffled = shuffle(@index);

	my $randomF = $shuffDir."\/shuffled_".$i."\.bed"; 
	## Select the lines from set1 bed file
	open SET1, $set1 or die "Could not open set1 bed file! $set1\n";
	open OUT, ">", $randomF or die "Could not open the output file! $randomF \n";

	while(<SET1>){
		
		print OUT $_ if($shuffled[$.]==1);
	}
	close(OUT);
	close(SET1);

	## Overlap the shuffled set1 with set2; using intersect from bedtools
	my $cmd = "intersectBed -a $randomF -b $set2 -wa -u > $overlDir/overlapped_$i.bed ";
	my $overlapJobID= submit_job($overlDir."\/sbatch_logs", $cmd, "overlap_".$i);
	
}
