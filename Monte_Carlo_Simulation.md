## Monte Carlo Simulation

Monte Carlo simulation was used to calculate the significance of overlap between BG4 peaks and high confidence DNMT1 peaks, defined by Irreproducible Discovery Rate (IDR) in ENCODE’s ChIP pipeline. 

We first counted how many BG4 peaks overlapped with OQSs in open chromatin (defined as all OQSs seen in potassium and/or PDS conditions (749,339 sequences), which overlap at least one DHS region (43,506 sequences)). We then randomly selected the same number of OQSs from all OQSs in open chromatin and counted how many overlapped with at least one high confidence DNMT1 peak. 

## Monte Carlo P-value Calculation

The Monte Carlo P-value was calculated as (N+1)/(M+1), where M is the number of iterations and N is the number of times the same or more overlaps were observed between randomised OQSs and high confidence DNMT1 peaks (compared to the number of overlaps observed between BG4 peaks and high confidence DNMT1 peaks). Randomisation was repeated for 8000 times and on average the number of overlaps between the shuffled OQSs and DNMT1 were two-fold less than those observed between BG4 and DNMT1 peaks.

## Code

Here is a more detailed explanation and the script used to calculate Monte Carlo P-value.

### Perl script 
[This](OQS_DNMT_BG4_permutation_github.pl) is the perl script used to shuffle and overlap the sets on a computer cluster with SLURM scheduler.

Let us first clarify the question we are trying to answer with Monte Carlo Simulation. We are asking whether the overlap obswrved between the BG4 peaks and DNMT1 peaks in K562 cell line is statistically significant. to do this we need to randomly senect a number of regiongs from the total G4 structures in open chromatin to which DNMT1 may potentially bind. So set1 is the G4-seq detected sequences which coincide with an open chromatin region provided by DHS in K562. We call this set open_Na_K_G4.bed. Total number of regions in this file is provided through the option - max. Set2 should be DNMT1 peaks in K562 cell line, DNMT1_K562.bed. 

The only option remained to run the above perl script is how many random regions to select from the whole set1, open_Na_K_G4.bed. This is basically the number of BG4 peaks which overlap any G4 sequences observed by G4-seq in open chromatin. That is 7491 peaks.

We just need to provide two directories to the script, one for the shuffled or randomly selected regions, randomised_sets, and another directory for the overlapping bed files, overlapping_bed_files.

So here is an example of how to run above perl script:
```
perl OQS_DNMT_BG4_permutation.pl -set1 open_Na_K_G4.bed -set2 DNMT1_K562.bed -number 7491 -max 43506 -shuffD randomised_sets -overlD overlapping_bed_files
```

Now that we have the overlapping bed files we shall count the number of overlaps observed in each randomised set:

```
wc -l overlapping_bed_files/* > number_of_overlaps.txt
```

To calculate Monte Carlo's P-value, we just need to find out how many of the randomised sets reach the same or more number of overlaps observed as it was observed between the BG4 and DNMT1 peaks files in K562 cell line, which was 400 overlaps:
```
awk '$1>399' number_of_overlaps.txt 
```

As out of 8000 iteration, we never observed equal or more overlaps than 400 then P-value = 1/8001 = 0.000125
