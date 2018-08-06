## Monte Carlo Simulation

Monte Carlo simulation was used to calculate the significance of overlap between BG4 peaks and high confidence DNMT1 peaks, defined by Irreproducible Discovery Rate (IDR) in ENCODE’s ChIP pipeline. 

We first counted how many BG4 peaks overlapped with OQSs in open chromatin (defined as all OQSs seen in potassium and/or PDS conditions (749,339 sequences), which overlap at least one DHS region (43,506 sequences)). We then randomly selected the same number of OQSs from all OQSs in open chromatin and counted how many overlapped with at least one high confidence DNMT1 peak. 

## Monte Carlo P-value calculation

The Monte Carlo P-value was calculated as (N+1)/(M+1), where M is the number of iterations and N is the number of times the same or more overlaps were observed between randomised OQSs and high confidence DNMT1 peaks (compared to the number of overlaps observed between BG4 peaks and high confidence DNMT1 peaks). Randomisation was repeated for 8000 times and on average the number of overlaps between the shuffled OQSs and DNMT1 were two-fold less than those observed between BG4 and DNMT1 peaks.
