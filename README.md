# phage-tRNAs
Calculating and comparing the GC content, effective number of codons (EnC), relative synonymous codon usage (RSCU), and tRNA adaptation index (tAI) for multiple phages infecting multiple hosts

This repository contains 4 scripts for comparing GC content, EnC, RSCU, and tAI between phages and their host bacteria. For all scripts, the phage(s) are compared to the host(s), with an additional option to compare the ∆ values of the indices (phage - host) to compare phages regardless of the genome composition of their hosts. Comparisons can be done both within a host species and between host species (globally).

The "FORMAL_GC" script calculates the GC content of phages and their hosts and compares each phage to its designated host. These values are then used to calculate ∆GC for the phages. ∆GC is used to compare phages based on tRNA group (High: > 5 tRNAs, Low: 1-5 tRNAs, None: 0 tRNAs) and lifestyle (virulent or temperate). Violin plots compare GC of phages to their hosts, heatmaps compare ∆GC across phage groups.

The "FORMAL_EnC" script calculates the effective number of codons in phages and their hosts. These values are then used to calculate ∆EnC for comparisons between phage groups (as above). The EnC/GC3 correlation is then calculated and plotted using Nc plots. Nc plots are provided for both phage vs. host, comparisons between phage tRNA groups, and comparisons between phage lifestyle groups.

The "FORMAL_RSCU" script calculates the relative synonymous codon usage of phages and their respective hosts. RSCU is calculated and compared on a per-codon basis so specific codon usage patterns can be compared. The ∆RSCU is then calculated per phage and can be used to compare phages by tRNA group or lifestyle group. Heatmaps are provided for both phage to host comparisons as well as phage group comparisons. A PERMANOVA can also be run on both lifestyle and tRNA group to determine if those affect RSCU patterns within a host. Global PERMANOVA can also be performed if multiple hosts are being investigated.

The "FORMAL_tAI" script calculates the tRNA adapation index of phages and their hosts using the host tRNA pool and compares the tAI of the phages to the tAI of the host(s). The ∆tAI is then calculated and can be used to compare phages based on tRNA group or lifestyle. Violin plots compare tAI between host and phages, as well as violin plots to compare ∆tAI of phages based on tRNA group or lifestyle.

Options for statistical analyses are provided in each script.
