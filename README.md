# DNA-rep
Code for DNA replication projects
## References
Das S et al. 2015. [Replication timing is regulated by the number of MCMs loaded at origins.](http://genome.cshlp.org/content/25/12/1886.abstract)

## BioNano
Trial filtering and segmenting protocol for BioNano data

Example: Sync_HeLA_1708

```bash
./direction_bed_seg_new_format.py -b input.bnx -o output.bed
```
Input:
* -b : bnx file WITHOUT HEADER (header lines start with #).  
 Assumes green labeling in Channel 1 and red labeling in Channel 2.  
 *\# Nickase Recognition Site 1:   gctcttc;green_01*  
 *\# Nickase Recognition Site 2:   cacgag*  
 
  If opposite use [direction_bed_multi_segment.py](BioNano/direction_bed_multie_segement.py).  
  Further information on bnx file format here: [30038-Rev-B-BNX-v1.2-File-Format-Specification-Sheet.pdf](https://bionanogenomics.com/wp-content/uploads/2017/03/30038-Rev-B-BNX-v1.2-File-Format-Specification-Sheet.pdf)
   




Output: 

* .bedGraphs for unfiltered red label tracks
* .bedGraphs for filtered red label tracks
* bedfile







