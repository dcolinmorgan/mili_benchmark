#!/bin/bash

cat camb_motif_wgbsA549/CG/* > camb_motif_wgbsA549/all_a549.txt
cut -f1,2,3,6 camb_motif_wgbsA549/all_a549.txt | uniq > camb_motif_wgbsA549/uni_all_a549.txt
cut -f4 camb_motif_wgbsA549/uni_all_a549.txt | sort| uniq -c |sort -nr  > camb_motif_wgbsA549/uni_all_a549_count.txt

cat camb_motif_wgbsK562/CG/* > camb_motif_wgbsK562/all_k562.txt
cut -f1,2,3,6 camb_motif_wgbsK562/all_k562.txt | uniq > camb_motif_wgbsK562/uni_all_k562.txt
cut -f4 camb_motif_wgbsK562/uni_all_k562.txt | sort| uniq -c |sort -nr  > camb_motif_wgbsK562/uni_all_k562_count.txt

cat camb_motif_wgbsHepG2/CG/* > camb_motif_wgbsHepG2/all_hepg2.txt
cut -f1,2,3,6 camb_motif_wgbsHepG2/all_hepg2.txt | uniq > camb_motif_wgbsHepG2/uni_all_hepg2.txt
cut -f4 camb_motif_wgbsHepG2/uni_all_hepg2.txt | sort| uniq -c |sort -nr  > camb_motif_wgbsHepG2/uni_all_hepg2_count.txt

cat camb_motif_wgbsHeLa/CG/* > camb_motif_wgbsHeLa/all_hela.txt
cut -f1,2,3,6 camb_motif_wgbsHeLa/all_hela.txt | uniq > camb_motif_wgbsHeLa/uni_all_hela.txt
cut -f4 camb_motif_wgbsHeLa/uni_all_hela.txt | sort| uniq -c |sort -nr  > camb_motif_wgbsHeLa/uni_all_hela_count.txt

cat camb_motif_wgbsSKNSH/CG/* > camb_motif_wgbsSKNSH/all_sknsh.txt
cut -f1,2,3,6 camb_motif_wgbsSKNSH/all_sknsh.txt | uniq > camb_motif_wgbsSKNSH/uni_all_sknsh.txt
cut -f4 camb_motif_wgbsSKNSH/uni_all_sknsh.txt | sort| uniq -c |sort -nr  > camb_motif_wgbsSKNSH/uni_all_sknsh_count.txt

cat camb_motif_wgbsGM12878/CG/* > camb_motif_wgbsGM12878/all_gm12878.txt
cut -f1,2,3,6 camb_motif_wgbsGM12878/all_gm12878.txt | uniq > camb_motif_wgbsGM12878/uni_all_gm12878.txt
cut -f4 camb_motif_wgbsGM12878/uni_all_gm12878.txt | sort| uniq -c |sort -nr  > camb_motif_wgbsGM12878/uni_all_gm12878_count.txt
