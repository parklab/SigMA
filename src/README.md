# Introduction

C implementation of [get_trinuc_norm.R](https://github.com/parklab/SigMA/blob/master/R/get_trinuc_norm.R) based up [kmer-cnt](https://github.com/lh3/kmer-cnt) and [cgranges](https://github.com/lh3/cgranges). The only library dependency is [zlib](https://www.zlib.net).

1. Features
    * Very fast
    * No region sampling for WES or WGS
    * WGS sig3 calculation on the fly

2. How to compile

    ```
    cd src
    make
    ```

3. Usage

    ```
    get_trinuc_counts -b panel_or_WES.bed -o sig3_norm_scale.tsv reference.fa
    ```

4. WES run

    * WES target region size calculated with [bedtk](https://github.com/lh3/bedtk)

        ```
        bedtk sort seqcap_capture_b37.bed | bedtk merge | bedtk sum
        63959731
        ```

    * Get trinucleotide counts (normalized)

        ```
        time get_trinuc_counts -b seqcap_capture_b37.bed -o seqcap_capture_b37.norm $hs37
        
        real	4m12.611s
        user	4m2.172s
        sys	0m7.739s
        ```

5. Output (seqcap_capture_b37.bed)

    ```
    ACA	0.772794
    ACC	1.052235
    ACG	1.660758
    ACT	0.770029
    CCA	1.019742
    CCC	1.335815
    CCG	2.549477
    CCT	1.060978
    GCA	1.045333
    GCC	1.337516
    GCG	2.529033
    GCT	1.098765
    TCA	0.817317
    TCC	1.066214
    TCG	1.801187
    TCT	0.823388
    ACA	0.772794
    ACC	1.052235
    ACG	1.660758
    ACT	0.770029
    CCA	1.019742
    CCC	1.335815
    CCG	2.549477
    CCT	1.060978
    GCA	1.045333
    GCC	1.337516
    GCG	2.529033
    GCT	1.098765
    TCA	0.817317
    TCC	1.066214
    TCG	1.801187
    TCT	0.823388
    ACA	0.772794
    ACC	1.052235
    ACG	1.660758
    ACT	0.770029
    CCA	1.019742
    CCC	1.335815
    CCG	2.549477
    CCT	1.060978
    GCA	1.045333
    GCC	1.337516
    GCG	2.529033
    GCT	1.098765
    TCA	0.817317
    TCC	1.066214
    TCG	1.801187
    TCT	0.823388
    ATA	0.520578
    ATC	0.803866
    ATG	0.774898
    ATT	0.572611
    CTA	0.621507
    CTC	1.005943
    CTG	1.114823
    CTT	0.870075
    GTA	0.757147
    GTC	1.083921
    GTG	0.965086
    GTT	0.759634
    TTA	0.569502
    TTC	0.847367
    TTG	0.741450
    TTT	0.584343
    ATA	0.520578
    ATC	0.803866
    ATG	0.774898
    ATT	0.572611
    CTA	0.621507
    CTC	1.005943
    CTG	1.114823
    CTT	0.870075
    GTA	0.757147
    GTC	1.083921
    GTG	0.965086
    GTT	0.759634
    TTA	0.569502
    TTC	0.847367
    TTG	0.741450
    TTT	0.584343
    ATA	0.520578
    ATC	0.803866
    ATG	0.774898
    ATT	0.572611
    CTA	0.621507
    CTC	1.005943
    CTG	1.114823
    CTT	0.870075
    GTA	0.757147
    GTC	1.083921
    GTG	0.965086
    GTT	0.759634
    TTA	0.569502
    TTC	0.847367
    TTG	0.741450
    TTT	0.584343
    ```
