The data was simuluated by translating "A Modest Proposal" by Jonathon Swift (project Gutenburg ID 1080) into nucleotides using a tetranucleotide substitution for ASCII characters. Degeneracies were inserted using letters worth more than 5 points in English Scrabble: Q, Z, J, X, and K. From this long nucelotide sequences, 50 nucleotide amplicon regions were extracted and used to build a custom database of sequences. Regions were selected to the combined reconstruction of the three sequences had properties described in Table 1, including sequences which were missing 1 or 2 regions, sequences which shared a single or multiple regions and sequences which were unique but differed by only a few nucleotides in certain regions. 

ASV tables for each region were built by randomly generating abundance weights for each sequence according to a power law distribution. For each individual, 25 sequences were selected randomly to be "present" with equal weight across all sequences; ref sequences were then drawn in a power law distribution where the weighting of a given sequence was the known weight if it was "present" or 1e-20 if it was not. In cases where degenerate sequences were present, a degenerate was chosen at random between the two options. There were 10,000 sequences per sample drawn at random. Error was not allowed in these sequences to minizime mismatching. (The goal of this was to see if we could reconstruct the table.) Representative Sequence and count tables were assembled for all samples, and any ASV with fewer than 10 counts across the test samples was removed to mimic the behavior of common denoising algorithms. 

Integration testing is performed using the following workflow.

1. 

For integration testing, sequence counts were reconstructed without accounting for degeneracy (`--no-count-degenerates` in `reconstruct-database-from-regions`). Additionally, the `--min-freq` was set to 1e-4 in reconstruction, corresponding to 

Table 1. Properties of sequences. 

| Sequence ID | Region 1 | Regions 2 | Region 3 | Weight | 
| --- | --- | --- | --- | --- |
| One Region | | | | | 
| `001` | -- | -- | `001` | 1e-4 |
| `002` | -- | -- | `002` | 1e-4 |
| `003` | -- | -- | `003` | 1e-4 |
| `004` | -- | -- | `004` | 1e-4 |
| `005` | -- | -- | `005` | 1e-4 |
| Two regions, both shared | | | | | 
| `006` | -- | `006\007` | `006\007` | 1e-4 |
| `007` | -- | `006/007` | `006/007` | 1e-4 |
| `008` | -- | `008/009` | `008/009` | 1e-3 |
| `009` | -- | `008/009` | `008/009` | 1e-4 |
| `010` | -- | `010/011` | `010/011` | 1e-3 |
| `011` | -- | `010/011` | `010/011` | 1e-3 |
| Two regions, one shared | | | | | 
| `012` | -- | `012/013` | `012` | 1e-4 |
| `013` | -- | `012/013` | `013` | 1e-4 |
| `014` | -- | `014/015` | `014` | 1e-4 |
| `015` | -- | `014/015` | `015` | 1e-4 |
| `016` | -- | `016/017` | `016` | 1e-4 |
| `017` | -- | `016/017` | `017` | 1e-4 |
| `018` | -- | `018/019` | `018` | 1e-4 |
| `019` | -- | `018/019` | `019` | 1e-4 |
| Two regions, both unique | | | | | 
| `020` | -- | `020` | `020` | 1e-4 |
| `021` | -- | `021` | `021` | 1e-4 |
| `022` | -- | `022` | `022` | 1e-3 |
| `023` | -- | `023` | `023` | 1e-4 |
| `024` | -- | `024` | `024` | 1e-2 |
| `025` | -- | `025` | `025` | 1e-4 |
| `026` | -- | `026` | `026` | 1e-4 |
| `027` | -- | `027` | `027` | 1e-4 |
| `028` | -- | `028` | `028` | 1e-4 |
| `029` | -- | `029` | `029` | 1e-4 |
| Three regions, all shared | | | | | 
| `030` | `030\031` | `030\031` | `030\031` | 1e-3 |
| `031` | `030/031` | `030\031` | `030\031` | 1e-4 |
| `032` | `032\033` | `032\033` | `032\033` | 1e-4 |
| `033` | `032\033` | `032\033` | `032\033` | 1e-4 |
| Three regions, two shared | | | | | 
| `034` | `034/035` | `034/035` | `034` | 1e-4 |
| `035` | `034/037` | `034/035` | `035` | 1e-4 |
| `036` | `036/037` | `036/037` | `036` | 1e-3 |
| `037` | `036/037` | `036/037` | `037` | 1e-3 |
| `038` | `038/039` | `038/039` | `038` | 1e-4 |
| `039` | `038/039` | `038/039` | `039` | 1e-4 |
| Three regions, one shared | | | | |
| `040` | `040/041` | `040` | `040` | 1e-4 |
| `041` | `040/041` | `041` | `041` | 1e-4 |
| `042` | `042/043` | `042` | `042` | 1e-4 |
| `043` | `042/043` | `043` | `043` | 1e-4 |
| `044` | `044/045` | `044` | `044` | 1e-4 |
| `045` | `044/045` | `045` | `045` | 1e-4 |
| `046` | `046/047` | `046` | `046` | 1e-4 |
| `047` | `046/047` | `047` | `047` | 1e-4 |
| `048` | `048/049` | `048` | `048` | 1e-4 |
| `049` | `048/049` | `049` | `049` | 1e-4 |
| Differ by 1 nucelotide in all three regions from the next sequence | | | | | 
| `050` | `050` | `050` | `050` | 1e-4 |
| `051` | `051` | `051` | `051` | 1e-4 |
| `052` | `052` | `052` | `052` | 1e-4 |
| `053` | `053` | `053` | `053` | 1e-4 |
| Differ by 2 nucelotides in one region from the next sequence | | | | | 
| `054` | `054` | `054` | `054` | 1e-4 |
| `055` | `055` | `055` | `055` | 1e-3 |
| `056` | `056` | `056` | `056` | 1e-4 |
| `057` | `057` | `057` | `057` | 1e-4 |
| `058` | `058` | `058` | `058` | 1e-4 |
| `059` | `059` | `059` | `059` | 1e-4 |
| All regions are unique with no build in redundency | | | | | 
| `060` | `060` | `060` | `060` | 1e-4 |
| `061` | `061` | `061` | `061` | 1e-4 |
| `062` | `062` | `062` | `062` | 1e-4 |
| `063` | `063` | `063` | `063` | 1e-4 |
| `064` | `064` | `064` | `064` | 1e-4 |
| `065` | `065` | `065` | `065` | 1e-4 |
| `066` | `066` | `066` | `066` | 1e-4 |
| `067` | `067` | `067` | `067` | 1e-3 |
| `068` | `068` | `068` | `068` | 1e-4 |
| `069` | `069` | `069` | `069` | 1e-4 |
| `070` | `070` | `070` | `070` | 1e-4 |
| `071` | `071` | `071` | `071` | 1e-4 |
| `072` | `072` | `072` | `072` | 1e-4 |
| `073` | `073` | `073` | `073` | 1e-3 |
| `074` | `074` | `074` | `074` | 1e-4 |
| `075` | `075` | `075` | `075` | 1e-4 |
| `076` | `076` | `076` | `076` | 1e-4 |
| `077` | `077` | `077` | `077` | 1e-4 |
| `078` | `078` | `078` | `078` | 1e-4 |
| `079` | `079` | `079` | `079` | 1e-4 |
| `080` | `080` | `080` | `080` | 1e-4 |
| `081` | `081` | `081` | `081` | 1e-4 |
| `082` | `082` | `082` | `082` | 1e-4 |
| `083` | `083` | `083` | `083` | 1e-4 |
| `084` | `084` | `084` | `084` | 1e-4 |
| `085` | `085` | `085` | `085` | 1e-4 |
| `086` | `086` | `086` | `086` | 1e-2 |
| `087` | `087` | `087` | `087` | 1e-4 |
| `088` | `088` | `088` | `088` | 1e-3 |
| `089` | `089` | `089` | `089` | 1e-3 |
| `090` | `090` | `090` | `090` | 1e-4 |
| `091` | `091` | `091` | `091` | 1e-4 |
| `092` | `092` | `092` | `092` | 1e-4 |
| `093` | `093` | `093` | `093` | 1e-4 |
| `094` | `094` | `094` | `094` | 1e-4 |
| `095` | `095` | `095` | `095` | 1e-3 |
| `096` | `096` | `096` | `096` | 1e-4 |
| `097` | `097` | `097` | `097` | 1e-4 |
| `098` | `098` | `098` | `098` | 1e-4 |
| `099` | `099` | `099` | `099` | 1e-4 |
| `100` | `100` | `100` | `100` | 1e-4 |





