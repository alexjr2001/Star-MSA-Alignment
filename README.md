# STAR-MSA-ALIGNMENT

This repository contains an implementation of the MSA Estrella alignment algorithm in C++. The algorithm uses dynamic programming with the Needleman-Wunsch algorithm to find the best alignment of DNA or RNA sequences.

## Usage

To run the program, make sure you have a C++ compiler installed on your system. Then, simply compile the `Star.cpp` file and execute the resulting binary.

```bash
g++ Star.cpp -o msa_star
./msa_star
```


## Input File
The BRCA1.txt file contains samples of segments of the BRCA1 gene, obtained by LR-PCR amplification. The program reads this file to obtain DNA sequences, which are then used for alignments.

## Contributions
Contributions are welcome! If you have suggestions for improvements to the code or documentation, feel free to open an issue or send a pull request.

## Credits
This code was developed as part of a project for the Computational Molecular Biology course.
