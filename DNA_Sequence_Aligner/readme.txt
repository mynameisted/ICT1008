
.==========================================================================================.
|																						   |
|    ___  _  _   _     ___                                    _   _ _                      |
|   |   \| \| | /_\   / __| ___ __ _ _  _ ___ _ _  __ ___    /_\ | (_)__ _ _ _  ___ _ _    |
|   | |) | .` |/ _ \  \__ \/ -_) _` | || / -_) ' \/ _/ -_)  / _ \| | / _` | ' \/ -_) '_|   |
|   |___/|_|\_/_/ \_\ |___/\___\__, |\_,_\___|_||_\__\___| /_/ \_\_|_\__, |_||_\___|_|     |
|                                 |_|                                |___/                 |
|																						   |
.==========================================================================================.


Every living organism’s cells have their own Deoxyribonucleic Acid (DNA) that carries genetic information. Using algorithms for pattern recognition and matching, biologists can learn more about new DNA, RNA and Protein sequences by comparing them to existing known sequence,learning more about a new sequence and its characteristics.

The basic concept involves pattern recognition or pattern matching for DNA, RNA and protein sequences. To put things into perspective, when a new gene is discovered, the biologists use information available from previously known genes to find some characteristics of the new gene. There are different algorithms which can be utilized for sequence alignment and our team have chosen the Needleman- Wunsch algorithm in this program.

The Needleman–Wunsch algorithm is an algorithm developed by Saul B. Needleman and Christian D. Wunsch and published in 1970, commonly used in bioinformatics to globally align protein or nucleotide sequences with the highest quality of alignment. The algorithm uses dynamic programming to compare the sequences, essentially divides a large problem (e.g. the full sequence) into a series of smaller problems (one pair of amino acid) and uses the solutions to the smaller problems to reconstruct a solution to the larger problem.

Assessing the relationships between two sequences can be done by counting the number of identical and similar amino acids. The number of identical and similar amino acids may then be compared to the total number of amino acids in the protein, giving the percentage of sequence identity and similarity. 

Gap penalty is a scoring system used to align a small portion of genetic code. The biological process of protein synthesis (transcription and translation or DNA replication) can produce errors resulting in mutations in the final nucleic acid sequence. In order to make more accurate decisions in aligning reads, mutations are annotated as gaps in the sequence. They are represented as contiguous dashes on a protein/DNA sequence alignment. The scoring that occurs in Gap Penalty allows for the optimisation of sequence alignment in order to obtain the best alignment possible based on the information available.

The program video is available at http://www.tinyurl.com/seqaligner

--==--==--==--==--==--==--

This program is developed to run on Python 2.7.14, using the built-in tkinter module for the GUI.
These dependencies can be downloaded at the following pages: 
[Python 2.7.14 Release] https://www.python.org/downloads/release/python-2714/

--==--==--==--==--==--==--

[Execution Instructions]
To compile and run the program via the Command Prompt/Terminal, please run main.py prefixed with your python intepreter.

--==--==--==--==--==--==--

The DNA Sequence Aligner is developed by a team of Undergraduates from Singapore Institute of Technology's Degree Programme as a group project in ICT1008-Data Structures and Algorithms AY/18.

Information and Communications Technology (Software Engineering), BEng (Hons)

- Lai Wei Fang Ted 		1802183
- Ang Zi Yan			1802851
- Chua Woi Zhao			1801961
- Yap Yong Sheng		1801233

--==--==--==--==--==--==--
Latest version:
[5 Apr 2019] V1 Release - Revamped GUI.
