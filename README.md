
***

Sequence Alignment
==================

Intro
-----
*SequenceAlignment* app contains solutions for two exercises from lessons of *Computional Biology*:

1. Optimal aligment of 2 sequences 
	* with penalty function 
	* without penalty function; but with linear memory complexity
2. Aligment of multiple sequences
	* consensus word, multialignment profile
	* alignment of 2 multialignments
	* progressive multi alignment of multiple sequences

Application has been written in F#. Source code can be found on *GitHub*: [https://github.com/theimowski/SequenceAlignment](https://github.com/theimowski/SequenceAlignment)

Application manual
------------------
This is a console app, to invoke it from command line:

	SequenceAlignment.exe <command> [-v|--verbose]

Available commands: 

* Gotoh

	Run algorithm for aligning 2 sequences with penalty function ``f(x) = -x -1``, known as *Gotoh*.

* NeedlemanWunsch

	Run Needleman-Wunsch algorithm for aligning 2 sequences without penalty function.

* Hirschberg

	Run Hirschberg algorithm for aligning 2 sequences without penalty function.

* profile

	Count profile of a multialignment.

* cons 

	Count consensus word for a multialignment.

* malign 

	Align two multialignments into one.

* UPGMA

	Run progressive multi align algorithm for multiple sequences.
	Unweighted Pair Grouping Method with Arithmetic Mean has been used for choosing two clusters to align in each step of the algorithm.

Flags:

* [*-v|--verbose*]
	
	Writes steps of current algorithm to the standard output.

Input data format
-----------------

The application operates only on the following alphabet of Nucletoids : [A;C;G;T]

Program requires similarity matrix on the standard input. Given:

	 1;
	 0; 1;
	 0; 0; 1;
	 0; 0; 0; 1;
	 0; 0; 0; 0; 1;

The similarity matrix is parsed to the following :

   | A | C | G | T | -
---|---|---|---|---|---
 A | 1 | 0 | 0 | 0 | 0
 C | 0 | 1 | 0 | 0 | 0
 G | 0 | 0 | 1 | 0 | 0
 T | 0 | 0 | 0 | 1 | 0
 - | 0 | 0 | 0 | 0 | 1

For *Gotoh*, there's no need to lookup similarity of ``-`` (a space) with other Nucleotide, because the algorithm uses the penalty function. Therefore for *Gotoh* command the **last line of the matrix is ignored**.

Likewise, for *profile* command, the similarity matrix is **completely ignored**.

Despite the fact that the lines are being ignored, they are required for the program to run.

--- 

Next, program expects an empty line on the standard input.
List of sequences / multialignment comes afterwards.

For *Gotoh*, *NeedlemanWunsch* and *Hirschberg* commands program reads two seqences:

	ACTATGGGGTTCC
	ACCCGGGTTC

For *profile* and *cons* commands program reads multialignment:

	A-TCCC
	AGTC-C
	-GTCCC
	AGTCC-

For *malign* command program reads two multialignments, separated by new line: 

	A-TCCC
	AGTC-C
	-GTCCC
	AGTCC-

	GTCC-
	GTCCC
	-TAC-

For *UPGMA* command program reads sequences until EOF:

	ACTATGGGGTTCC
	ACCCGGGTTC
	ACCCGGGTTGTC
	ACCGGCTTAG
	CCCCGTTAGCC

---

Example of correct program execution:

	.\SequenceAlignment.exe Hirschberg

Standard input:

	 2
	 0; 2
	 0; 0; 2
	 0; 0; 0; 2
	-1;-1;-1;-1; 2

	ACTGATTAA
	ATAA

Standard output:

	ACTGATTAA
	A-----TAA

Implementation details
----------------------

* Gotoh

	The algorithm is based on Dynamic Programming. It has been implemented with the help of 4 matrices:

	* Matrix A assumes that the last element from first sequence is aligned with a space.
	* Matrix B assumes that the last element from second sequence is aligned with a space.
	* Matrix C assumes that the last element from first and second sequence are aligned together.
	* Matrix S contains best alignments from matrices A;B;C

	Each matrix cell (i,j) contains a pair of elements:

	* optimal alignment score of aligning subwords u[0..i], w[0..j],
	* trace with reference to previous matrix and indices.

	Break penalty function used in the algorithm is ``f(x) = -x -1``.

	Traceback collects the suballigments into the final optimal solution.

* NeedlemanWunsch

	The algorithm operates on a single matrix, where cell with indices (i,j) contains:

	* optimal solution for subwords u[0..i], w[0..j],
	* trace: Up, Left or Diagonal

	As in Gotoh, a traceback function builds the optimal solution.

	In ``NeedlemanWunsch`` module, there are also helper functions:

	* ``runScoreLastRow`` which is used by Hirschberg algorithm to get the last row of matrix, without the trace,
	* ``runGeneric`` which is used by malign algorithm to align 2 sequences of indices based on multialginment profiles.

* Hirschberg

	Hirschberg algorithm is based on Divide-and-Conquer technique. 
	Optimal alignment of 2 sequences is a concatenation of optimal alignments of specific pairs of subsequences.
	Based on that property, the algorithm splits the 2 sequences into smaller ones, solves the problem for those and concatenates the results.

	To split the sequences, ``runScoreLastRow`` function is used from ``NeedlemanWunsch`` module.
	The function is implemented recursivelly, using continuations technique to ensure tail recursive calls and prevent stack overflow.

* profile

	For given multialignment, the ``profile`` function computes a probability matrix [i,j] of occuring a specific Nucletoide i in column j.

* cons 

	For given multialignment, the ``consensusWord`` function computes a word that is most similar to the multialignment.

* malign 

	``alignByProfiles`` function aligns two given multialignments into one, 
	based on their profiles and the ``runGeneric`` function of ``NeedlemanWunsch`` module.

* UPGMA

	UPGMA algorithm starts with creating seperate cluster for each of the sequence from the input. 
	For each pair of the sequences, it calculates the similarity matrix using ``Hirschberg`` module.
	Then, the similarity is converted linearly to distance, so that two most similar sequences have the smallest distance.

	In each step, the UPGMA algorithm chooses two closest clusters and joins them (using ``alignByProfiles``), and recreates the distance matrix.
	The algorithm stops when there's only one cluster left.