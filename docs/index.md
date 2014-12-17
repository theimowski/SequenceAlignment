
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

Description of the algorithms
-----------------------------