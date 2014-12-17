*Tomasz Heimowski, 131517*

*Filip Kucharczyk, 131548*

*Mateusz Szymański, 131633*

***

Dokumentacja
============

Wstęp
-----
Aplikacja *SequenceAlignment* zawiera rozwiązanie dwóch zadań projektowych z przedmiotu *Elementy Bioinformatyki*:

1. Zestawienie optymalne 2 sekwencji
	* algorytm z funkcją kary 
	* algorytm bez funkcji kary, ze złożonością czasową O(n)
2. Zestawienia wielu sekwencji
	* słowo konsensusowe, profil wielodopasowania
	* złożenie 2 wielodopasowań
	* progressive multi alignment wielu sekwencji

Aplikacja została napisana w języku F#. Jej kod źródłowy został udostępniony na witrynie *GitHub*, 
pod adresem [https://github.com/theimowski/SequenceAlignment](https://github.com/theimowski/SequenceAlignment)

Obsługa aplikacji
-----------------
Aplikacja jest typu konsolowego, obsługa wywołania z lini poleceń jest następująca:

	SequenceAlignment.exe <command> [-v|--verbose]

Dostępne komendy: 

* Gotoh

	Uruchamia algorytm zestawienia 2 sekwencji z funkcją kary ``f(x) = -x -1``, znany pod nazwą *Gotoh*.

* NeedlemanWunsch

	Uruchamia algorytm Needleman-Wunsch zestawienia 2 sekwencji bez funkcji kary.

* Hirschberg

	Uruchamia algorytm Hirschberg zestawienia 2 sekwencji bez funkcji kary.

* profile

	Wylicza profil wielodopasowania.

* cons 

	Wylicza słowo konsensusowe wielodopasowania.

* malign 

	Składa dwa wielodopasowania w jedno.

* UPGMA

	Uruchamia algorytm progressive multialignment dla wielu sekwencji.
	Funkcja wybierająca dwa klastry do złączenia to Unweighted Pair Grouping Method with Arithmetic Mean.

Flagi:

* [*-v|--verbose*]
	
	Wypisuje kolejne kroki algotymu na standardowe wyjście

Format danych
-------------

Opis algorytmów
---------------