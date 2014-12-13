module SequenceAlignment.MultiAlign

type MultiAlignment = Nucleotide'[,]
type ConsensusWord = Nucleotide'[]
type Similarity' = Nucleotide' * Nucleotide' -> float

let consensusWord 
    (malign : MultiAlignment, sim : Similarity') : ConsensusWord = 

    [|0 .. Array2D.length2 malign - 1|]
    |> Array.map (fun j -> malign.[*,j])
    |> Array.map (fun col -> 
                      [Break;Nucl A;Nucl C;Nucl G;Nucl T]
                      |> List.maxBy (fun n -> 
                                         col 
                                         |> Array.sumBy (fun c -> sim(n,c))))