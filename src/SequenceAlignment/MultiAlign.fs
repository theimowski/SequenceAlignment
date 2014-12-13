module SequenceAlignment.MultiAlign

type MultiAlignment = Nucleotide'[,]
type ConsensusWord = Nucleotide'[]
type Similarity' = Nucleotide' * Nucleotide' -> float
type MultiAlignmentProfile = Map<Nucleotide',float>[]

let consensusWord (malign : MultiAlignment, sim : Similarity') : ConsensusWord = 
    [| 0..Array2D.length2 malign - 1 |]
    |> Array.map (fun j -> malign.[*, j])
    |> Array.map (fun col -> 
           [ Break
             Nucl A
             Nucl C
             Nucl G
             Nucl T ]
           |> List.maxBy (fun n -> col |> Array.sumBy (fun c -> sim (n, c))))

let profile (malign : MultiAlignment) : MultiAlignmentProfile = 
    let len1 = Array2D.length1 malign

    let emptyMap = 
        [ Nucl A, 0.
          Nucl C, 0.
          Nucl G, 0.
          Nucl T, 0.
          Break, 0. ]
        |> Map.ofList

    [| 0..Array2D.length2 malign - 1 |]
    |> Array.map (fun j -> malign.[*, j])
    |> Array.map (fun col -> 
           col
           |> Seq.countBy id
           |> Seq.map (fun (k, c) -> k, float c / float len1)
           |> Seq.fold (fun s (n,value) -> s |> Map.add n value) emptyMap)