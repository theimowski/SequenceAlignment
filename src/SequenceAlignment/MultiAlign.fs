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

open NeedlemanWunsch

module Array2D =
    let transpose(a : _[,]) = 
        [0..Array2D.length2 a - 1]
        |> List.map (fun j -> a.[*,j])
        |> array2D

let alignByProfiles 
    (malign1 : MultiAlignment, malign2 : MultiAlignment, sim : Similarity') 
    : MultiAlignment = 
    let p,q = profile malign1, profile malign2

    let alphabet = [Nucl A;Nucl C;Nucl G;Nucl T]
    let alphabet' = alphabet @ [Break]

    let ops = {
        LeftIndelCost = 
            fun j -> alphabet |> List.sumBy (fun n -> sim(n,Break) * q.[j].[n])
        UpIndelCost = 
            fun i -> alphabet |> List.sumBy (fun n -> sim(n,Break) * p.[i].[n])
        DiagonalCost = 
            fun (i,j) -> 
                alphabet' |> List.sumBy (fun n1 ->
                    alphabet' |> List.sumBy (fun n2 -> sim(n1,n2) * p.[i].[n1] * q.[j].[n2]))
    }

    let _, list = 
        NeedlemanWunsch.runGeneric([|0..Array2D.length2 malign1 - 1|],
                             [|0..Array2D.length2 malign2 - 1|],
                             ops)

    let len1,len2 = Array2D.length1 malign1, Array2D.length1 malign2

    list
    |> List.map (fun (f,s) ->
                    Array.append
                        (match f with
                        | Some f -> malign1.[*,f]
                        | None -> Array.create len1 Break)
                        (match s with
                        | Some s -> malign2.[*,s]
                        | None -> Array.create len2 Break))
    |> array2D
    |> Array2D.transpose

let score (malign : MultiAlignment, sim : Similarity') : float = 
    [0..Array2D.length2 malign - 1]
    |> List.map (fun j -> malign.[*,j])
    |> List.map (fun col -> 
                     col |> Array.mapi (fun i n -> 
                                            col.[i+1..col.Length-1]
                                            |> Array.sumBy (fun n2 -> sim(n,n2) )))
    |> List.sumBy Array.sum


/// http://stackoverflow.com/a/19406376/1397724
let concat (a1: 'a[,]) (a2: 'a[,]) =
    let a1l1,a1l2,a2l1,a2l2 = (Array2D.length1 a1),(Array2D.length2 a1),(Array2D.length1 a2),(Array2D.length2 a2)
    if a1l2 <> a2l2 then failwith "arrays have different column sizes"
    let result = Array2D.zeroCreate (a1l1 + a2l1) a1l2
    Array2D.blit a1 0 0 result 0 0 a1l1 a1l2
    Array2D.blit a2 0 0 result a1l1 0 a2l1 a2l2
    result

let UPGMA(seqs : Sequence[], sim' : Similarity') : MultiAlignment = 
    
    let sim(x,y) = sim'(Nucl x, Nucl y)
    //////////// !!!!!!!!!!
    let magicValueDoSthWithIT = -2.
    let distances = 
        Array2D.init seqs.Length seqs.Length (fun i j -> 
            if j < i then Hirschberg.run (seqs.[i], seqs.[j], sim, magicValueDoSthWithIT) |> fst
            else System.Double.PositiveInfinity)

    let maxSim = distances |> Seq.cast |> Seq.max

    distances |> Array2D.iteri (fun i j value -> 
                     if j < i then distances.[i, j] <- maxSim - value)

    let rec clusterize (maligns,distances) =
        match maligns with
        | [||] -> failwith "maligns should not be empty"
        | [|malign|] -> malign
        | _ -> 
            let i, j = 
                distances 
                |> Array2D.mapi (fun i j v-> (i,j),v) 
                |> Seq.cast 
                |> Seq.maxBy snd 
                |> fst

            let malign = alignByProfiles(maligns.[i], maligns.[j],sim')
            let maligns' = 
                Array.concat [maligns.[..i-1];[|malign|];maligns.[i..j-1];maligns.[j..]]
            let lenI, lenJ = Array2D.length1 maligns.[i] |> float,Array2D.length1 maligns.[j] |> float
            let distances' =
                seq { 
                    for row in 0..maligns.Length - 1 do
                        yield seq {
                            for col in 0..maligns.Length - 1 do
                                match row,col with
                                | row,col when row <= col -> yield System.Double.PositiveInfinity
                                | _, col when col = j -> ()
                                | row,_ when row = j -> ()
                                | row,_ when row = i ->
                                    yield (distances.[i,col] * lenI + distances.[j,col] * lenJ) /
                                            (lenI + lenJ)
                                | _,col when col = i ->
                                    yield (distances.[row,i] * lenI + distances.[row,j] * lenJ) /
                                            (lenI + lenJ)
                                | _ -> yield distances.[row,col]
                        }
                } |> array2D
            
            clusterize(maligns',distances')

    let initialMaligns = seqs |> Array.map (fun s  -> [|s |> Array.map Nucl|] |> array2D)
    clusterize (initialMaligns, distances)

let seqs = [|
    [|G;C;T;T|]
    [|A;C;T;T|]
|]
let sim (x,y) = if x = y then 1. else 0.
UPGMA(seqs,sim)