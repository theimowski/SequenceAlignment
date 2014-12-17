module SequenceAlignment.MultiAlign

open NeedlemanWunsch

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


let UPGMA(seqs : Sequence[], sim' : Similarity') : MultiAlignment = 
    
    let distances = 
        Array2D.init seqs.Length seqs.Length (fun i j -> 
            if j < i then Hirschberg.run (seqs.[i], seqs.[j], sim') |> fst
            else System.Double.NegativeInfinity)

    let maxSim = distances |> Seq.cast |> Seq.max

    distances |> Array2D.iteri (fun i j value -> 
                      distances.[i, j] <- if j < i then maxSim - value else System.Double.PositiveInfinity)

    let rec clusterize (maligns,distances) =
        match maligns with
        | [||] -> failwith "maligns should not be empty"
        | [|malign|] -> malign
        | _ -> 
            let i, j = 
                distances 
                |> Array2D.mapi (fun i j v-> (i,j),v) 
                |> Seq.cast<(int*int)*float>
                |> Seq.minBy snd 
                |> fst

            let malign = alignByProfiles(maligns.[j], maligns.[i],sim')
            let maligns' = 
                Array.concat [maligns.[0..j-1];[|malign|];maligns.[j+1..i-1];maligns.[i+1..]]
            let lenI, lenJ = Array2D.length1 maligns.[i] |> float,Array2D.length1 maligns.[j] |> float
            let dist (i,j) = let i,j = max i j, min i j in distances.[i,j]
            let distances' =
                seq { 
                    for row in 0..maligns.Length - 1 do
                        if row <> i then yield seq {
                            for col in 0..maligns.Length - 1 do
                                if col <> i then 
                                    match row,col with
                                    | row,col when row <= col -> yield System.Double.PositiveInfinity
                                    | row,_ when row = j ->
                                        yield (dist(i,col) * lenI + dist(j,col) * lenJ) / (lenI + lenJ)
                                    | _,col when col = j ->
                                        yield (dist(i,row) * lenI + dist(j,row) * lenJ) / (lenI + lenJ)
                                    | _ -> yield distances.[row,col]
                        }
                } |> array2D
            
            clusterize(maligns',distances')

    let initialMaligns = seqs |> Array.map (fun s  -> [|s |> Array.map Nucl|] |> array2D)
    clusterize (initialMaligns, distances)