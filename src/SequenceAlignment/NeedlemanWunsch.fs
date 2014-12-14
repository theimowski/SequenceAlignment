module SequenceAlignment.NeedlemanWunsch

open System

let runScore
    (fstSeq : Sequence, sndSeq: Sequence, sim : Similarity, indelCost : float) 
    : Alignment[,] =
    let len1,len2 = fstSeq.Length, sndSeq.Length
    let array = Array2D.zeroCreate (len1+1) (len2+1)
    [1..len1] |> List.iter (fun i -> array.[i,0] <- float i * indelCost)
    [1..len2] |> List.iter (fun j -> array.[0,j] <- float j * indelCost)
    for i in 1..len1 do
        for j in 1..len2 do
            array.[i,j] <- [array.[i-1,j-1] + sim(fstSeq.[i-1],sndSeq.[j-1])
                            array.[i-1,j]+indelCost
                            array.[i,j-1]+indelCost
                            ] |> List.max
    array

let runScoreLastRow
    (fstSeq : Sequence, sndSeq: Sequence, sim : Similarity, indelCost : float) 
    : Alignment[] =
    let len1,len2 = fstSeq.Length, sndSeq.Length
    let prev_row = Array.init (len2+1) (fun j -> float j * indelCost)
    let cur_row = Array.copy prev_row
    for i in 1..len1 do 
        cur_row.[0] <- float i * indelCost
        for j in 1..len2 do
            cur_row.[j] <- [prev_row.[j-1] + sim(fstSeq.[i-1],sndSeq.[j-1])
                            prev_row.[j] + indelCost
                            cur_row.[j-1]+indelCost
                            ] |> List.max
        Array.blit cur_row 0 prev_row 0 (len2+1)
            
    cur_row

type Trace = NoTrace | Left | Diagonal | Up
type ArrayCell = Alignment * Trace

let printState(a) =
    let formatT = function | NoTrace -> "." | Up -> "^" | Left -> "<" | Diagonal  -> "\\"
    if Types.verbose then
        Console.Clear()
        logV "Needleman-Wunsch state:"
        logV "%A" (a |> Array2D.map (fun (d,t) -> formatT t, d))
        Threading.Thread.Sleep(300)

type RunOperations = {
    LeftIndelCost : int -> float
    UpIndelCost   : int -> float
    DiagonalCost : int * int -> float
}

let runGeneric
    (fstSeq : 'a[], sndSeq : 'a[], ops : RunOperations)
    : float * list<option<'a> * option<'a>> = 
    let len1,len2 = fstSeq.Length, sndSeq.Length
    let array = Array2D.create (len1+1) (len2+1) (0., NoTrace)
    [1..len1] |> List.iter (fun i -> array.[i,0] <- fst array.[i-1,0] + ops.UpIndelCost(i-1), Up)
    [1..len2] |> List.iter (fun j -> array.[0,j] <- fst array.[0,j-1] + ops.LeftIndelCost(j-1), Left)
    for i in 1..len1 do
        for j in 1..len2 do
            array.[i,j] <- [fst array.[i-1,j-1] + ops.DiagonalCost(i-1,j-1), Diagonal
                            fst array.[i-1,j] + ops.UpIndelCost(i-1), Up
                            fst array.[i,j-1]+ ops.LeftIndelCost(j-1), Left
                            ] |> List.maxBy fst

    let rec traceBack (i,j,acc) =
        let trace = snd array.[i,j]
        match trace with
        | NoTrace -> acc
        | Diagonal -> traceBack(i-1, j-1, (Some fstSeq.[i-1], Some sndSeq.[j-1]) :: acc)
        | Up ->       traceBack(i-1, j, (Some fstSeq.[i-1], None) :: acc)
        | Left -> traceBack(i, j-1, (None, Some sndSeq.[j-1]) :: acc)

    fst array.[len1,len2], traceBack(len1, len2, [])

let run
    (fstSeq : Sequence, sndSeq: Sequence, sim : Similarity, indelCost : float) 
    : Alignment * list<Nucleotide' * Nucleotide'> = 

    let ops = {
        LeftIndelCost = fun _ -> indelCost
        UpIndelCost = fun _ -> indelCost
        DiagonalCost = fun (i,j) -> sim(fstSeq.[i],sndSeq.[j])
    }

    let translate = function
        | Some n -> Nucl n
        | None -> Break
    let a,l = runGeneric(fstSeq, sndSeq, ops)
    a,
    l |> List.map (fun (x,y) -> translate x, translate y)

//    let len1,len2 = fstSeq.Length, sndSeq.Length
//    let array = Array2D.create (len1+1) (len2+1) (0., NoTrace)
//    [1..len1] |> List.iter (fun i -> array.[i,0] <- float i * indelCost, Up)
//    [1..len2] |> List.iter (fun j -> array.[0,j] <- float j * indelCost, Left)
//    for i in 1..len1 do
//        for j in 1..len2 do
//            array.[i,j] <- [fst array.[i-1,j-1] + sim(fstSeq.[i-1],sndSeq.[j-1]), Diagonal
//                            fst array.[i-1,j] + indelCost, Up
//                            fst array.[i,j-1]+ indelCost, Left
//                            ] |> List.maxBy fst
//            printState(array)
//    
//    let rec traceBack (i,j,acc) =
//        let trace = snd array.[i,j]
//        match trace with
//        | NoTrace -> acc
//        | Diagonal -> traceBack(i-1, j-1, (Nucl fstSeq.[i-1], Nucl sndSeq.[j-1]) :: acc)
//        | Up ->       traceBack(i-1, j, (Nucl fstSeq.[i-1], Break) :: acc)
//        | Left -> traceBack(i, j-1, (Break, Nucl sndSeq.[j-1]) :: acc)
//
//    fst array.[len1,len2], traceBack(len1, len2, [])