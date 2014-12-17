module SequenceAlignment.NeedlemanWunsch

open System

let runScore
    (fstSeq : Sequence, sndSeq: Sequence, sim' : Similarity') : Alignment[,] =
    let len1,len2 = fstSeq.Length, sndSeq.Length
    let array = Array2D.zeroCreate (len1+1) (len2+1)
    [1..len1] |> List.iter (fun i -> array.[i,0] <- array.[i-1,0] + sim'(Nucl fstSeq.[i-1], Break))
    [1..len2] |> List.iter (fun j -> array.[0,j] <- array.[0,j-1] + sim'(Break, Nucl sndSeq.[j-1]))
    for i in 1..len1 do
        for j in 1..len2 do
            array.[i,j] <- [array.[i-1,j-1] + sim'(Nucl fstSeq.[i-1],Nucl sndSeq.[j-1])
                            array.[i-1,j] + sim'(Nucl fstSeq.[i-1], Break)
                            array.[i,j-1] + sim'(Break, Nucl sndSeq.[j-1])
                            ] |> List.max
    array

let runScoreLastRow
    (fstSeq : Sequence, sndSeq: Sequence, sim' : Similarity') : Alignment[] =
    let len1,len2 = fstSeq.Length, sndSeq.Length
    let prev_row = Array.zeroCreate (len2+1)
    [1..len2] |> List.iter (fun j -> prev_row.[j] <- prev_row.[j-1] + sim'(Break, Nucl sndSeq.[j-1]))
    let cur_row = Array.copy prev_row
    for i in 1..len1 do 
        cur_row.[0] <- prev_row.[0] + sim'(Nucl fstSeq.[i-1], Break)
        for j in 1..len2 do
            cur_row.[j] <- [prev_row.[j-1] + sim'(Nucl fstSeq.[i-1],Nucl sndSeq.[j-1])
                            prev_row.[j] + sim'(Nucl fstSeq.[i-1], Break)
                            cur_row.[j-1]+ sim'(Break, Nucl sndSeq.[j-1])
                            ] |> List.max
        Array.blit cur_row 0 prev_row 0 (len2+1)
            
    cur_row

type Trace = NoTrace | Left | Diagonal | Up
type ArrayCell = Alignment * Trace
type RunOperations = {
    LeftIndelCost : int -> float
    UpIndelCost   : int -> float
    DiagonalCost : int * int -> float
}

let printState(a) =
    let formatT = function | NoTrace -> "." | Up -> "^" | Left -> "<" | Diagonal  -> "\\"
    if Types.verbose then
        Console.Clear()
        logV "Needleman-Wunsch state:"
        logV "%A" (a |> Array2D.map (fun (d,t) -> formatT t, d))
        Threading.Thread.Sleep(300)

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
            printState(array)

    let rec traceBack (i,j,acc) =
        let trace = snd array.[i,j]
        match trace with
        | NoTrace -> acc
        | Diagonal -> traceBack(i-1, j-1, (Some fstSeq.[i-1], Some sndSeq.[j-1]) :: acc)
        | Up ->       traceBack(i-1, j, (Some fstSeq.[i-1], None) :: acc)
        | Left -> traceBack(i, j-1, (None, Some sndSeq.[j-1]) :: acc)

    fst array.[len1,len2], traceBack(len1, len2, [])

let run
    (fstSeq : Sequence, sndSeq: Sequence, sim' : Similarity') 
    : Alignment * list<Nucleotide' * Nucleotide'> = 

    let ops = {
        LeftIndelCost = fun j -> sim'(Break, Nucl sndSeq.[j])
        UpIndelCost = fun i -> sim'(Nucl fstSeq.[i], Break)
        DiagonalCost = fun (i,j) -> sim'(Nucl fstSeq.[i],Nucl sndSeq.[j])
    }

    let translate = function
        | Some n -> Nucl n
        | None -> Break
    let a,l = runGeneric(fstSeq, sndSeq, ops)
    a,
    l |> List.map (fun (x,y) -> translate x, translate y)