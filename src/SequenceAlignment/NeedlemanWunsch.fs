module SequenceAlignment.NeedlemanWunsch

type Trace = Left | Diagonal | Up
type ArrayCell = Alignment * option<Trace>

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