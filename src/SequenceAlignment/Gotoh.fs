﻿module SequenceAlignment.Gotoh

open System

type Trail<'a> = Trail of ArrayCell<'a>[,] * int * int
and ArrayCell<'a> = 'a * option<Trail<'a>>
type AlignmentCell = ArrayCell<Alignment>

let highest(arrays: list<ArrayCell<'a>[,]>, i, j) : ArrayCell<'a> =
    arrays
    |> List.map (fun a -> fst a.[i,j], a)
    |> List.maxBy fst
    |> (fun (value,array) -> value, Some(Trail (array,i,j)))

let inline addFst x (f,s) = f + x, s

let private printState (a,b,c,s) =
    let format = fst >> string >> (function | "-Infinity" -> "-i" | s -> s) >> sprintf "%5s"
    if Types.verbose then
        Console.Clear()
        logV "Gotoh arrays state"
        logV "A:"
        logV "%A" (a |> Array2D.map format)
        logV "B:"
        logV "%A" (b |> Array2D.map format)
        logV "C:"
        logV "%A" (c  |> Array2D.map format)
        logV "S:"
        logV "%A" (s  |> Array2D.map format)
        Threading.Thread.Sleep(300)

let run 
    (fstSeq : Sequence, sndSeq : Sequence, sim : Similarity, p : BreakPenalty)
    : Alignment * list<Nucleotide' * Nucleotide'> =
        
    let len1, len2 = fstSeq.Length, sndSeq.Length
    let initA2D() : AlignmentCell[,] = 
        Array2D.create (len1+1) (len2+1) (Double.NegativeInfinity, None)
    let a,b,c,s = initA2D(), initA2D(), initA2D(), initA2D()
    let countA(i,j) : AlignmentCell = 
        [0..j-1] 
        |> List.map (fun k -> highest([b;c], i, k) |> addFst (p(j-k)))
        |> List.maxBy fst
    let countB(i,j) : AlignmentCell = 
        [0..i-1] 
        |> List.map (fun k -> highest([a;c], k, j) |> addFst (p(i-k)))
        |> List.maxBy fst
    let countC(i,j) : AlignmentCell =
        fst s.[i-1,j-1] + sim(fstSeq.[i-1], sndSeq.[j-1]), Some(Trail(c,i,j))

    s.[0,0] <- 0., None
    for i in 1..len1 do 
        b.[i,0] <- p(i), None
        s.[i,0] <- p(i), Some(Trail(b,i,0))
    for i in 1..len2 do
        a.[0,i] <- p(i), None
        s.[0,i] <- p(i), Some(Trail(a,0,i))
    for i in 1..len1 do
        for j in 1..len2 do
            a.[i,j] <- countA(i,j)
            b.[i,j] <- countB(i,j)
            c.[i,j] <- countC(i,j)
            s.[i,j] <- highest([a;b;c],i,j)
            printState(a,b,c,s)

    let rec unfoldTrail (acc,trail : option<Trail<Alignment>>) =
        match trail with
        | Some(Trail(array,i,j)) when array.Equals(a) -> 
            let nextTrail = snd a.[i,j]
            let from = match nextTrail with Some(Trail(_,_,j)) -> j+1 | None -> 1
            let indels = [from..j] |> List.map (fun j -> Break, Nucl sndSeq.[j-1])
            unfoldTrail (indels @ acc, nextTrail)
        | Some(Trail(array,i,j)) when array.Equals(b) -> 
            let nextTrail = snd b.[i,j]
            let from = match nextTrail with Some(Trail(_,i,_)) -> i+1 | None -> 1
            let indels = [from..i] |> List.map (fun i -> Nucl fstSeq.[i-1], Break)
            unfoldTrail (indels @ acc, nextTrail)
        | Some(Trail(array,i,j)) when array.Equals(c) -> 
            unfoldTrail ((Nucl fstSeq.[i-1],Nucl sndSeq.[j-1])::acc, snd s.[i-1,j-1])
        | None -> acc
        | _ -> failwith "unexpected trail"
    
    let sequence = unfoldTrail([], snd s.[len1,len2])
    fst s.[len1,len2], sequence