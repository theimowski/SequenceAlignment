module SequenceAlignment.Tests

open Xunit

open FsCheck
open FsCheck.Xunit

let p = fun (x:int) -> -1. - (1. * float x)
let sim (x,y) = if x = y then 2. else 0.

[<Property>]
let ``removing breaks gives input`` (fstSeq : Sequence, sndSeq : Sequence) =
    let _,sequence =  alignTwoWithPenalty((fstSeq,sndSeq), sim, p)
    let f,s = sequence |> List.unzip
    f |> List.filter ((<>) Break) = (fstSeq |> Array.toList |> List.map Nucl) &&
    s |> List.filter ((<>) Break) = (sndSeq |> Array.toList |> List.map Nucl)

[<Property>]
let ``alignment is correct`` (fstSeq : Sequence, sndSeq : Sequence) =
    let alignment,sequence =  alignTwoWithPenalty((fstSeq,sndSeq), sim, p)
    let rec count (sequence,breakLength,sum) =
        match sequence with 
        | [] -> 
            let penalty = if breakLength > 0 then p(breakLength) else 0.
            sum + penalty
        | (Break,_) :: t | (_,Break) :: t ->
            count(t,breakLength+1,sum)
        | (f,s) :: t ->
            let penalty = if breakLength > 0 then p(breakLength) else 0.
            count(t,0,sum + penalty + sim(f,s))

    count(sequence,0,0.) = alignment

[<Property>]
let ``no double breaks`` (fstSeq : Sequence, sndSeq : Sequence) =
    let _,sequence =  alignTwoWithPenalty((fstSeq,sndSeq), sim, p)
    sequence
    |> List.forall (fun (f,s) -> not (f = Break && s = Break))