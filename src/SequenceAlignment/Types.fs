[<AutoOpen>]
module SequenceAlignment.Types

open System

type Nucleotide = A | C | G | T
type Sequence = Nucleotide[]
type Nucleotide' = Nucl of Nucleotide | Break
type Alignment = float

type MultiAlignment = Nucleotide'[,]
type ConsensusWord = Nucleotide'[]
type Similarity' = Nucleotide' * Nucleotide' -> float
type MultiAlignmentProfile = Map<Nucleotide',float>[]

let mutable verbose = false
let logV fmt = Microsoft.FSharp.Core.Printf.kprintf (fun s -> printfn "%s" s) fmt 


let formatNucl = function
| Break -> "-"
| Nucl x -> 
    match x with
    | A -> "A"
    | C -> "C"
    | G -> "G"
    | T -> "T"

let formatSeq (x: seq<Nucleotide'>) = x |> Seq.map formatNucl |> String.concat ""

let formatMAlign (malign: MultiAlignment) = 
    [0..Array2D.length1 malign - 1]
    |> List.map (fun i -> malign.[i,*] |> formatSeq)
    |> String.concat Environment.NewLine
    |> printfn "%s"

let formatA2D (a: float[,]) = 
    [0..Array2D.length2 a - 1]
    |> List.map (fun j -> sprintf "%13d" j)
    |> List.append ["   "]
    |> String.concat ""
    |> printfn  "%s"

    [0..Array2D.length1 a - 1]
    |> List.map (fun i -> 
        a.[i,*] 
        |> Array.map (sprintf "|%10f|") 
        |> Array.append [|sprintf "%3d:" i|]
        |> String.concat " ")
    |> String.concat Environment.NewLine
    |> printfn "%s"