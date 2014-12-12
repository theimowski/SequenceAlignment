[<AutoOpen>]
module SequenceAlignment.Types

type Nucleotide = A | C | G | T
type Sequence = Nucleotide[]
type Nucleotide' = Nucl of Nucleotide | Break
type BreakPenalty = int -> float
type Similarity = Nucleotide * Nucleotide -> float
type Alignment = float

let mutable verbose = false
let logV fmt = Microsoft.FSharp.Core.Printf.kprintf (fun s -> if verbose then printfn "%s" s) fmt 