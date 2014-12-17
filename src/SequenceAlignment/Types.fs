[<AutoOpen>]
module SequenceAlignment.Types

type Nucleotide = A | C | G | T
type Sequence = Nucleotide[]
type Nucleotide' = Nucl of Nucleotide | Break
type BreakPenalty = int -> float
type Similarity = Nucleotide * Nucleotide -> float
type Alignment = float

type MultiAlignment = Nucleotide'[,]
type ConsensusWord = Nucleotide'[]
type Similarity' = Nucleotide' * Nucleotide' -> float
type MultiAlignmentProfile = Map<Nucleotide',float>[]

let mutable verbose = false
let logV fmt = Microsoft.FSharp.Core.Printf.kprintf (fun s -> printfn "%s" s) fmt 