[<AutoOpen>]
module SequenceAlignment.Types

type Nucleotide = A | C | G | T
type Sequence = Nucleotide[]
type Nucleotide' = Nucl of Nucleotide | Break
type BreakPenalty = int -> float
type Similarity = Nucleotide * Nucleotide -> float
type Alignment = float