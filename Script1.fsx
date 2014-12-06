open System

type Nucleotide = A | C | G | T
type Sequence = Nucleotide[]
type Nucleotide' = Nucl of Nucleotide | Break
type Sequence' = Nucleotide'[]
type BreakPenalty = int -> float
type Similarity = Nucleotide * Nucleotide -> float
type Alignment = float

type private Trail = char * int * int

let alignTwoWithPenalty 
    ((fstSeq : Sequence, sndSeq : Sequence), sim : Similarity, p : BreakPenalty)
    : (Sequence' * Sequence') * Alignment =
        
    let len1, len2 = fstSeq.Length, sndSeq.Length
    let initA2D() = Array2D.create (len1+1) (len2+1) Double.NegativeInfinity
    let initA2DwithTrail() : (float*option<Trail>)[,] = 
        Array2D.create (len1+1) (len2+1) (Double.NegativeInfinity, None)
    let a,b = initA2DwithTrail(), initA2DwithTrail()
    let c,s = initA2D(), initA2D()
    let countA(i,j) = 
        [0..j-1] 
        |> List.map (fun k -> 
            let b,c = fst b.[i,k], c.[i,k]
            if b > c then b + p(j-k), Some('b',i,k)
            else c + p(j-k), Some('c',i,k))
        |> List.maxBy fst
    let countB(i,j) = 
        [0..i-1] 
        |> List.map (fun k -> 
            let a,c = fst a.[k,j], c.[k,j]
            if a > c then a + p(i-k), Some('a',k,j)
            else c + p(i-k), Some('c',k,j))
        |> List.maxBy fst
    let countC(i,j) = s.[i-1,j-1] + sim(fstSeq.[i-1], sndSeq.[j-1])

    s.[0,0] <- 0.
    for i in 1..len1 do 
        b.[i,0] <- p(i), None
        s.[i,0] <- p(i)
    for i in 1..len2 do
        a.[0,i] <- p(i), None
        s.[0,i] <- p(i)
    for i in 1..len1 do
        for j in 1..len2 do
            a.[i,j] <- countA(i,j)
            b.[i,j] <- countB(i,j)
            c.[i,j] <- countC(i,j)
            s.[i,j] <- List.max [fst a.[i,j]; fst b.[i,j];c.[i,j]]

    ([||],[||]), s.[len1,len2]

let fstSeq = [|C;A;G;C;C|]
let sndSeq = [|C;C;T;G|]
let p = fun (x:int) -> -1. - (1. * float x)
let sim (x,y) = if x = y then 2. else 0.

alignTwoWithPenalty((fstSeq,sndSeq), sim, p)