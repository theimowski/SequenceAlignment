open System

type Nucleotide = A | C | G | T
type Sequence = Nucleotide[]
type Nucleotide' = Nucl of Nucleotide | Break
type BreakPenalty = int -> float
type Similarity = Nucleotide * Nucleotide -> float
type Alignment = float

type private Trail = char * int * int

let alignTwoWithPenalty 
    ((fstSeq : Sequence, sndSeq : Sequence), sim : Similarity, p : BreakPenalty)
    : Alignment * list<Nucleotide' * Nucleotide'> =
        
    let len1, len2 = fstSeq.Length, sndSeq.Length
    let initA2D() : (float*option<Trail>)[,] = 
        Array2D.create (len1+1) (len2+1) (Double.NegativeInfinity, None)
    let a,b,c,s = initA2D(), initA2D(), initA2D(), initA2D()
    let countA(i,j) = 
        [0..j-1] 
        |> List.map (fun k -> 
            let b,c = fst b.[i,k], fst c.[i,k]
            if b > c then b + p(j-k), Some('b',i,k)
            else c + p(j-k), Some('c',i,k))
        |> List.maxBy fst
    let countB(i,j) = 
        [0..i-1] 
        |> List.map (fun k -> 
            let a,c = fst a.[k,j], fst c.[k,j]
            if a > c then a + p(i-k), Some('a',k,j)
            else c + p(i-k), Some('c',k,j))
        |> List.maxBy fst
    let countC(i,j) = fst s.[i-1,j-1] + sim(fstSeq.[i-1], sndSeq.[j-1]), Some('c',i,j)

    s.[0,0] <- 0., None
    for i in 1..len1 do 
        b.[i,0] <- p(i), None
        s.[i,0] <- p(i), None
    for i in 1..len2 do
        a.[0,i] <- p(i), None
        s.[0,i] <- p(i), None
    for i in 1..len1 do
        for j in 1..len2 do
            a.[i,j] <- countA(i,j)
            b.[i,j] <- countB(i,j)
            c.[i,j] <- countC(i,j)
            
            let (v,t) = 
                [fst a.[i,j], 'a'
                 fst b.[i,j], 'b'
                 fst c.[i,j], 'c' 
                ] 
                |> List.maxBy fst

            s.[i,j] <- v, Some(t,i,j)

    let rec x (acc,trail : option<Trail>) =
        printfn "%O" (trail)  
        match trail with
        | Some('a',i,j) -> 
            let j_prev = match snd a.[i,j] with Some(_,_,j) -> j | None -> j
            let indels = 
                [j_prev+1..j]
                |> List.map (fun j -> Break, Nucl sndSeq.[j-1])
            x (indels @ acc, snd a.[i,j])
        | Some('b',i,j) -> 
            let i_prev = match snd b.[i,j] with Some(_,i,_) -> i | None -> i
            let indels = 
                [i_prev+1..i]
                |> List.map (fun i -> Nucl fstSeq.[i-1], Break)
            x (indels @ acc, snd b.[i,j])
        | Some('c',i,j) -> 
            x ((Nucl fstSeq.[i-1],Nucl sndSeq.[j-1])::acc, snd s.[i-1,j-1])
        | None -> acc
        | _ -> failwith "unexpected trail"
    
    
    let sequence = x([], snd s.[len1,len2])
    fst s.[len1,len2], sequence

let fstSeq = [|C;A;G;C|]//;C|]
let sndSeq = [|C|]//;C;T;G|]
let p = fun (x:int) -> -1. - (1. * float x)
let sim (x,y) = if x = y then 2. else 0.

alignTwoWithPenalty((fstSeq,sndSeq), sim, p)