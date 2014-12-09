﻿module SequenceAlignment

open System

type Nucleotide = A | C | G | T
type Sequence = Nucleotide[]
type Nucleotide' = Nucl of Nucleotide | Break
type BreakPenalty = int -> float
type Similarity = Nucleotide * Nucleotide -> float
type Alignment = float

type Trail<'a> = Trail of ArrayCell<'a>[,] * int * int
and ArrayCell<'a> = 'a * option<Trail<'a>>
type AlignmentCell = ArrayCell<Alignment>

let highest(arrays: list<ArrayCell<'a>[,]>, i, j) : ArrayCell<'a> =
    arrays
    |> List.map (fun a -> fst a.[i,j], a)
    |> List.maxBy fst
    |> (fun (value,array) -> value, Some(Trail (array,i,j)))

let inline addFst x (f,s as tuple) = f + x, s

let alignTwoWithPenalty 
    ((fstSeq : Sequence, sndSeq : Sequence), sim : Similarity, p : BreakPenalty)
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

let formatNucl = function
    | Break -> "-"
    | Nucl x -> 
        match x with
        | A -> "A"
        | C -> "C"
        | G -> "G"
        | T -> "T"

let formatOutput (alignment, sequence) = 
    let f,s = sequence |> List.unzip
    f |> List.map formatNucl |> String.concat "" |> printfn "%s"
    s |> List.map formatNucl |> String.concat "" |> printfn "%s"
    printfn "similarity: %f" alignment

let readInputSequence() = 
    let parse = function
    | 'A' -> Some A
    | 'C' -> Some C
    | 'G' -> Some G
    | 'T' -> Some T
    | _ -> None
    let line = Console.ReadLine()
    line |> Seq.choose parse |> Seq.toArray

let readSimilarity() : Similarity =
    let parseLine (s:string) = s.Split([|';'|]) |> Array.map float
    let lookup = 
        [A;C;G;T] 
        |> List.map (fun n -> n, Console.ReadLine() |> parseLine)
        |> Map.ofList
    (fun (f,s) ->
        let f',s' = max f s, min f s
        let index = match s' with A -> 0 | C -> 1 | G -> 2 | T -> 3
        lookup.[f'].[index])


let NeedlemanWunsch
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

let sim (a,b) = if a = b  then 2. else -1.
NeedlemanWunsch([|A;G;T;A;C;G;C;A|],[|T;A;T;G;C|],sim,-2.)

let inline lastRow (x: _[,]) = x.[Array2D.length1 x-1,*]

let inline splitBefore i (x: _[]) = 
    match i with 
    | _ when i < 1 -> [||], x
    | _ when i > x.Length-1 -> x, [||]
    | _ -> x.[0..i-1], x.[i..(x.Length-1)]

let Hirschberg
    (fstSeq : Sequence, sndSeq: Sequence,sim : Similarity, indelCost : float)
    : Alignment * list<Nucleotide' * Nucleotide'> =

    let rec Hirschberg' (fstSeq,sndSeq,cont) = 

        match fstSeq,sndSeq with
        | [|f|], [|s|] -> cont(sim(f,s), [Nucl f,Nucl s])
        | [||], _ -> 
            cont(indelCost * float sndSeq.Length, 
                 sndSeq |> Array.map (fun s -> Break, Nucl s) |> List.ofArray)
        | _, [||] -> 
            cont(indelCost * float fstSeq.Length, 
                 fstSeq |> Array.map (fun f -> Nucl f, Break) |> List.ofArray)
        | _ -> 
            let imid = fstSeq.Length / 2
            let f1,f2 = fstSeq |> splitBefore imid
        
            let upper = 
                NeedlemanWunsch(f1,sndSeq,sim,indelCost) 
                |> lastRow
            let lower =
                NeedlemanWunsch(Array.rev f2,
                                Array.rev sndSeq, sim, indelCost)
                |> lastRow

            let jmid = 
                (upper,lower |> Array.rev) 
                ||> Array.map2 (+)
                |> Array.mapi (fun i v -> i,v)
                |> Array.maxBy snd
                |> fst
            
            let s1,s2 = sndSeq |> splitBefore jmid
            
            Hirschberg' (f1,s1,(fun (lA,lL) ->
                Hirschberg' (f2,s2,(fun (rA,rL) ->
                    cont(lA+rA,lL@rL)))))

    Hirschberg' (fstSeq,sndSeq,id)

Hirschberg([|A;G;T;A;C;G;C;A|],[|T;A;T;G;C|],sim,-2.) |> formatOutput

[<EntryPoint>]
let main argv = 
    let p = fun (x:int) -> -1. - (1. * float x)
    let f, s = readInputSequence(), readInputSequence()
    let sim = readSimilarity()
    alignTwoWithPenalty((f,s), sim, p) |> formatOutput

    0