module SequenceAlignment.Gotoh

open System

type Trace<'a> = Trace of ArrayCell<'a>[,] * int * int
and ArrayCell<'a> = 'a * option<Trace<'a>>
type AlignmentCell = ArrayCell<Alignment>

let p x = -1. - (1. * float x)

let highest(arrays: list<ArrayCell<'a>[,]>, i, j) : ArrayCell<'a> =
    arrays
    |> List.map (fun a -> fst a.[i,j], a)
    |> List.maxBy fst
    |> (fun (value,array) -> value, Some(Trace (array,i,j)))

let inline addFst x (f,s) = f + x, s

let private printState (a,b,c,s) =
    if Types.verbose then
        Console.Clear()
        printfn "Gotoh arrays state"
        printfn "A:"
        a |> Array2D.map fst |> formatA2D
        printfn "B:"
        b |> Array2D.map fst |> formatA2D
        printfn "C:"
        c  |> Array2D.map fst |> formatA2D
        printfn "S:"
        s  |> Array2D.map fst |> formatA2D
        printfn "%s" Environment.NewLine
        do Threading.Thread.Sleep 100

let run 
    (fstSeq : Sequence, sndSeq : Sequence, sim : Similarity')
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
        fst s.[i-1,j-1] + sim(Nucl fstSeq.[i-1], Nucl sndSeq.[j-1]), Some(Trace(c,i,j))

    s.[0,0] <- 0., None
    for i in 1..len1 do 
        b.[i,0] <- p(i), None
        s.[i,0] <- p(i), Some(Trace(b,i,0))
    for i in 1..len2 do
        a.[0,i] <- p(i), None
        s.[0,i] <- p(i), Some(Trace(a,0,i))
    for i in 1..len1 do
        for j in 1..len2 do
            a.[i,j] <- countA(i,j)
            b.[i,j] <- countB(i,j)
            c.[i,j] <- countC(i,j)
            s.[i,j] <- highest([a;b;c],i,j)
            printState(a,b,c,s)

    let rec traceBack (acc,trail : option<Trace<Alignment>>) =
        match trail with
        | Some(Trace(array,i,j)) when array.Equals(a) -> 
            let nextTrail = snd a.[i,j]
            let from = match nextTrail with Some(Trace(_,_,j)) -> j+1 | None -> 1
            let indels = [from..j] |> List.map (fun j -> Break, Nucl sndSeq.[j-1])
            traceBack (indels @ acc, nextTrail)
        | Some(Trace(array,i,j)) when array.Equals(b) -> 
            let nextTrail = snd b.[i,j]
            let from = match nextTrail with Some(Trace(_,i,_)) -> i+1 | None -> 1
            let indels = [from..i] |> List.map (fun i -> Nucl fstSeq.[i-1], Break)
            traceBack (indels @ acc, nextTrail)
        | Some(Trace(array,i,j)) when array.Equals(c) -> 
            traceBack ((Nucl fstSeq.[i-1],Nucl sndSeq.[j-1])::acc, snd s.[i-1,j-1])
        | None -> acc
        | _ -> failwith "unexpected trail"
    
    let sequence = traceBack([], snd s.[len1,len2])
    fst s.[len1,len2], sequence