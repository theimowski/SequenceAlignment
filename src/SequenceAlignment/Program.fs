module SequenceAlignment.Program

open System

let formatNucl = function
    | Break -> "-"
    | Nucl x -> 
        match x with
        | A -> "A"
        | C -> "C"
        | G -> "G"
        | T -> "T"

let formatSeq = List.map formatNucl >> String.concat ""

let formatOutput (alignment, sequence) = 
    let f,s = sequence |> List.unzip
    f |> formatSeq |> printfn "%s"
    s |> formatSeq |> printfn "%s"
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



[<EntryPoint>]
let main argv = 
    let p = fun (x:int) -> -1. - (1. * float x)
    let f, s = readInputSequence(), readInputSequence()
    let sim = readSimilarity()
    Gotoh.run((f,s), sim, p) |> formatOutput

    0