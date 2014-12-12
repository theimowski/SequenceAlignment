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

let parse = function
    | 'A' -> Some A
    | 'C' -> Some C
    | 'G' -> Some G
    | 'T' -> Some T
    | _ -> None

let parseLine : string -> Nucleotide[] = Seq.choose parse >> Seq.toArray

let readSequence() = Console.ReadLine() |> parseLine

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
    let penaltyBreak x = -1. - (1. * float x)
    let indelCost = -2.

#if DEBUG
    use _in = new IO.StreamReader("input")
    Console.SetIn(_in)
#endif
    
    Types.verbose <- 
        argv |> Array.exists ((=) "-v") ||
        argv |> Array.exists ((=) "--verbose")

    match argv |> Array.toList with
    | "Gotoh" :: _ -> 
        Gotoh.run(readSequence(), readSequence(), readSimilarity(), penaltyBreak) |> formatOutput
    | "NeedlemanWunsch" :: _ ->
        NeedlemanWunsch.run(readSequence() ,readSequence(), readSimilarity(), indelCost) |> formatOutput
    | "Hirschberg" :: _ ->
        Hirschberg.run(readSequence(), readSequence(), readSimilarity(), indelCost) |> formatOutput
    | _ -> printfn "usage: SequenceAlignment.exe [Gotoh|NeedlemanWunsch|Hirschberg] [-v|--verbose]"

    0