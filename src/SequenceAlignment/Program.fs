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

let parse' = 
    function 
    | '-' -> Break
    | 'A' -> Nucl A
    | 'C' -> Nucl C
    | 'G' -> Nucl G
    | 'T' -> Nucl T
    | x -> failwithf "parse error: %c" x 

let parseLine : string -> Nucleotide[] = Seq.choose parse >> Seq.toArray

let readSequence() = Console.ReadLine() |> parseLine

let readSequences() = 
    let line =  ref (Console.ReadLine())
    seq { 
        while !line |> String.IsNullOrEmpty |> not do 
            yield !line |> parseLine
            line := Console.ReadLine()
    } |> Seq.toArray

let readSimilarity() : Similarity' =
    let parseLine (s:string) = s.Split([|';'|]) |> Array.map float
    let lookup = 
        [Nucl A;Nucl C;Nucl G;Nucl T;Break] 
        |> List.map (fun n -> n, Console.ReadLine() |> parseLine)
        |> Map.ofList
    (fun (f,s) ->
        let f',s' = max f s, min f s
        let index = match s' with Nucl A -> 0 | Nucl C -> 1 | Nucl G -> 2 | Nucl T -> 3 | Break -> 4
        lookup.[f'].[index])

let readMultiAlignment() : MultiAlignment = 
    let line =  ref (Console.ReadLine())
    seq { 
        while !line |> String.IsNullOrEmpty |> not do 
            yield line.Value.ToCharArray() |> Array.map parse'
            line := Console.ReadLine()
    } |> array2D

[<EntryPoint>]
let main argv = 

#if DEBUG
    use _in = new IO.StreamReader("in_2")
    Console.SetIn(_in)
#endif
    
    Types.verbose <- 
        argv |> Array.exists ((=) "-v") ||
        argv |> Array.exists ((=) "--verbose")

    let sim = readSimilarity()
    Console.ReadLine() |> ignore
    match argv |> Array.toList with
    | "Gotoh" :: _ -> 
        Gotoh.run(readSequence(), readSequence(), sim) |> formatOutput
    | "NeedlemanWunsch" :: _ ->
        NeedlemanWunsch.run(readSequence() ,readSequence(), sim) |> formatOutput
    | "Hirschberg" :: _ ->
        Hirschberg.run(readSequence(), readSequence(), sim) |> formatOutput
    | "profile" :: _ -> 
        MultiAlign.profile(readMultiAlignment()) |> printfn "%A"
    | "cons" :: _ ->
        MultiAlign.consensusWord(readMultiAlignment(), sim) |> printfn "%A"
    | "malign" :: _ -> 
        MultiAlign.alignByProfiles(readMultiAlignment(), readMultiAlignment(), sim) |> printfn "%A"
    | "UPGMA" :: _ ->
        MultiAlign.UPGMA(readSequences(), sim) |> printfn "%A"

    | _ -> printfn "usage: SequenceAlignment.exe <command> [-v|--verbose]"
           printfn """
    available commands:
    - Gotoh 
    - NeedlemanWunsch 
    - Hirschberg 
    - profile 
    - cons 
    - malign 
    - UPGMA
"""

    0