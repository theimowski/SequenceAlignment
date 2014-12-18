module SequenceAlignment.Hirschberg
    
let inline splitBefore i (x: _[]) = 
    match i with 
    | _ when i < 1 -> [||], x
    | _ when i > x.Length-1 -> x, [||]
    | _ -> x.[0..i-1], x.[i..(x.Length-1)]

let split 
    (fstSeq : Sequence, sndSeq: Sequence,sim' : Similarity') =
    let imid = fstSeq.Length / 2
    let f1,f2 = fstSeq |> splitBefore imid

    let upper = NeedlemanWunsch.runScoreLastRow(f1,sndSeq,sim') 
    let lower = NeedlemanWunsch.runScoreLastRow(Array.rev f2, Array.rev sndSeq, sim')

    let jmid = 
        (upper,lower |> Array.rev) 
        ||> Array.map2 (+)
        |> Array.mapi (fun i v -> i,v)
        |> Array.maxBy snd
        |> fst
                        
    let s1,s2 = sndSeq |> splitBefore jmid
    f1,f2,s1,s2

let printState(f,s) = 
    if verbose then
        printfn "Hirschberg running for sequences:"
        printfn "%s" <| formatSeq (f |> Seq.map Nucl)
        printfn "%s" <| formatSeq (s |> Seq.map Nucl)
        printfn ""

let run
    (f : Sequence, s: Sequence,sim' : Similarity')
    : Alignment * list<Nucleotide' * Nucleotide'> =
    let rec run' (fstSeq : Sequence,sndSeq : Sequence,cont) = 

        printState(fstSeq,sndSeq)

        match fstSeq,sndSeq with
        | [||], _ -> 
            let al = sndSeq |> Array.map (fun s -> Break, Nucl s)
            cont(al |> Array.sumBy sim', List.ofArray al)
        | _, [||] -> 
            let al = fstSeq |> Array.map (fun f -> Nucl f, Break)
            cont(al |> Array.sumBy sim', List.ofArray al)
        | [|_|], _ | _, [|_|] -> cont(NeedlemanWunsch.run(fstSeq,sndSeq,sim'))
        | _ -> 
            let f1,f2,s1,s2 = split(fstSeq,sndSeq,sim')
            
            run' (f1,s1,(fun (lA,lL) ->
                run' (f2,s2,(fun (rA,rL) ->
                    cont(lA+rA,lL@rL)))))

    run' (f,s,id)