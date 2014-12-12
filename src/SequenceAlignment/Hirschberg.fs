module SequenceAlignment.Hirschberg
    
let inline splitBefore i (x: _[]) = 
    match i with 
    | _ when i < 1 -> [||], x
    | _ when i > x.Length-1 -> x, [||]
    | _ -> x.[0..i-1], x.[i..(x.Length-1)]

let split 
    (fstSeq : Sequence, sndSeq: Sequence,sim : Similarity, indelCost : float) =
    let imid = fstSeq.Length / 2
    let f1,f2 = fstSeq |> splitBefore imid

    let upper = NeedlemanWunsch.runScoreLastRow(f1,sndSeq,sim,indelCost) 
    let lower = NeedlemanWunsch.runScoreLastRow(Array.rev f2, Array.rev sndSeq, sim, indelCost)

    let jmid = 
        (upper,lower |> Array.rev) 
        ||> Array.map2 (+)
        |> Array.mapi (fun i v -> i,v)
        |> Array.maxBy snd
        |> fst
                        
    let s1,s2 = sndSeq |> splitBefore jmid
    f1,f2,s1,s2

let run
    (f : Sequence, s: Sequence,sim : Similarity, indelCost : float)
    : Alignment * list<Nucleotide' * Nucleotide'> =
    let rec run' (fstSeq,sndSeq,cont) = 

        match fstSeq,sndSeq with
        | [|_|], _ | _, [|_|] -> cont(NeedlemanWunsch.run(fstSeq,sndSeq,sim,indelCost))
        | [||], _ -> 
            cont(indelCost * float sndSeq.Length, 
                 sndSeq |> Array.map (fun s -> Break, Nucl s) |> List.ofArray)
        | _, [||] -> 
            cont(indelCost * float fstSeq.Length, 
                 fstSeq |> Array.map (fun f -> Nucl f, Break) |> List.ofArray)
        | _ -> 
            let f1,f2,s1,s2 = split(fstSeq,sndSeq,sim,indelCost)
            
            run' (f1,s1,(fun (lA,lL) ->
                run' (f2,s2,(fun (rA,rL) ->
                    cont(lA+rA,lL@rL)))))

    run' (f,s,id)