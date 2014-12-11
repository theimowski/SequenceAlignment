module SequenceAlignment.Hirschberg

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

let inline lastRow (x: _[,]) = x.[Array2D.length1 x-1,*]

let NeedlemanWunschLast
    (fstSeq : Sequence, sndSeq: Sequence, sim : Similarity, indelCost : float) 
    : Alignment[] =
    let len1,len2 = fstSeq.Length, sndSeq.Length
    let prev_row = Array.init (len2+1) (fun j -> float j * indelCost)
    let cur_row = Array.copy prev_row
    for i in 1..len1 do 
        cur_row.[0] <- float i * indelCost
        for j in 1..len2 do
            cur_row.[j] <- [prev_row.[j-1] + sim(fstSeq.[i-1],sndSeq.[j-1])
                            prev_row.[j] + indelCost
                            cur_row.[j-1]+indelCost
                            ] |> List.max
        Array.blit cur_row 0 prev_row 0 (len2+1)
            
    cur_row

let inline splitBefore i (x: _[]) = 
    match i with 
    | _ when i < 1 -> [||], x
    | _ when i > x.Length-1 -> x, [||]
    | _ -> x.[0..i-1], x.[i..(x.Length-1)]

let run
    (fstSeq : Sequence, sndSeq: Sequence,sim : Similarity, indelCost : float)
    : Alignment * list<Nucleotide' * Nucleotide'> =
    let rec run' (fstSeq,sndSeq,cont) = 

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
            printfn "before"

            let upper = 
                NeedlemanWunschLast(f1,sndSeq,sim,indelCost) 
            let lower =
                NeedlemanWunschLast(Array.rev f2,
                                Array.rev sndSeq, sim, indelCost)

            printfn "%A" upper
            printfn "%A" lower

            let jmid = 
                (upper,lower |> Array.rev) 
                ||> Array.map2 (+)
                |> Array.mapi (fun i v -> i,v)
                |> Array.maxBy snd
                |> fst
            
            let s1,s2 = sndSeq |> splitBefore jmid
            
            run' (f1,s1,(fun (lA,lL) ->
                run' (f2,s2,(fun (rA,rL) ->
                    cont(lA+rA,lL@rL)))))

    run' (fstSeq,sndSeq,id)
//#time
//run([|A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A;G;T;A;C;G;C;A|],
//           [|T;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;A;T;G;C;G;C;T;A;T;G;C;G;C;T;A|]
//           ,sim,-2.) |> formatOutput
//#time