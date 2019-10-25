open System
open System.Collections.Generic
open Argu

//////////////////////////////////////////////////////////////////////////////////
// Multisets, to represent complexes
//////////////////////////////////////////////////////////////////////////////////
type entry<'a> = {
  element:'a;
  multiplicity:int
}

type Mset<'a> = entry<'a> list
type 'a t = entry<'a> list // List is assumed to be without duplicates

let empty = [] : 'a t


let emptyMap () = new Dictionary<_,_>(HashIdentity.Structural)
let add (ht:Dictionary<'a,'b>) (key:'a) (value:'b) = ht.[key] <- value

let tryFind (ht:Dictionary<'a,'b>) (key:'a) =
  if ht.ContainsKey(key) then
    Some(ht.[key])
  else
    None

let rec rev_map_acc f l acc =
    match l with 
    | [] -> acc
    | h::t -> rev_map_acc f t (f h :: acc)

let rev_map f l = rev_map_acc f l []

let from_list xs =
  let counts = emptyMap ()
  let count_and_remove l element =
    match tryFind counts element with
      | None   -> add counts element 1
                  element::l
      | Some n -> add counts element (n+1)
                  l
  let uniques = List.fold count_and_remove [] xs
  rev_map (fun element -> {element=element; multiplicity=counts.[element]}) uniques 

let fold_left_m f acc (mset: 'a t) = List.fold f acc mset

(* returns the elements with multiplicity *)
let to_list (mset: 'a t) = 
  let rec add_duplicates l (x: 'a) = function
    | 0 -> l
    | n -> add_duplicates (x::l) x (n-1) in
  let add_mult l (e:entry<_>) = add_duplicates l e.element e.multiplicity in
  let l = fold_left_m add_mult [] mset in
  List.rev l

let Mset_add (eq:'a -> 'a -> bool) (xs:'a entry list) (xm:'a entry) : entry<'a> list =
  let rec m_add new_entry = function
    | [] -> [new_entry]
    | entry::l -> if eq entry.element new_entry.element then
                    {element=entry.element; multiplicity=entry.multiplicity + new_entry.multiplicity}::l
                  else
                    entry::m_add new_entry l
  m_add xm xs

let string_of_list (f:'a -> string) (sep:string) (xs:'a list) : string =
  let rec inner acc xs =
    match xs with
    | [] -> acc
    | [x] -> acc + (f x)
    | (x::xs) -> inner (acc + (f x) + sep) xs
  inner "" xs

let to_string e_to_string sep mset =
  let sn_to_string = function
    | {multiplicity=0} -> ""
    | {multiplicity=1; element=s} -> e_to_string s
    | {multiplicity=n; element=s} -> (string n) (*^ "*"*) + (e_to_string s) in  
  string_of_list sn_to_string sep mset

let Mset_map f mset = List.map (fun entry -> {element=f entry.element; multiplicity=entry.multiplicity}) mset
let mapm (f: entry<'a> -> entry<'b>) mset = List.map f mset

let msetNaught              = []
let msetHomomer sp          = [{element = sp;  multiplicity = 1}]
let msetHomodimer sp        = [{element = sp;  multiplicity = 2}]
let msetHeterodimer sp1 sp2 = [{element = sp1; multiplicity = 1};
                               {element = sp2; multiplicity = 1}]
//////////////////////////////////////////////////////////////////////////////////

let alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ" |>Array.ofSeq |>Array.map string

let rec pairs l = 
  [ match l with
    | h::t -> 
        yield! t |> Seq.map (fun elem -> (h, elem))
        yield! t |> pairs
    | _ -> () ]

let renameCRN (rxns : (Mset<string> * Mset<string>) list) (sub : Map<string, string>) : (Mset<string> * Mset<string>) list = 
  let f x = sub.[x]
  let g = Mset_map f >> List.sort
  rxns
  |> List.map (fun (r, p) -> g r, g p)
  |> List.sort
    
(** generate n-species permutations **)
let rec distribute e = function
  | [] -> [[e]]
  | x::xs' as xs -> (e::xs)::[for xs in distribute e xs' -> x::xs]

let rec do_permute = function
  | [] -> [[]]
  | e::xs -> List.collect (distribute e) (do_permute xs)

let permute n = alphabet.[0..n-1] 
                |> List.ofArray 
                |> do_permute
      
let makeSubs n = permute n
                  |> List.map (List.mapi (fun i x -> alphabet.[i], x))
                  |> List.map Map.ofList

let rec makeClassPermutations acc = function
| [] -> [[]]
| x :: xs -> [acc .. acc + x - 1]
              |> do_permute 
              |> List.map (List.mapi (fun i x -> alphabet.[acc+i], alphabet.[x]))
              |> fun p -> let p'= makeClassPermutations (acc + x) xs
                          List.allPairs p p'
                          |> List.map (fun (x,y) -> x @ y)
                   
let makeClassSubs classes = makeClassPermutations 0 classes |> List.map Map.ofList

// Printing methods
let print (r, p) = sprintf "%s->%s" (to_string id "+" r |> fun x -> if x = "" then "0" else x) 
                                    (to_string id "+" p |> fun x -> if x = "" then "0" else x)
let printCRN crn = crn 
                    |> List.map print
                    |> String.concat " | "

(** load a bimolecular CRN file in LBS format **)
// We use F#'s active patterns, from: https://stackoverflow.com/questions/3722591/pattern-matching-on-the-beginning-of-a-string-in-f
let (|Prefix|_|) (p:string) (s:string) =
  if s.StartsWith(p) then
      Some(s.Substring(p.Length))
  else
      None

let loadCRNFile filePath =
  let text = System.IO.File.ReadAllLines filePath
  let parseComplex (x:string) =
    match x with 
    | ""                   -> msetNaught, ""
    | Prefix "->" complex2 -> msetNaught, complex2
    | Prefix "2" sp1AndRest -> 
      let sp1 = sp1AndRest.[0] |> Char.ToString 
      match sp1AndRest.Substring(1) with 
      | ""                   -> msetHomodimer sp1, ""
      | Prefix "->" complex2 -> msetHomodimer sp1, complex2
      | _ -> failwithf "Parsing failed at: %s" x
    | _               -> 
      let sp1     = x.[0] |> Char.ToString
      match x.Substring(1) with 
      | ""                   -> msetHomomer sp1, ""
      | Prefix "->" complex2 -> msetHomomer sp1, complex2
      | Prefix "+"  sp2AndRest -> 
        let sp2 = sp2AndRest.[0] |> Char.ToString
        match sp2AndRest.Substring(1) with 
        | ""                   -> msetHeterodimer sp1 sp2, ""
        | Prefix "->" complex2 -> msetHeterodimer sp1 sp2, complex2
        | _ -> failwithf "Parsing failed at: %s" x
      | _   -> failwithf "Parsing failed at: %s" x
  let parseLine (x:string) =
    match x with 
    | Prefix "--" _ -> None
    | Prefix "|" reaction -> 
      let c1, rest = parseComplex reaction
      let c2, _    = parseComplex rest
      Some (c1, c2)
    | Prefix "CRNs generated" _ -> None
    | Prefix "Time elapsed" _ -> None
    | _   -> failwithf "Parsing failed at: %s" x
  text 
  |> Array.fold (fun (reactionsAcc, crnsAcc) line ->
        match parseLine line with 
        | None   -> ([], if not reactionsAcc.IsEmpty 
                          then reactionsAcc :: crnsAcc
                          else crnsAcc)
        | Some r -> (r::reactionsAcc, crnsAcc)) 
      ([], [])
  |> snd 
  |> List.rev


/// The test method
let testClass classes expectedPath actualPath =
  classes |> List.map (sprintf "%d") |> String.concat ";" |> printfn "Testing species classes [%s]" 
  let maxSpecies = List.sum classes
  if maxSpecies > 26 then failwith "Unsupported total number of species (max 26 species)"     

  (** load test data**)
  // compute expected CRNs with diffusibles from non-isomorphic CRNs 
  let subs = makeSubs (List.sum classes)
  let classSubs = makeClassSubs classes
  printfn "- Generating expected CRNs by permuting species"
  let loadedCRNs = loadCRNFile expectedPath 
  let expectedCRNs = 
    loadedCRNs 
    |> List.map List.sort
    |> List.collect (fun crn ->
        // expand CRN into all its permutations
        subs 
        |> List.map (renameCRN crn)
            
        // filter by class permutations
        |> List.fold (fun acc crn -> 
              if classSubs 
                  |> List.exists (fun perm -> 
                      renameCRN crn perm 
                      |> fun x -> List.contains x acc)
                  |> not
                then crn :: acc
                else acc) []
        |> List.distinct
        |> List.rev)
    |> Set.ofList
      
  printfn "- Generating actual CRNs from file"
  let actual     = loadCRNFile actualPath
                      |> List.map List.sort
  let actualCRNs = actual |> Set.ofList
  // remove identical CRNs first
  
  printf "- Computing the set difference (expected \ actual): "
  let missingExpected = Set.difference expectedCRNs actualCRNs
  printfn "found %d" missingExpected.Count
  printf "- Computing the set difference (actual \ expected): "
  let missingActual   = Set.difference actualCRNs expectedCRNs 
  printfn "found %d" missingActual.Count
  
  let finalExpected = 
    if missingExpected.Count > 0
    then 
      printfn "- Check if the missing CRNs are actually present as an isomorph"
      missingExpected
      |> Set.toList 
      |> List.filter (fun crn -> 
        classSubs 
        |> List.map (renameCRN crn)
        |> List.exists (fun isomorph -> actualCRNs |> Set.contains isomorph)
        |> not)
    else []
  
  let finalActual          = 
    if missingActual.Count > 0
    then 
      missingActual
      |> Set.toList 
      |> List.filter (fun crn -> 
        classSubs 
        |> List.map (renameCRN crn)
        |> List.exists (fun isomorph -> expectedCRNs |> Set.contains isomorph)
        |> not)
    else []

  if not (finalExpected.IsEmpty && finalActual.IsEmpty)
  then 
    // prepare error msg
    let provenance, counterExample = 
      if not finalExpected.IsEmpty 
      then "actual",   finalExpected
      else "expected", finalActual
      |> fun (x, y) ->  x, y |> List.head |> printCRN
    printfn "The following CRN was not found in the %s list: \n%s\n" provenance counterExample
    1
  else 
    printfn "Test passed!"
    0

// Parser
type Arguments =
  | [<Mandatory;AltCommandLine("-s")>] Species_Classes of int list
  | [<Mandatory;AltCommandLine("-e")>] Expected of string
  | [<Mandatory;AltCommandLine("-a")>] Actual of string
  | Test_Directory of path:string
  with
  interface IArgParserTemplate with
      member s.Usage =
          match s with
          | Species_Classes _ -> "Species classes to test"
          | Expected _        -> "File path for genCRN results for no species classes"
          | Actual _          -> "File path for genCRN results with species classes"
          | Test_Directory _  -> "Directory where genCRN results can be read from"

[<EntryPoint>]
let main args = 
    let parser = ArgumentParser.Create<Arguments>(programName = "TestGenCRN.exe")
    let parse_results = parser.Parse args

    let species_classes = parse_results.GetResult Species_Classes
    let test_dir = parse_results.GetResult (Test_Directory, defaultValue = "")
    let path_expected = System.IO.Path.Combine (test_dir, parse_results.GetResult Expected)
    let path_actual   = System.IO.Path.Combine (test_dir, parse_results.GetResult Actual)

    testClass species_classes path_expected path_actual
