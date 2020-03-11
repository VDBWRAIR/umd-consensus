<<<<<<< HEAD
from functools import partial
from typing import List, Tuple, Dict, NamedTuple
from pathlib import Path
import glob
import os
import jsonschema
from collections import OrderedDict
from toolz import dicttoolz as dtz
from jsonschema.exceptions import ValidationError
import argparse
import os.path

Base = str


class VCFRow(NamedTuple):
    POS: int
    DP: int
    # def get_depth(self) -> int: ...  #  def get_ref(self) -> Base: ...
    REF: Base


class VCFCall(VCFRow):
    AF: float
    INDEL: bool
    ALT: Base


class VCFNoCall(VCFRow): ...


class PercentAnalysis(NamedTuple):
    mind: int
    majority: int


Schema = dict
# MAJORITY = 80
Log = str
IUPAC_AMBIG = ['A', 'C', 'T', 'G', 'R', 'Y',
    'S', 'W', 'R', 'K', 'V', 'H', 'D', 'B']


def main():
 parser = argparse.ArgumentParser()
 parser.add_argument('-p', '--bam_input_file', type=str,
                     required=True, help='Path to input file to be read')
 parser.add_argument('-i', '--ref_input_file', type=str,
                     required=True, help='Path to input file to be read')
 parser.add_argument('-o', '--output_file', type=str,
                     required=True, help='Path to input file to be read')
 parser.add_argument('-n', '--mind', type=int, required=True,
                     help='input value for min depth')
 parser.add_argument('-r', '--vcf_input_file', type=str,
                     required=True, help='Path to input file to be refernce file')
 parser.add_argument('-m', '--majority', type=int,
                     required=True, help='input value for majoirty')
 args = parser.parse_args()
 process = consensus(args.n, args.i, args.r, args.m)
 return process

#Write process into args.o
# args.vcf_input_file, args.majority
 

if __name__ == '__main__':
    main()

Base = str
RefSeq = List[Base]
SamDepth = int

 #The consensus program is designed to take iterate through VCFrows to curate and generate consensus sequences.
 #Define Classses for the consensus program include VCFRow, VCFCall, VCFNoCall, PercentAnalysis

''' The percentN(PercentAnalysis) function takes in two arguements i.e. minimum depth(eg. 10, 500 or 1000) and mamjotiry cal i.e.(80, 95 or 99.
 percentN returns four possible outcomes outcomes as an OrderedDict based on preconditions of depth, allel frequency and percent anbnalysis.
 all outcomes are returned as a sorted orederedDict based on N, ALT, MAJORITY, and Ref.'''


def percentN(MIND, MAJORITY):  # TODO: requied
    return OrderedDict({"callN": {"DP": {"exclusiveMaximum": MIND},  "<RESULT>": lambda x: "N"},
         # TODO: ambiguous
         "callAlt": {"AF": {"inclusiveMinimum": MAJORITY}, "<RESULT>": lambda x: x['ALT']},
         "callAmbiguous": {"AF?": {"inclusiveMinimum": 100 - MAJORITY, "exclusiveMaximum": MAJORITY}, "<RESULT>": "???"},
         "callRef": {"anyOf": [{"PRC": 100 - MAJORITY},
                     {"AF": {"inclusiveMaximum": 100 - MAJORITY}},
                     {"REF": {"anyOf":  IUPAC_AMBIG}}]},  # calling ref is the default case
                     "<RESULT>": lambda x: x["REF"]})
# , "required" : ["AF", "ALT"]
# , "required" : True
# def consensus(depths: List[SamDepth], ref: RefSeq, alts: List[VCFRow], mind: int, majority: int) -> List[Base]:


''' The consensus function takes in four  arguements i.e. the depths as a List of inetegers, and the reference as a list of ref bases,
alts as a list of VCFRows and percentanalysis,It returns a list of Tuples with a log file.
Precoditions include checking that there is an equal number of depths correlating with the lenght of the reference,
the depths is not equal to zero. Mechanism: There is an iterative function with defining condtions that  loops through a range of positions
to return consensus baseses as a list of tuples'''


def consensus(depths: List[int], ref: List[Base], alts: List[VCFRow], pa: PercentAnalysis) -> Tuple[List[Base], List[Log]]:
     assert len(depths) == len(ref), f"{depths}, {ref}"
     assert len(depths) > 0
     just_alts = {row.POS: row for row in alts}
     ordered_alts = []
     with_pos = zip(range(1_000_000_000), depths, ref)
     for pos, depth, ref in with_pos:
           existing_alt = just_alts.get(pos, None)
     if existing_alt:
           ordered_alts.append(existing_alt)
     else:
           ordered_alts.append(VCFNoCall(pos, depth, ref))
           as_list_tuples = consensus_impl(ordered_alts, pa.mind, pa.majority)
           return tuple(zip(*as_list_tuples))


''' The consensus_impl takes 3 arguments i.e. (List of VCFRow, minimum depth and majority call analysis eg(80, 95 or 99).
The keys specify the four calling options based on the schema or preconditions adressed by the percentN function.
the return is of base calls the utilizes a callbase function'''


def consensus_impl(sorted_alts: List[VCFRow], mind: int, majority: int) -> List[Tuple[Base, Log]]:
    keys = ["callN", "callAlt", "callAmbiguous", "callRef"]
    schema = percentN(mind, majority)
    result = list(map(lambda alt:  call_base(schema, alt), sorted_alts))
    return result


'''The callbase function takes in two arguements i.e. VCFRow and OrderedDict(Dictionary).
It utilizes the jsonschema to validate base calls, read up on jason schema for more details'''


def call_base(schemaNodes: OrderedDict, alt: VCFRow) -> Tuple[Base, Log]:
  vcf_dict = alt.__dict__
  log = []
  for key, node in schemaNodes.items():
    try:
      justSchema = {"type": "object",  "properties": dtz.dissoc(node, '<RESULT>')}
      print(vcf_dict, justSchema)
      jsonschema.validate(vcf_dict, justSchema)
      result = node["<RESULT>"](vcf_dict)
      return result, log
    except ValidationError as e:  # TODO: handle`jsonschema.exceptions.SchemaError
      log.append((vcf_dict, e))
      log.append(ValueError(f"FAIL"))
      result = None
      return result, log




 # The work flow or call stack for these fucntions
  #consensus_iml cals on the callbase and percentN function,
 #the consensus function calls on consensus_iml.
#The main entry point fucntion is consensus and the base case fucntion are callbase and percentN
#with cosensus_impl being the middle function. Validate preconditions through assertion statements'''

def write_validate(output_file: List[VCFRow])-> bool:   
    try:
        if os.access(output_file, os.W_OK):
            return True 
    except IOError as err:
        print("could not write file {0}: {1}".format(output_file, err))
        return None

def read_validate(inputfiles )-> bool:
    try:
        for input_file in inputfiles:
          if os.access(input_file, os.R_OK):
              return True
    except IOError as err:
        print("could not read file {0}: {1}".format(input_file, err))
        return None

def inputArgs(bam_input_file, vcf_input_file, ref_input_file,out_file)-> bool:
    #Input files include 3 input files form user and a name outputile
    inputfiles = {bam_input_file, vcf_input_file, ref_input_file}
    
    output= write_validate(out_file)

    for file in inputfiles:
        done = read_validate(file)
    return done and output
    #Postcodition: returns a boolean for read_validate function


      
def bound(mind: int, majority: int):
    try:
        for mind in range(0,20) and majority in range(80,100):
            assert mind >= 20 and majority >= 80
            print("valid input")
    except IOError as err:
        print("invalid input{0}: {1}".format(mind,majority))
        return None
    
      








#  return   OrderedDict({ "callN" : {"DP" : { "exclusiveMaximum" : MIND},  "<RESULT>" : lambda x: "N" } ,
 #                  "callAlt" : { "AF" : { "inclusiveMinimum" : MAJORITY} , "<RESULT>" : lambda x: x['ALT'] },  #TODO: ambiguous
#                   "callAmbiguous" : {"AF?" : { "inclusiveMinimum" : 100 - MAJORITY, "exclusiveMaximum" : MAJORITY }, "<RESULT>" : "???" },
 #                  "callRef" : { "anyOf" : [{ "PRC" : 100 - MAJORITY },
  #                             { "AF" : { "inclusiveMaximum" : 100 - MAJORITY } }] },
   #                            "<RESULT>" : lambda x: x["REF"] })

=======
from vcf.model import _Record
import vcf
from typing import Dict,Optional,Sequence,Callable,NamedTuple,Tuple,Iterable,List, Union, TypeVar
# from typing_extenions import T
from toolz.itertoolz import first, second
from toolz.dicttoolz import valfilter
from typing import NamedTuple
import sys
#from plumbum.cmd import samtools
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from itertools import groupby
from functools import partial
from Bio import SeqIO
import argparse
import subprocess
DNA = 'ACTGactg'
T = TypeVar('T')

class PercentAnalysis:
    def __init__(self, mind: int, majority: int, minbq: int = 0, sample: str = None): 
        self.mind, self.majority, self.minbq, self.sample = mind, majority, minbq, sample

Result = Tuple[Sequence[str], Sequence[str]]
Base = str

def samtoolsDepth(ref: SeqRecord, bam: str, minbq: int):
    output = subprocess.check_output(['samtools', 'depth', bam, '-r', ref.id, '-q', str(minbq)], universal_newlines=True)
    #lines = samtools['depth'][bam, '-r', ref.id, '-q', minbq]().split('\n')#lines = open('depth.txt')
    lines = filter(lambda x: x.strip(), output.split('\n'))
    lines = map(lambda x: x.split('\t'), lines)
    byPos = { int(x[1]) - 1 :  int(x[2])  for x in lines }
    depths = []
    for i in range(0, len(ref.seq)):
        depths.append(byPos.get(i, 0))
    return depths
        
def group_muts_by_refs(references, muts): # type: (List[SeqRecord], List[VCFRow]) -> List[List[VCFRow]]
    '''group and sort the mutations so that they match the order of the references.'''
    #NOTE: muts will already be "sorted" in that they are grouped together in the vcf
    #fix the groupby so it doesn't incidentally drain the first object of the group
    unzip = lambda x: zip(*x)
    chroms, groups = unzip(map(lambda kv: (kv[0], list(kv[1])), groupby(muts, lambda x: x.CHROM)))
    def index_of_ref(chrom_first): # type: (Tuple[str, List[SeqRecord]]) -> int
        return  list(map(lambda x: x.id, references)).index(chrom_first[0])
    _, muts_by_ref = unzip(sorted(zip(chroms, groups), key=index_of_ref))
    return muts_by_ref

def io_runner(ref_fasta: str, vcf_path: str, bam_path: str, out_path: str, log_path: str, pa: PercentAnalysis) -> None:
  _refs = list(SeqIO.parse(ref_fasta, 'fasta'))
  with open(vcf_path, 'r') as vcf_handle:
    recs =  vcf.Reader(vcf_handle)
    muts_by_ref = group_muts_by_refs(_refs, recs)
    single_dispatch = partial(consensus_wrapper, pa, bam_path)
    results = list(map(single_dispatch, _refs, muts_by_ref))
  with open(out_path, 'w') as outf:
    for seqr in results:
      SeqIO.write(outf, seqr, 'fasta')


def consensus_wrapper(pa: PercentAnalysis, bam_path: str, ref: SeqRecord, vars: Sequence[_Record], sample=None) -> SeqRecord:
    depths = samtoolsDepth(ref, bam_path, pa.minbq)
    result_seq = consensus(depths, list(ref.seq), vars, pa)
    id_str = f">{ref.id}_{sample or ''}:Consensus"
    return SeqRecord(seq=Seq(''.join(result_seq)), id=id_str)


        
def consensus_str(ref: str, consensus: str, sample: str): # type: (SeqRecord, str) -> str
    return ">{0}_{1}:Consensus\n{2}".format(sample if sample else '', ref.id, consensus)
    # return SeqRecord(seq=Seq(''.join(new_ref)), id=ref_seq.id)

def consensus(depths: Sequence[int], ref: List[str], alts: Sequence[_Record], pa: PercentAnalysis) -> Result :
    result= simplest(pa.mind, pa.majority/100.0, ref, depths, alts)
    print(result)
    return result

def simplest(mind: int, majority: float, ref: List[str], depths: Sequence[int], vars: Sequence[_Record]) -> Result:
    assert 0 < majority < 1
    withn_consensus: List[str] = ref[:]
    DEFAULT_UNDER_DEPTH = DEFAULT_IS_INDEL = DEFAULT_NO_CALL = 'N'
    under_depths, msgs = set(), []
    for i, dp in enumerate(depths):
        if dp < mind:
            under_depths.add(i)
            withn_consensus[i] = DEFAULT_UNDER_DEPTH
    msgs.append(f"The following positions were replaced with `{DEFAULT_UNDER_DEPTH}` becauase they were under minimum depth {mind} according to the bam file: {list(under_depths)}")
    out_consensus = withn_consensus[:]
    must_be_list: Callable[[Union[List[T],T]],List[T]] = lambda xs: xs if hasattr(xs, '__iter__') else [xs]
    base: str
    for var in vars:
        var.__dict__['af']  =  must_be_list(var.INFO['AF'])
        var.__dict__['alt'] = must_be_list(var.ALT) #  var.af, var.alt = must_be_list(var.AF), must_be_list(var.ALT)
    for var in vars:
        pos = var.POS - 1
        assert var.REF == ref[pos]
        consider_var = any(map(lambda af: af > (1-majority), var.af)) and not pos in under_depths
        if not consider_var:
            continue
        if var.is_indel(): # includes deletion
            base, msg = DEFAULT_IS_INDEL, f"called POS={var.POS} the default base {DEFAULT_IS_INDEL} because it was called an indel by pyvcf: {var}."
        alt_base = try_call_alt(var, majority)
        if alt_base:
            base, msg = alt_base, f"called POS={var.POS} alternate because base {alt_base} had AF >= {majority}"
        else:
            ambig_base = try_call_ambiguous(majority, withn_consensus[pos], var)
            if ambig_base:
                base, msg = ambig_base, f"called POS={var.POS} 'ambiguous' base {ambig_base} because there was no ALT with AF >= {majority} but at least one ALT was >= {1-majority}"
            else:
                base, msg = DEFAULT_NO_CALL, f"called POS={var.POS} the default base {DEFAULT_NO_CALL} because calling it as an alternate and calling it as an ambiguous base failed, and the depth was still over minimum {mind}"
        # base should only have length greater than 1 if it was an alternate call.
        out_consensus[pos:pos+len(base)] = base
        msgs.append(msg)
        # for s in msgs: sys.stderr.write(s)
    return out_consensus, msgs

def try_call_alt(var: _Record, majority: float) -> Optional[str]:
   called = [ alt for af, alt in zip(var.af, var.alt) 
                     if af > majority and alt in DNA]
   return None if not called else called[0]

def try_call_ambiguous(majority: float, ref_base: Base, var: _Record) -> Optional[str]:
    alt_over_ambig_thresh = [ alt for af, alt in zip(var.af, var.alt) 
                 if af > (1 - majority) ]
    only_first_alt_bases = [ x[0] for x in alt_over_ambig_thresh ]
    prc = 1-sum(var.af) 
    ref_contrib = [ref_base] if prc >= 1-majority else []
    # we ignore the ref base if it is ambiguous, and we ignore the alt base if it is weird (like '*')
    # this can result in an "ambiguous call" giving a non-ambiguous base (i.e. any of A,T,C,G)
    valid_ambig_combo = list(filter(DNA.__contains__, only_first_alt_bases + ref_contrib))
    return make_ambig(valid_ambig_combo)

def make_ambig(bs_: Sequence[str]) -> Optional[str]:
    bs = list(map(str.upper, bs_))
    AMBIGUITY_TABLE = { 'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B', 'ACGT': 'N' }
    assert all(map(AMBIGUITY_TABLE.keys().__contains__, bs))
    cleaned = ''.join(sorted(bs))
    return AMBIGUITY_TABLE.get(cleaned)



# as opposed to using trim_ref, just fill with `N`s, and if required add 
# feature of turning end-facing `N`s into '-'

#def uncoveredPositions(mind, minbq, bam, ref):
#    depthStats = samtoolsDepth(str(ref.id), bam, minbq) # use ref string
#    allPositions = range(1, len(ref.seq)+1)
#    underStats = filter(lambda x: x['depth'] < mind, depthStats)
#    underPositions = map(lambda x: x['pos'], underStats)
#    allDepthPositions = map(lambda x: x['pos'], depthStats)
#    zeroPositions = set(allPositions) - set(allDepthPositions)
#    return sorted( map( lambda x: x - 1, set(underPositions) | zeroPositions ) )
#
>>>>>>> origin/fixes
