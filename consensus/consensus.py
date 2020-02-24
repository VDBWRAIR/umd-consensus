from dataclasses import dataclass
from typing import List,Tuple,Dict
import jsonschema
from collections import OrderedDict
from toolz import dicttoolz as dtz
from jsonschema.exceptions import ValidationError
import argparse
import os.path

@dataclass
class VCFRow:
POS: int
DP: int
REF: Base #  def get_depth(self) -> int: ...  #  def get_ref(self) -> Base: ...

class VCFCall(VCFRow):
AF: float
INDEL: bool
ALT: Base


class VCFNoCall(VCFRow): ...

@dataclass
class PercentAnalysis:
mind: int
majority: int

Schema = dict
# MAJORITY = 80
Log = str
from functools import partial
IUPAC_AMBIG = ['A', 'C', 'T', 'G', 'R', 'Y', 'S', 'W', 'R', 'K', 'V', 'H', 'D', 'B']


def ():
 parser = argparse.ArgumentParser()

 parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to input file to be read')
 parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to input file to be read')
 parser.add_argument('-n', '--min_depth', type=int, required=True, help='input value for min depth')
 parser.add_argument('-r', '--input_file', type=str, required=True, help='Path to input file to be refernce file')
 parser.add_argument('-m', '--majority', type=int, required=True, help='input value for majoirty')
 Base = str
 RefSeq = List[Base]
 SamDepth = int

  '''The consensus program is designed to take iterate through VCFrows to curate and generate consensus sequences.
  Define Classses for the consensus program include VCFRow, VCFCall, VCFNoCall, PercentAnalysis'''




''' The percentN(PercentAnalysis) function takes in two arguements i.e. minimum depth(eg. 10, 500 or 1000) and mamjotiry cal i.e.(80, 95 or 99.
percentN returns four possible outcomes outcomes as an OrderedDict based on preconditions of depth, allel frequency and percent anbnalysis.
all outcomes are returned as a sorted orederedDict based on N, ALT, MAJORITY, and Ref.'''
def percentN(MIND, MAJORITY): #TODO: requied
    return   OrderedDict({ "callN" : {"DP" : { "exclusiveMaximum" : MIND},  "<RESULT>" : lambda x: "N" } ,
         "callAlt" : { "AF" : { "inclusiveMinimum" : MAJORITY} , "<RESULT>" : lambda x: x['ALT'] },  #TODO: ambiguous
         "callAmbiguous" : {"AF?" : { "inclusiveMinimum" : 100 - MAJORITY, "exclusiveMaximum" : MAJORITY }, "<RESULT>" : "???"},
         "callRef" : { "anyOf" : [{ "PRC" : 100 - MAJORITY },
                     { "AF" : { "inclusiveMaximum" : 100 - MAJORITY } },
                     { "REF" : {"anyOf" :  IUPAC_AMBIG  }}] }, # calling ref is the default case
                     "<RESULT>" : lambda x: x["REF"] })
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
     just_alts = { row.POS : row for row in alts }
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
    result = list( map(lambda alt:  call_base(schema, alt), sorted_alts) )
    return result


'''The callbase function takes in two arguements i.e. VCFRow and OrderedDict(Dictionary).
It utilizes the jsonschema to validate base calls, read up on jason schema for more details'''
def call_base(schemaNodes: OrderedDict, alt: VCFRow) -> Tuple[Base, Log]:
    vcf_dict = alt.__dict__
    log = []
    for key, node in schemaNodes.items():
       try:
          justSchema = {"type" : "object", "properties" : dtz.dissoc(node, '<RESULT>') }
          print( vcf_dict, justSchema)
          jsonschema.validate( vcf_dict, justSchema)
          result = node["<RESULT>"](vcf_dict)
          return result, log
          except ValidationError as e:
          log.append((vcf_dict, e)) # TODO: handle`jsonschema.exceptions.SchemaError
          log.append( ValueError(f"FAIL"))
          result = None
          return result, log


''' The work flow or call stack for these fucntions:
consensus_iml cals on the callbase and percentN function,
the consensus function calls on consensus_iml.
The main entry point fucntion is consensus and the base case fucntion are callbase and percentN
with cosensus_impl being the middle function'''


def validate(sorted_alts: List[VCFRow], mind: int, majority:int, ref: list[bases]):

    '''  Validate preconditions through assertion statements '''
    assert ref == list[bases]
    assert sorted_alts == List[VCFRow]
    assert mind == int
    assert majoirty == int
    return null




#  return   OrderedDict({ "callN" : {"DP" : { "exclusiveMaximum" : MIND},  "<RESULT>" : lambda x: "N" } ,
 #                  "callAlt" : { "AF" : { "inclusiveMinimum" : MAJORITY} , "<RESULT>" : lambda x: x['ALT'] },  #TODO: ambiguous
#                   "callAmbiguous" : {"AF?" : { "inclusiveMinimum" : 100 - MAJORITY, "exclusiveMaximum" : MAJORITY }, "<RESULT>" : "???" },
 #                  "callRef" : { "anyOf" : [{ "PRC" : 100 - MAJORITY },
  #                             { "AF" : { "inclusiveMaximum" : 100 - MAJORITY } }] },
   #                            "<RESULT>" : lambda x: x["REF"] })
