from dataclasses import dataclass
from typing import List,Tuple,Dict
import jsonschema
from collections import OrderedDict
from toolz import dicttoolz as dtz
from jsonschema.exceptions import ValidationError
Base = str
RefSeq = List[Base]
SamDepth = int
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

def consensus_impl(sorted_alts: List[VCFRow], mind: int, majority: int) -> List[Tuple[Base, Log]]:
  keys = ["callN", "callAlt", "callAmbiguous", "callRef"]
  schema = percentN(mind, majority)
  result = list( map(lambda alt:  call_base(schema, alt), sorted_alts) )
  return result

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



#  return   OrderedDict({ "callN" : {"DP" : { "exclusiveMaximum" : MIND},  "<RESULT>" : lambda x: "N" } ,
 #                  "callAlt" : { "AF" : { "inclusiveMinimum" : MAJORITY} , "<RESULT>" : lambda x: x['ALT'] },  #TODO: ambiguous
#                   "callAmbiguous" : {"AF?" : { "inclusiveMinimum" : 100 - MAJORITY, "exclusiveMaximum" : MAJORITY }, "<RESULT>" : "???" },
 #                  "callRef" : { "anyOf" : [{ "PRC" : 100 - MAJORITY },
  #                             { "AF" : { "inclusiveMaximum" : 100 - MAJORITY } }] },
   #                            "<RESULT>" : lambda x: x["REF"] })
