from dataclasses import dataclass
from typing import List,Tuple
import jsonschema
from collections import OrderedDict
from toolz import dicttoolz as dtz
from jsonschema.exceptions import ValidationError
Base = str
RefSeq = List[Base]
SamDepth = int
@dataclass
class VCFRow:
  AF: float
  DP: int
  INDEL: bool
  ALT: Base
  POS: int
  REF: Base
SamDepth = int
Schema = dict
# MAJORITY = 80 
Log = str 
from functools import partial
def percentN(MIND, MAJORITY): #TODO: requied
  return   OrderedDict({ "callN" : {"DP" : { "exclusiveMaximum" : MIND},  "<RESULT>" : lambda x: "N" } ,
                   "callAlt" : { "AF" : { "inclusiveMinimum" : MAJORITY} , "<RESULT>" : lambda x: x['ALT'] },  #TODO: ambiguous
                   "callAmbiguous" : {"AF?" : { "inclusiveMinimum" : 100 - MAJORITY, "exclusiveMaximum" : MAJORITY }, "<RESULT>" : "???" }, 
                   "callRef" : { "anyOf" : [{ "PRC" : 100 - MAJORITY }, 
                               { "AF" : { "inclusiveMaximum" : 100 - MAJORITY } }] },
                               "<RESULT>" : lambda x: x["REF"] })
def consensus(depths: List[SamDepth], ref: RefSeq, alts: List[VCFRow], mind: int, majority: int) -> List[Base]:
  keys = ["callN", "callAlt", "callAmbiguous", "callRef"]
  schema = percentN(mind, majority)
  result =  list(map(partial(call_base, schema), depths, ref, alts))
  return result

def call_base(schemaNodes: OrderedDict, depth: SamDepth, refBase: Base, alt: VCFRow) -> Tuple[Base, Log]:
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
