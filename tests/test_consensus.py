from typing import List
from consensus.consensus import Base, PercentAnalysis
from consensus.consensus import consensus
from dataclasses import dataclass
import unittest
RefSeq = List[Base]
Base = str
depths = ()
dumbRows = ()


tc = unittest.TestCase()
SamDepth= int

@dataclass
class VCFCall:
    AF: float
    DP: int
    INDEL: bool
    ALT: Base
    POS: int
    REF: Base




def test_simplest_consensus20():
    refSeq = ["A","M","Y","C", "T"]
    depths = [1000,5, 150,100,250]
    #dumbRows = [VCFCall(AF=1, DP=1, INDEL=False, ALT='X', POS=i+1, REF='M') for i in range(3)]
    alts = [VCFCall(AF=50, DP=5, INDEL=False, ALT='A', POS=1, REF='M')] #+dumbRows
    expected = ["A","N","Y","C","T"]
    result, log =  consensus(depths=[100, 100, 100, 100, 100], ref = refSeq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    assert expected == list(result)
    print(f"consensus({depths}, {refSeq}, {alts})")

#test_simplest_consensus20(i)



def test_filter_ambigs():
    refseq = ["G","T","S","T"]
    alts= [VCFCall(AF=75, DP=500 , INDEL=False, ALT="C", POS=2, REF="S" )]
    expected = ["G","T","S","T"]
    result, log =  consensus(depths=[100, 100, 100, 100], ref = refseq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    assert expected == list(result)
    print(f"consensus({depths}, {refseq}, {alts})")


def test_AF():
    refseq = ["G","T","S","T"]
    alts= [VCFCall(AF=15, DP=500 , INDEL=False, ALT="C", POS=2, REF = "S" )]
    expected = ["G","T","G","T"]
    result, log =  consensus(depths=[100, 100, 100, 100], ref = refseq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    assert expected == list(result)
    print(f"consensus({depths}, {refseq}, {alts})")


def test_variant():
    refseq = ["G","T","S","T"]
    alts= [VCFCall(AF=85, DP=500 , INDEL=False, ALT="A", POS=2, REF = "S" )]
    expected = ["G","T","A","T"]
    result, log =  consensus(depths=[100, 100, 100, 100], ref = refseq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    assert expected == list(result)
    print(f"consensus({depths}, {refseq}, {alts})")


def test_variant_ambig():
    refseq = ["G","T","S","T"]
    alts= [VCFCall(AF=50, DP=500 , INDEL=False, ALT="A", POS=2, REF = "S" )]
    expected = ["G","T","R","T"]
    result, log =  consensus(depths=[100, 100, 100, 100], ref = refseq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    print(f"log: {log}")
    assert expected == list(result), f"{result} != {expected}"
   # print(f"consensus({depths}, {refseq}, {alts})")


def test_simplest_consensus20():
  refSeq = [ 'A', 'C', 'C', 'C' ]
  depths =  [ 100 ] * 4
  dumbRows = [VCFCall(AF=1, DP=1, INDEL=False, ALT='X', POS=i+1, REF='C') for i in range(3)]
  alts = [VCFCall(AF=81, DP=100, INDEL=False, ALT='G', POS=0, REF='A')] + dumbRows
  actual, log = consensus(depths, refSeq, alts, pa=PercentAnalysis(10, 80))

  # actual, log = zip(*actual) # [x[0] for x in actual]
  print(f"log: {log}")
  expected = [ 'G', 'N', 'N', 'N' ]
  assert  expected == list(actual) , f"{actual} != {expected}"





