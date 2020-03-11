from typing import List, NamedTuple, Optional
from consensus.impl import consensus, io_runner, Base, PercentAnalysis

import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class VCFCall:
    def __init__(self, AF: float ,DP: int ,INDEL: bool ,ALT: Base ,POS: int ,REF: Base ,CHROM: Optional[str] = None):
        self.DP, self.INDEL, self.ALT, self.POS, self.REF, self.CHROM = DP, INDEL, ALT, POS, REF, CHROM
        self.INFO = { 'AF' : AF}
    def is_indel(self) -> bool:
       return self.INDEL


# def VCFCall(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample_indexes
# NOTE: undefined behavior when ALT matches REF
def test_calls_ambig() -> None:
    refSeq = ["A","M","Y","C", "T"]
    depths = [1000,5, 150,100,250]
    alts = [VCFCall(AF=.50, DP=100, INDEL=False, ALT='G', POS=1, REF='A')] #+dumbRows
    expected = ["R","N","Y","C","T"]
    result, log =  consensus(depths=depths, ref = refSeq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    assert expected == list(result)
    print(f"consensus({depths}, {refSeq}, {alts})")

def test_calls_N_no_alt() -> None:
    refSeq = ["A","M","Y","C", "T"]
    depths = [1000,5, 150,100,250]
    alts = [VCFCall(AF=.05, DP=5, INDEL=False, ALT='G', POS=2, REF='M')] #+dumbRows
    expected = ["A","N","Y","C","T"]
    result, log =  consensus(depths=depths, ref = refSeq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    assert expected == list(result)
    print(f"consensus({depths}, {refSeq}, {alts})")

def test_simplest_consensus20() -> None:
    refSeq = ["A","M","Y","C", "T"]
    depths = [1000,5, 150,100,250]
    alts = [VCFCall(AF=.50, DP=5, INDEL=False, ALT='G', POS=2, REF='M')] #+dumbRows
    expected = ["A","N","Y","C","T"]
    result, log =  consensus(depths=depths, ref = refSeq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    assert expected == list(result)
    print(f"consensus({depths}, {refSeq}, {alts})")
#test_simplest_consensus20(i)


def test_under_majority_single_ambig() -> None:
    refseq = ["G","T","S","T"]
    alts= [VCFCall(AF=75, DP=500 , INDEL=False, ALT="C", POS=2, REF="T" )] # , VCFCall(AF=25, DP=500 , INDEL=False, ALT="T", POS=2, REF="S" )]
    expected = ["G","C","S","T"]
    result, log =  consensus(depths=[100, 100, 100, 100], ref = refseq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    assert expected == list(result)


def test_keep_ref_under_one_minus_majority() -> None:
    refseq = ["G","T","A","T"]
    alts= [VCFCall(AF=.15, DP=500 , INDEL=False, ALT="C", POS=3, REF = "A" ), VCFCall(AF=.15, DP=500 , INDEL=False, ALT="C", POS=3, REF = "A" )]
    expected = ["G","T","A","T"]
    depths = [100, 100, 100, 100]
    result, log =  consensus(depths=depths, ref = refseq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    assert expected == list(result)
    print(f"consensus({depths}, {refseq}, {alts})")


def test_variant() -> None:
    refseq = ["G","T","S","T"]
    alts= [VCFCall(AF=85, DP=500 , INDEL=False, ALT="A", POS=3, REF = "S" )]
    expected = ["G","T","A","T"]
    depths = [100, 100, 100, 100]
    result, log =  consensus(depths=depths, ref = refseq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    assert expected == list(result)
    print(f"consensus({depths}, {refseq}, {alts})")


def test_variant_ambig_ref_not_considered_in_ambig_call() -> None:
    refseq = ["G","T","S","T"]
    alts= [VCFCall(AF=.50, DP=500 , INDEL=False, ALT="A", POS=3, REF = "S" )]
    expected = ["G","T","A","T"]
    result, log =  consensus(depths=[100, 100, 100, 100], ref = refseq, alts = alts, pa=PercentAnalysis(mind = 10, majority = 80))
    print(f"log: {log}")
    assert expected == list(result), f"{result} != {expected}"
   # print(f"consensus({depths}, {refseq}, {alts})")


def test_simplest_consensus20_again() -> None:
  refSeq = [ 'A', 'C', 'C', 'C' ]
  depths =  [ 100, 3, 0, 3]
  dumbRows = [VCFCall(AF=.01, DP=1, INDEL=False, ALT='X', POS=i+2, REF='C') for i in range(3)]
  alts = [VCFCall(AF=.81, DP=100, INDEL=False, ALT='G', POS=1, REF='A')] + dumbRows
  actual, log = consensus(depths, refSeq, alts, pa=PercentAnalysis(10, 80))

  # actual, log = zip(*actual) # [x[0] for x in actual]
  print(f"log: {log}")
  expected = [ 'G', 'N', 'N', 'N' ]
  assert  expected == list(actual) , f"{actual} != {expected}"



# def test_runner():
#   io_runner('testdata/Den1__WestPac__1997.fasta', 'testdata/fullsample.bam.vcf', 
#   'testdata/fullsample.bam', 'out_io_test.fasta', 'io_test.log', PercentAnalysis(10, 80, 30))
   



