from vcf.model import _Record
from typing import Dict,Optional,Sequence,T, Callable,NamedTuple,Tuple
from toolz.itertoolz import first, second
from dataclasses import dataclass
import sys
DNA = 'ACTGactg'

def simplest(mind: int, majority: float, ref: str, depths: Sequence[int], vars: Sequence[_Record]):
    assert 0 < majority < 1
    withn_consensus = ref[:]
    under_depths, msgs = [], []
    DEFAULT_UNDER_DEPTH = 'N'
    for i, dp in enumerate(depths):
        if dp < mind:
            under_depths.append(i, dp)
            withn_consensus[i] = DEFAULT_UNDER_DEPTH
    def consider_var(var):
        return any(lambda af: af > (100-majority), var.AF)
    out_consensus = withn_consensus[:]
    filtered_vars = filter(consider_var, vars)
    DEFAULT_BASE = 'N'
    msgs.append(f"The following positions were replaced with `N` becauase they were under minimum depth {mind} according to the bam file: {under_depths}")
    for var in filtered_vars:
        pos = var.POS - 1
        if var.is_indel(): # includes deletion
            base, msg = DEFAULT_BASE, f"called POS={var.POS} the default base {DEFAULT_BASE} because it was called an indel by pyvcf: {var}."
        alt_base = try_call_alt(var, majority)
        if alt_base:
            base, msg = alt_base, f"called POS={var.POS} alternate because base {base} had AF >= {majority}"
        else:
            ambig_base = try_call_ambiguous(majority, withn_consensus[pos], var)
            if ambig_base:
                base, msg = ambig_base, f"called POS={var.POS} 'ambiguous' base {base} because there was no ALT with AF >= {majority} but at least one ALT was >= {1-majority}"
            else:
                base, msg = DEFAULT_BASE, f"called POS={var.POS} the default base {DEFAULT_BASE} because calling it as an alternate and calling it as an ambiguous base failed, and the depth was still over minimum {mind}"
        # base should only have length greater than 1 if it was an alternate call.
        out_consensus[pos:pos+len(base)] = base
        msgs.append(msg)
        # sys.stderr.write(f"At position {i}, evaluating variant {var}, attempted to make an ambiguity from cobined bases and failed! Calling `N`")

def try_call_alt(var, majority) -> Optional[str]:
   called = [ var.ALT[i] for i in range(var.ALT) 
                     if var.AF[i] > majority and var.ALT[i] in DNA]
   return None if not called else called

def try_call_ambiguous(majority, ref_base, var):
    alt_over_ambig_thresh = [ var.ALT[i] for i in range(var.ALT) 
                 if var.AF[i] > (1 - majority) ]
    only_first_alt_bases = [ x[0] for x in alt_over_ambig_thresh ]
    prc = 1-sum(var.AF) 
    ref_contrib = [ref_base] if prc >= 1-majority else []
    # we ignore the ref base if it is ambiguous, and we ignore the alt base if it is weird (like '*')
    # this can result in an "ambiguous call" giving a non-ambiguous base (i.e. any of A,T,C,G)
    valid_ambig_combo = filter(DNA.__contains__, only_first_alt_bases + ref_contrib)
    return make_ambig(valid_ambig_combo)

def make_ambig(bs_: Sequence[str]) -> Optional[str]:
    bs = list(map(str.upper, bs_))
    AMBIGUITY_TABLE = { 'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B', 'ACGT': 'N' }
    assert all(AMBIGUITY_TABLE.keys().__contains__, bs)
    cleaned = ''.join(sorted(bs))
    return AMBIGUITY_TABLE.get(cleaned)

'''def main(mind: int, majority: float, ref: str, 
         depths: Sequence[int], varMap: Dict[int,_Record]):
    assert 0 < majority < 1
    consensus, i = [], 0
    while i < len(ref):
        var: Optional[_Record] = varMap.get(i)
        if not var:
            if depths[i] < mind:
                bs=['N']
            else:
                bs = [ref[i]]
        else:
            must_be_list = lambda xs: xs if hasattr(xs, '__iter__') else [xs]
            # DO I need to account for base_caller.py? Considering that tool has already done its call?
            DumbAlt = NamedTuple('DumbAlt', [('alt', str), ('af', float)])
            aaf_alt_pairs = map(DumbAlt, must_be_list(var.ALT), must_be_list(var.AF))
            candidates = filter(lambda x: x.af > (100-majority), aaf_alt_pairs)
            if not candidates:
                bs = [ref[i]]
            else:
                called_alt = find(lambda x: x.af >= majority, candidates)
                if called_alt:
                    bs = called_alt.alt
                else:
                    get_contribution = lambda s: s[ : var.affected_end - var.affected_start]
                    candidate_seqs = [v.alt for v in candidates]
                    alt_contributions = list(map(get_contribution, candidate_seqs))
                    prc = 1-sum(var.aaf) # TODO: a given ALT could be len(1) but I think that's okay
                    if prc >= 1-majority:                   # ref_contribution = ref[i] if prc >= 1-majority else ''
                        ref_contributions = ref[var.affected_start:var.affected_end]
                    else:
                        ref_contributions = []
                    results = list(map(make_ambig, ref_contributions + alt_contributions))
                    bad_results = filter(second, results)
                    if bad_results:
                         sys.stderr.write(f"At position {i}, evaluating variant {var}, \
                             attempted to make an ambiguity from cobined bases and failed! Calling `N` \n")
                         bs = 'N'
        i += len(bs)
        consensus.extend(bs)
    return consensus
def find(pred: Callable, xs: Sequence[T]) -> Optional[T]:
    return next(filter(pred, xs)) # camelCase
def make_ambig(bs: Sequence[str]) -> Tuple[str,Optional[str]]:  # Either-like error-handling with tuple
    AMBIGUITY_TABLE = { 'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B', 'ACGT': 'N' }
    assert all(AMBIGUITY_TABLE.keys().__contains__, bs)
    DE_AMBIG_TABLE = dict(map(reversed, AMBIGUITY_TABLE.items()))
    de_ambiged = map(DE_AMBIG_TABLE.get, set(bs.upper()))
    cleaned = ''.join(sorted(de_ambiged))
    AMBIGUITY_TABLE[cleaned]
    '''