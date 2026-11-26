"""Verify one rational sweep certificate from saved interpolant nodes."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

from sweep_polynomial import verify


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',required=True,type=Path)
    parser.add_argument('--min-hat',required=True)
    args=parser.parse_args()
    raw=args.input.read_bytes()
    data=json.loads(gzip.decompress(raw))
    result=verify(data,Fraction(args.min_hat))
    result['input_sha256']=hashlib.sha256(raw).hexdigest()
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
