#!/usr/bin/env python3
"""Reuse the independently gradient-checked GPU engine with the current 19D/21D baselines."""
import argparse,sys
from pathlib import Path
from baselines import baseline,DATA

sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'compute17_20_20261104'))
import gpu_deform
gpu_deform.DATA=DATA;gpu_deform.baseline=baseline
p=argparse.ArgumentParser();p.add_argument('--dimension',type=int,choices=(19,21),required=True)
p.add_argument('--seconds',type=float,default=180);p.add_argument('--seed',type=int,default=1)
p.add_argument('--device',type=int,default=0);a=p.parse_args()
gpu_deform.run(a.dimension,a.seconds,a.seed,a.device)
