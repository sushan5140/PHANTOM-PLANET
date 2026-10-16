import argparse,csv,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).parent/'src'))
from phantom.sim import CASES,write_csv
from phantom.baseline import periodic_candidate

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--seeds',type=int,default=10)
    p.add_argument('--n',type=int,default=90)
    p.add_argument('--output',default='results')
    args=p.parse_args()
    if args.seeds<1: p.error('seeds must be positive')
    out=Path(args.output); (out/'synthetic').mkdir(parents=True,exist_ok=True)
    rows=[]
    for case in CASES:
        for seed in range(args.seeds):
            path=out/'synthetic'/f'{case}_{seed:03d}.csv'
            write_csv(path,case,seed,args.n)
            data=np.genfromtxt(path,delimiter=',',names=True)
            candidate=periodic_candidate(data['time_days'],data['rv_ms'])
            rows.append({'case':case,'seed':seed,**candidate,'path':str(path)})
    with (out/'periodogram_candidates.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print(json.dumps({'cases':len(rows),'output':str(out),'note':'Candidates only; no validated planet classifications and no juliet fits'}))
if __name__=='__main__':main()
