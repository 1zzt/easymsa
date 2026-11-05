"""Independent small-input tests for the selected derivative core."""
import sys, random, itertools, subprocess

binary, mode = sys.argv[1], int(sys.argv[2])
rng=random.Random(20261105)

def align(sequences):
    source=''.join(f'>s{i}\n{s}\n' for i,s in enumerate(sequences))
    p=subprocess.run([binary,'-i','-','-o','-','--quiet'],input=source,
                     capture_output=True,text=True,timeout=10,check=True)
    rows=[]
    for line in p.stdout.splitlines():
        if line.startswith('>'): rows.append('')
        else: rows[-1]+=line
    assert len(rows)==len(sequences)
    assert [r.replace('-','') for r in rows]==sequences
    assert len({len(r) for r in rows})==1
    assert all(any(c!='-' for c in col) for col in zip(*rows))
    return rows

for _ in range(60):
    sequences=[''.join(rng.choices('ACGTUN',k=rng.randint(1,30))) for _ in range(rng.randint(1,8))]
    rows=align(sequences)
    assert rows==align(sequences)
    if mode==2: assert rows==list(reversed(align(list(reversed(sequences)))))

def substitution(a,b):
    if 'N' in (a,b): return 0
    return 2 if a.replace('U','T')==b.replace('U','T') else -1

def oracle(a,b):
    def walk(i,j,state,score):
        if i==len(a) and j==len(b): return score
        alternatives=[]
        if i<len(a) and j<len(b): alternatives.append(walk(i+1,j+1,0,score+substitution(a[i],b[j])))
        if i<len(a): alternatives.append(walk(i+1,j,1,score-(1 if state==1 else 3)))
        if j<len(b): alternatives.append(walk(i,j+1,2,score-(1 if state==2 else 3)))
        return max(alternatives)
    return walk(0,0,0,0)

def score(rows):
    value=0;state=0
    for a,b in zip(*rows):
        if '-' in (a,b):
            next_state=2 if a=='-' else 1
            value-=1 if next_state==state else 3;state=next_state
        else: value+=substitution(a,b);state=0
    return value

if mode==3:
    words=[''.join(x) for n in range(1,4) for x in itertools.product('AC',repeat=n)]+['N','TU','AN']
    for a,b in itertools.product(words,repeat=2):
        assert score(align([a,b]))==oracle(a,b)
print(f'PASS mode={mode}: independent algorithm tests')
