import sys
import os
import random

SCALE = 10**6

if len(sys.argv) <= 1:
    print(f"Usage: python3 {sys.argv[0]} DIRECTORY NUM_OF_CASES [MAX_N=50] [SEED]")
    exit()

dir = sys.argv[1]
cases = int(sys.argv[2])
max_n = 50
seed = 1243
if len(sys.argv) >= 4:
    max_n = int(sys.argv[3])

if len(sys.argv) >= 5:
    seed = int(sys.argv[4])

random.seed(seed)

if not os.path.exists(dir):
    os.makedirs(dir)

for j in range(cases):
    path = os.path.join(dir, f"{j}")
    with open(path, 'w') as f:
        n = random.randint(3, max_n)
        f.write(f'{n}\n')
        for i in range(n):
            x = SCALE * random.random()
            y = SCALE * random.random()
            f.write(f'{x} {y}\n')
