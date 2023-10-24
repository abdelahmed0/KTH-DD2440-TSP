import subprocess
from subprocess import Popen, PIPE, STDOUT
import sys
import os
from math import sqrt

PYTHON = "python3"
SHELL = "bash"
NAIVE_PATH = "naive.py"

if len(sys.argv) <= 2:
    print(f"Usage: {PYTHON} {sys.argv[0]} FILE TESTS [IGNORE_NAIVE=0]")
    exit()

executable = sys.argv[1]
test_dir = sys.argv[2]

ignore_naive = False
if len(sys.argv) >= 4:
    ignore_naive = int(sys.argv[3]) != 0

def distance(a, b):
    ax, ay = a
    bx, by = b
    dx = ax - bx
    dy = ay - by
    return round(sqrt(dx * dx + dy * dy))

def compute_length(points, tour):
    length = 0
    last_point = points[tour[0]]
    for i in range(1, len(tour)):
        next_point = points[tour[i]]
        length += distance(last_point, next_point)
        last_point = next_point
    return length

def tour_from_out(outs):
    tour = []
    for l in outs.decode().split('\n'):
        if len(l) > 0:
            tour.append(int(l))
    return tour

def length_from_out(points, outs):
    tour = tour_from_out(outs)
    return compute_length(points, tour)

tests = os.listdir(test_dir)
for test in tests:
    test_path = os.path.join(test_dir, test)
    case = open(test_path, 'r').read()
    points = []
    case_split = case.split('\n')
    n = int(case_split[0])
    for i in range(1, n + 1):
        args = case_split[i].split(' ')
        points.append((float(args[0]), float(args[1])))

    p = Popen([executable], stdout=PIPE, stdin=PIPE, stderr=subprocess.STDOUT)
    outs, errs = p.communicate(input=case.encode(), timeout=2)

    print(f'{test}: {length_from_out(points, outs)}', end='')
    
    if not ignore_naive:
        naive = Popen([PYTHON, NAIVE_PATH], stdout=PIPE, stdin=PIPE, stderr=subprocess.STDOUT)
        outs, errs = naive.communicate(input=case.encode(), timeout=2)
        print(f' ({length_from_out(points, outs)})', end='')
    
    print('')
