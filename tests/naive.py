import sys
from math import sqrt
n = int(sys.stdin.readline())
points = []

for i in range(n):
    line = sys.stdin.readline()
    args = line.split(' ')
    points.append((float(args[0]), float(args[1])))

def distance(a, b):
    ax, ay = a
    bx, by = b
    dx = ax - bx
    dy = ay - by
    return round(sqrt(dx * dx + dy * dy))

tour = [-1 for i in range(n)]
used = [False for i in range(n)]

tour[0] = 0
used[0] = True

for i in range(1, n):
    best = -1
    for j in range(n):
        tour_v = points[tour[i - 1]]
        j_v = points[j]
        best_v = points[best]
        if (not used[j]) and (best == -1 or distance(tour_v, j_v) < distance(tour_v, best_v)):
            best = j
    tour[i] = best
    used[best] = True

for t in tour:
    print(t)
