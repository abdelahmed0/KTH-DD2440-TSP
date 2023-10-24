from math import sqrt

def distance(a, b):
    ax, ay = a
    bx, by = b
    dx = ax - bx
    dy = ay - by
    return round(sqrt(dx * dx + dy * dy))

def find_min(ls):
    result = ls[0]
    for x in ls:
        if x < result:
            result = x
    return result

def find_max(ls):
    result = ls[0]
    for x in ls:
        if x > result:
            result = x
    return result

def mean(ls):
    n = len(ls)
    sum = 0
    for x in ls:
        sum += x
    return sum / n
