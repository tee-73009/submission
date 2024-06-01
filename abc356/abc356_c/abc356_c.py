def solve(N, M, K, tests):
    from itertools import combinations

    valid_count = 0

    for mask in range(1 << N):
        valid = True
        for test in tests:
            C, keys, result = test
            count = sum((mask >> (key - 1)) & 1 for key in keys)
            if (result == 'o' and count < K) or (result == 'x' and count >= K):
                valid = False
                break
        if valid:
            valid_count += 1

    return valid_count

import sys
input = sys.stdin.read
data = input().split()

N = int(data[0])
M = int(data[1])
K = int(data[2])

tests = []
index = 3
for _ in range(M):
    C = int(data[index])
    keys = list(map(int, data[index+1:index+1+C]))
    result = data[index+1+C]
    tests.append((C, keys, result))
    index += C + 2

print(solve(N, M, K, tests))
