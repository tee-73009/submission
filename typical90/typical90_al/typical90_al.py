from math import lcm

a, b = map(int, input().split())
c = lcm(a, b) 

print(c if c <= 10**18 else "Large")