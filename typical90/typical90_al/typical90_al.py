from math import lcm
a,b=map(int,input().split())
ans=lcm(a,b)
if ans>10**18:
    print("Large")
else:
    print(ans)