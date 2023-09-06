N = int(input())
ans = 0
for a in range(1, 38):
    for b in range(1, 26):
        if 3 ** a + 5 ** b == N :
            print(a, b)
            ans = 1
if ans == 0:
  print(-1)