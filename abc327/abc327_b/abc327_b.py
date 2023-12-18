b = int(input())
for a in range(1,10000):
  if a ** a == b:
    print(a)
    exit()
  if a ** a > b:
    print(-1)
    exit()