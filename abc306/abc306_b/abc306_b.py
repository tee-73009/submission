a = list(map(int, input().split()))
e = 0
for b in range(64):
  c = a[-b-1]
  if b == 0:
    d = 0
  else:
    d = a[-b]
  e =  2 * e + c
print(e)