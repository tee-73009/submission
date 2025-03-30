import sys
sys.set_int_max_str_digits(0)
n=int(input())
if n==1:
  print(1)
else:
  ans=(1+n)*pow(2,n-2,998244353)
  print(ans%998244353)