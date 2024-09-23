import matplotlib.pyplot as plt


# complexity of my implementation of strassen-winograd algorithm
def T(n: int, dim_limit: int) -> int:
  # make sure n is a power of 2
  assert bin(n).count('1') == 1
  if n == 1:
    return 1
  elif n <= dim_limit:
    return 2 * n ** 3
  else:
    return 7 * T(n // 2, dim_limit) + (15 * n ** 2) // 4
  
n_values = (64, 128, 256, 512, 1024, 2048)
dim_limits = (1, 16, 64, 256)

mega = 10 ** 9

for dim_limit in dim_limits:
  plt.plot(n_values, [T(n, dim_limit) / mega for n in n_values], label=f"dim_limit={dim_limit}")
plt.plot(n_values, [2 * n ** 3 / mega for n in n_values], label="naive")

plt.xlabel("Matrix dimension")
plt.ylabel("MFLOPs")
plt.title("Strassen-Winograd algorithm complexity")
plt.legend()
# plt.show()
plt.savefig("./visualization_scripts/strassen_comp.png", dpi=300)