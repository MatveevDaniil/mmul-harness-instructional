import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def cblas_mflop(n: int) -> int:
  return 2 * n ** 3

def basic_mflop(n: int) -> int:
  return 2 * n ** 3

def blocked_mflop(n: int, block_size: int) -> int:
  return 2 * n ** 3 + n ** 3 // block_size

def strassen_mflop(n: int, dim_limit: int) -> int:
  if n == 1:
    return 1
  elif n <= dim_limit:
    return 2 * n ** 3
  else:
    return 7 * strassen_mflop(n // 2, dim_limit) + (15 * n ** 2) // 4

for plot_mode in 'mflops', 'time':
  ### blas
  blas_benchmark = 'blas'
  ylabel = 'MFLOPs' if plot_mode == 'mflops' else 'Time (s)'

  def plot_blas():
    df = pd.read_csv(f"results/blas.csv").drop(0)
    runtime = df['Time']
    mflops = df['n'].apply(cblas_mflop) / df['Time'] / 10 ** 6
    y_value = mflops if plot_mode == 'mflops' else runtime
    plt.scatter(df['n'], y_value, label=blas_benchmark, s=10, marker='^', color='red')
    plt.plot(df['n'], y_value, color='red')

  def set_axes_and_save(fname: str, title: str):
    plt.xlabel("Matrix dimension")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.savefig(f'visualization_scripts/{fname}.png', dpi=300)
    plt.clf()

  ### basic
  basic_benchmarks = [f'basic-{benchmark}' for benchmark in ('ijk', 'ikj', 'jik', 'jki', 'kij', 'kji')]
  for benchmark in basic_benchmarks:
    df = pd.read_csv(f"results/{benchmark}.csv").drop(0)
    runtime = df['Time']
    mflops = df['n'].apply(cblas_mflop) / df['Time'] / 10 ** 6
    y_value = mflops if plot_mode == 'mflops' else runtime
    plt.scatter(df['n'], y_value, label=benchmark, s=10)
    plt.plot(df['n'], y_value)
  plot_blas()
  set_axes_and_save(f'basic_{plot_mode}', "Three-loop algorithms' performance")

  ### blocked
  blocked_benchmarks = ['blocked', 'blocked-templated']
  markers = {'blocked': 'o', 'blocked-templated': 's'}
  for benchmark in blocked_benchmarks:
    for block_size in 2, 16, 32, 64:
      df = pd.read_csv(f"results/{benchmark}.csv").drop([0, 1, 2, 3])
      df = df[df['Block Size'] == block_size] 
      runtime = df['Time']
      mflops = df['n'].apply(lambda n: blocked_mflop(n, block_size)) / df['Time'] / 10 ** 6
      y_value = mflops if plot_mode == 'mflops' else runtime
      plt.scatter(df['n'], y_value, label=f'{benchmark}-{block_size}', s=10, marker=markers[benchmark])
      plt.plot(df['n'], y_value)
  plot_blas()
  set_axes_and_save(f'blocked_{plot_mode}', "Fixed-block algorithms' performance")

  ### strassen
  strassen_benchmarks = ['strassen-1', 'strassen-16', 'strassen-64', 'strassen-256']
  for benchmark in strassen_benchmarks:
    dim_limit = int(benchmark.split('-')[1])
    df = pd.read_csv(f"results/{benchmark}.csv").drop(0)
    runtime = df['Time']
    mflops = df['n'].apply(lambda n: strassen_mflop(n, dim_limit)) / df['Time'] / 10 ** 6
    y_value = mflops if plot_mode == 'mflops' else runtime
    plt.scatter(df['n'], y_value, label=benchmark, s=10)
    plt.plot(df['n'], y_value)
  plot_blas()
  set_axes_and_save(f'strassen_{plot_mode}', "Strassen-Winograd algorithms' performance") 

  ### compare winners
  if plot_mode == 'time':
    ijk_time = pd.read_csv(f"results/basic-ijk.csv").drop(0)['Time']
    blas_time = pd.read_csv(f"results/blas.csv").drop(0)['Time']

    # ikj
    df = pd.read_csv(f"results/basic-ikj.csv").drop(0)
    plt.scatter(df['n'], df['Time'], label='basic-ikj', s=10)
    ijk_ratio = np.array(ijk_time) / np.array(df['Time'])
    print(ijk_ratio.mean(), ijk_ratio.max(), ijk_ratio[-1])
    plt.plot(df['n'], df['Time'])

    # blocked-templated-32
    df = pd.read_csv(f"results/blocked-templated.csv").drop([0, 1, 2, 3])
    df = df[df['Block Size'] == 32]
    plt.scatter(df['n'], df['Time'], label='blocked-templated-32', s=10)
    plt.plot(df['n'], df['Time'])
    ijk_ratio = np.array(ijk_time) / np.array(df['Time'])
    print(ijk_ratio.mean(), ijk_ratio.max(), ijk_ratio[-1])

    # strassen-64
    df = pd.read_csv(f"results/strassen-64.csv").drop(0)
    plt.scatter(df['n'], df['Time'], label='strassen-64', s=10)
    plt.plot(df['n'], df['Time'])
    set_axes_and_save(f'best_run_time', "Run time of best methods among from each category") 
    ijk_ratio = np.array(ijk_time) / np.array(df['Time'])
    print(ijk_ratio.mean(), ijk_ratio.max(), ijk_ratio[-1])
    blas_ratio = np.array(blas_time) / np.array(df['Time'])
    print(blas_ratio.mean(), blas_ratio.max(), blas_ratio[-1])