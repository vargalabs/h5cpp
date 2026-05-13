---
hide:
  - toc
---

## Motivating Example: Parallel HDF5 for Strategy Research

A common quantitative workflow is to evaluate many independent strategies over the same historical dataset. For example:

* multiple mean-reversion variants
* different execution models
* alternative feature sets
* parameter sweeps for thresholds, lookback windows, and risk constraints

This workload is naturally **embarrassingly parallel**. Each worker can process its own strategy or parameter set with minimal communication: $R_{i,j} = \mathcal{S}_i(D,\theta_j)$ where:

* $D$ is the shared historical dataset
* $\mathcal{S}_i$ is the $i$-th strategy family
* $\theta_j$ is the $j$-th parameter vector
* $R_{i,j}$ is the resulting performance record, signal matrix, or trade log

The important point is that the compute is independent, but the input and output formats should remain consistent and reusable.

### Why parallel HDF5 helps

In this setup, HDF5 acts as a shared research container:

* **input datasets** hold historical prices, returns, features, or market event streams
* **workers** read the same common data or pre-assigned slices
* **output datasets** collect PnL curves, summary statistics, factor exposures, or structured trade logs

This gives quantitative teams a clean path from:

* local prototype
* to workstation batch runs
* to cluster-scale parameter sweeps

without changing the storage model.

### Typical research layout

A single HDF5 container might look like this:

```text
/research/input/prices
/research/input/returns
/research/input/features
/research/input/events

/research/output/strategy_0001/pnl
/research/output/strategy_0001/stats
/research/output/strategy_0001/trades

/research/output/strategy_0002/pnl
/research/output/strategy_0002/stats
/research/output/strategy_0002/trades
```

Each MPI rank can be assigned one strategy, one parameter block, or one date range.

### Why this is embarrassingly parallel

This pattern is attractive because it avoids the hard distributed-systems problems. Each worker does roughly:

1. open the shared HDF5 research container
2. read required inputs
3. run one strategy or parameter block
4. write results to its own output path

There is little or no synchronization beyond startup and final completion. That makes the workload ideal for:

* parameter sweeps
* walk-forward studies
* feature ablations
* cross-sectional factor tests
* Monte Carlo scenario batches
* multi-strategy backtest farms

### Example sketch

The following simplified sketch shows the idea: each MPI rank reads the same input return matrix, runs a different strategy variant, and writes its own result vector.

```cpp
#include <mpi.h>
#include <armadillo>
#include <h5cpp/all>

struct stats_t {
    double sharpe;
    double pnl;
    double max_drawdown;
};

auto run_strategy(const arma::mat& returns, int rank) -> std::pair<arma::vec, stats_t> {
    arma::vec pnl = arma::sum(returns, 1) * (0.1 * (rank + 1));

    stats_t stats{
        .sharpe = arma::mean(pnl) / std::max(1e-12, arma::stddev(pnl)),
        .pnl = arma::accu(pnl),
        .max_drawdown = 0.0
    };
    return {pnl, stats};
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    auto fd = h5::open(
        "research.h5",
        H5F_ACC_RDWR,
        h5::mpiio({MPI_COMM_WORLD, MPI_INFO_NULL})
    );

    arma::mat returns = h5::read<arma::mat>(fd, "/research/input/returns");

    auto [pnl, stats] = run_strategy(returns, rank);

    std::string base = "/research/output/strategy_" + std::to_string(rank);

    h5::write(fd, base + "/pnl", pnl, h5::independent);
    h5::write(fd, base + "/stats", stats, h5::independent);

    MPI_Finalize();
}
```

### What this buys you

From a quant research perspective, this gives you:

* a single structured container for both inputs and outputs
* reproducible storage for large batch experiments
* direct interoperability with Python, Julia, MATLAB, and R
* easy post-processing of all results after the run
* a clean separation between **strategy logic** and **data plumbing**

In other words, the cluster is used for what matters — evaluating ideas — instead of wasting effort on ad hoc file conventions and result collection scripts.

### Practical payoff

For trading and research teams, parallel HDF5 is valuable when:

* the dataset is too large to duplicate carelessly
* many independent experiments must run over the same inputs
* output needs to remain structured and queryable afterward
* the same data should feed both backtests and later analysis

This is especially natural for quants with backgrounds in physics, engineering, optimization, or scientific computing, where HDF5 is already familiar territory and matrix-oriented workflows are standard. The storage model stays stable while the compute scales out.

