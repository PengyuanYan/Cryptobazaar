# Cryptobazaar

This repository provides a re-implementation of Cryptobazaar Crypotbazaar with ICICLE library for efficiency. The benchmarks are same as the origianl Cryptobazaar. 

Please use "cargo test -- --nocapture --test-threads=1" to run the tests.
If the code is tested in multi-thread mod the code will crash becasue the conflict between ICICLE domains.

The microbenchmarks can be re-run as in the original code.

In the default setting, the code cannot run with the GPU.
If you want it run with GPU, just remove "//" before the first line of code of the function load_backend() in the utils.rs.

In the end of the README file, an instruction of how to test the code on Google Colab is provided.

## Rust

To setup Rust, please follow the [official installation instructions](https://www.rust-lang.org/tools/install).

## Benchmarks

There are six microbenchmarks, namely for the computation of the four validity proofs, the AV protocol, and the results vector:

- **Benchmark 1: Validity proof $\pi_{x_i}$ (Table 1a)**
    ```
    cargo bench --bench veceq
    ```

- **Benchmark 2: Validity proof $\pi_{r_i}$ (Table 1a)**
    ```
    cargo bench --bench nonzero
    ```

- **Benchmark 3: Validity proof $\pi_{b_i}$ (Table 1a)**
    ```
    cargo bench --bench lderivative
    ```

- **Benchmark 4: Validity proof $\pi_{Z_i}$ (Table 1a)**
    ```
    cargo bench --bench ipa
    ```

- **Benchmark 5: AV matrix $Y$ (Table 1b)**
    ```
    cargo bench --bench auctioneer_r1
    ```

- **Benchmark 6: Results vector $R$ (Table 1c)**
    ```
    cargo bench --bench auctioneer_r2
    ```

## Colab Instruction
1. Click this link "https://colab.research.google.com"
2. Close the pop up window.
3. Click the triangle on the top right corner and choose the second option.
4. Choose a T4 GPU and save.
5. Click connect.
6. Click the terminal at the bottom left corner.
7. Copy and paste the following command into the terminal "curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh"
8. Input 1 into the terminal.
9. Input exit.
10. Close the terminal and open it again.
11. Copy and paste "wget https://github.com/ingonyama-zk/icicle/releases/download/v3.9.2/icicle_3_9_2-ubuntu20-cuda122.tar.gz" into terminal.
12. Copy and paste "tar -xvf icicle_3_9_2-ubuntu20-cuda122.tar.gz -C /opt" into terminal.
13. Copy and paste "git clone https://github.com/PengyuanYan/Cryptobazaar.git" into terminal.
14. Copy and paste "cd Cryptobazaar".
15. Run any benchmark you want.



