# Branch Predictor Simulator in C++

This repository contains a C++ simulation of various branch prediction algorithms used in CPU architectures. Branch predictors are essential for enhancing CPU pipeline efficiency, minimizing stalls, and maximizing instruction throughput by predicting the outcome of conditional branches.


## Overview
Branch prediction is critical for reducing the performance impact of conditional branch instructions. This project simulates various branch prediction strategies, providing insights into their accuracy and efficiency.

## Features
- Implements multiple branch prediction algorithms, including static and dynamic predictors.
- Supports evaluation metrics such as hit/miss ratio and prediction accuracy.
- Extensible design to easily add or modify prediction algorithms.

## Branch Prediction Algorithms
The simulator includes the following branch prediction algorithms:
1. **Static Predictor**: Always predicts either "taken" or "not taken."
2. **One-Bit Predictor**: Maintains a single-bit history for each branch instruction.
3. **Two-Bit Predictor**: Uses a two-bit saturating counter to improve prediction accuracy by tracking recent branch behavior.
4. **GShare Predictor**: Combines global branch history with the program counter (PC) to make predictions, reducing aliasing and increasing accuracy.

## File Structure
- **bp_api.h**: Contains the branch predictor API and classes for different prediction algorithms.
- **bp_main.cpp**: Main entry point for the simulator, handling input, initializing predictors, and calculating performance metrics.
- **test1.in**: Sample input file that provides a sequence of branch instructions and outcomes.
- **test1.out**: Expected output file for comparison with the actual results from the simulation.

## Sample Input Format
The input format for this simulator consists of branch instructions, each with the following fields:
- **Program Counter (PC)**: The address of the branch instruction.
- **Outcome**: Indicates whether the branch was actually taken (`T`) or not taken (`N`).
- **Target Address**: The address to which the branch would jump if taken.

### Example Input
Each line in the input file corresponds to a branch instruction, formatted as follows:

0x3a44 N 0x3a48 
0x5c44 T 0x36f30 
0x2244 N 0x2248 
0x1444 T 0x1fc4c


In this example:
- `0x3a44 N 0x3a48` means that the branch at address `0x3a44` was **not taken** and would have jumped to `0x3a48` if taken.
- `0x5c44 T 0x36f30` means that the branch at address `0x5c44` was **taken** and jumped to `0x36f30`.

## Usage
1. **Clone the repository**:
    ```bash
    git clone https://github.com/foziD123/Branch-Predictor-CPP.git
    cd Branch-Predictor-CPP
    ```

2. **Compile the code**:
    ```bash
    g++ -o branch_predictor bp_main.cpp
    ```

3. **Run the simulation**:
    ```bash
    ./branch_predictor < input_file > output_file
    ```
   Replace `input_file` with a test file to evaluate the prediction accuracy.




